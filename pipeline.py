#!/usr/bin/env python3
"""DeepDNA drug-gene pipeline.

Pulls genes, drug-gene pairs, variants, and their sources from the DeepDNA API
into a SQLite database (deepdna.db), one timestamped snapshot per run.

    python pipeline.py run                      # pull everything and run quality checks
    python pipeline.py run --genes CYP2D6,DPYD  # only these genes
    python pipeline.py run --skip-variants      # drug-gene pairs only
    python pipeline.py run --strict             # exit with code 2 if a quality check warns
    python pipeline.py summary                  # show the run history
    python pipeline.py inspect genes            # print a raw API response
    python pipeline.py inspect variant rs4244285
    python pipeline.py reprocess                # rebuild tables from stored raw data (no API calls)
    python pipeline.py prune --keep 52 --keep-raw 2   # keep 52 runs, raw API data for the newest 2

Exit codes: 0 = ok, 1 = the run failed, 2 = finished but --strict and a check warned.
The functions can also be imported from a notebook (import pipeline).

Educational information only, not medical advice. Check deepdna.ai/terms before
scheduling frequent calls.
"""
import argparse
import json
import logging
import os
import re
import sqlite3
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote, unquote

import requests

log = logging.getLogger("pipeline")

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
BASE = "https://deepdna.ai/api/v1"
DEFAULT_DB = os.environ.get("DEEPDNA_DB", str(Path(__file__).parent / "deepdna.db"))
DELAY_SECONDS = 1.0  # pause between calls; the API is free and shared
USER_AGENT = "deepdna-pipeline/1.0 (learning project)"

# Possible field names for a gene symbol inside the GET /genes response.
GENE_SYMBOL_KEYS = ("symbol", "gene", "id", "name")

# Field names inside GET /pgx/{gene}.
FIELD_MAP = {
    "drug_list_key": "pairs",      # list of drug-gene pairs
    "drug_name_key": "drug",
    "guideline_key": "guideline",  # e.g. "CPIC" or "DPWG; FDA label"
    "level": "level",              # e.g. "CPIC A"
    "summary": "summary",
    "phenotypes": "phenotypes",    # list of strings (can be null)
    "sources": "sources",          # list of citations (can be null)
}

# Field names inside GET /variants/{rsID}.
VARIANT_MAP = {
    "name": "name",
    "alleles": "alleles",
    "effect": "effect",
    "frequency": "frequency_note",
    "sources": "sources",
}

TABLES = ("raw_responses", "drug_gene_pairs", "phenotypes", "pair_sources",
          "variants", "variant_sources", "run_summary")

PHENOTYPE_PATTERN = re.compile(r"^\s*(?P<name>[^(:]+?)\s*(?:\((?P<code>[^)]*)\))?\s*:\s*(?P<desc>.*)$")
RSID_PATTERN = re.compile(r"\brs\d+\b", re.IGNORECASE)
SOURCE_LABEL_KEYS = ("label", "name", "title", "source", "citation", "id")
SOURCE_URL_KEYS = ("url", "link", "href")

_session = requests.Session()
_session.headers["User-Agent"] = USER_AGENT


# --------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------
def now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def as_text(value):
    """SQLite can't store lists or dicts, so turn them into JSON text."""
    if value is None or isinstance(value, (str, int, float)):
        return value
    return json.dumps(value)


def as_list(value):
    """Always return a list: None -> [], a single value -> [value]."""
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


def parse_source(item):
    """Turn one source (a string or an object) into (label, url)."""
    if isinstance(item, str):
        return item, (item if item.startswith("http") else None)
    if isinstance(item, dict):
        label = next((str(item[k]) for k in SOURCE_LABEL_KEYS if item.get(k)), None)
        url = next((str(item[k]) for k in SOURCE_URL_KEYS if item.get(k)), None)
        return (label or url or json.dumps(item)), url
    return str(item), None


def parse_phenotype(text):
    """Split 'Poor metabolizer (PM): two no-function alleles' into (name, code, description)."""
    if not isinstance(text, str):
        return None, None, None
    m = PHENOTYPE_PATTERN.match(text)
    if not m:
        return text.strip(), None, None
    return m["name"].strip(), (m["code"] or "").strip() or None, m["desc"].strip() or None


def find_rsids(body):
    """Find every rsID anywhere in a response, whatever its structure."""
    return sorted(set(m.lower() for m in RSID_PATTERN.findall(json.dumps(body or {}))))


def parse_gene_list(body):
    """Pull gene symbols out of a /genes response, whether items are strings or objects."""
    if not body:
        return []
    if isinstance(body, list):
        items = body
    else:  # use the first list-valued field (skips the "meta" block)
        items = next((v for v in body.values() if isinstance(v, list)), [])
    symbols = []
    for item in items:
        if isinstance(item, str):
            symbols.append(item)
        elif isinstance(item, dict):
            for key in GENE_SYMBOL_KEYS:
                if item.get(key):
                    symbols.append(str(item[key]))
                    break
    return sorted(set(symbols))


# --------------------------------------------------------------------------
# Database
# --------------------------------------------------------------------------
def init_db(con):
    con.executescript("""
    CREATE TABLE IF NOT EXISTS raw_responses (
        id INTEGER PRIMARY KEY,
        run_ts TEXT, endpoint TEXT, status INTEGER, body TEXT
    );
    CREATE TABLE IF NOT EXISTS drug_gene_pairs (
        run_ts TEXT, gene TEXT, drug TEXT, guideline TEXT,
        level TEXT, summary TEXT
    );
    CREATE TABLE IF NOT EXISTS phenotypes (
        run_ts TEXT, gene TEXT, drug TEXT, phenotype TEXT,
        phenotype_name TEXT, code TEXT, description TEXT
    );
    CREATE TABLE IF NOT EXISTS pair_sources (
        run_ts TEXT, gene TEXT, drug TEXT, source TEXT, url TEXT
    );
    CREATE TABLE IF NOT EXISTS variants (
        run_ts TEXT, rsid TEXT, gene TEXT, name TEXT, alleles TEXT,
        effect TEXT, frequency TEXT
    );
    CREATE TABLE IF NOT EXISTS variant_sources (
        run_ts TEXT, rsid TEXT, gene TEXT, source TEXT, url TEXT
    );
    CREATE TABLE IF NOT EXISTS run_summary (
        run_ts TEXT, genes_found INTEGER, genes_missing INTEGER, pairs INTEGER,
        genes_listed INTEGER, variants INTEGER
    );
    """)
    # Older databases were created without some columns; add them if missing.
    for table, col, kind in (("drug_gene_pairs", "level", "TEXT"),
                             ("drug_gene_pairs", "summary", "TEXT"),
                             ("run_summary", "genes_listed", "INTEGER"),
                             ("run_summary", "variants", "INTEGER"),
                             ("phenotypes", "phenotype_name", "TEXT"),
                             ("phenotypes", "code", "TEXT"),
                             ("phenotypes", "description", "TEXT")):
        try:
            con.execute(f"ALTER TABLE {table} ADD COLUMN {col} {kind}")
        except sqlite3.OperationalError:
            pass  # column already exists
    con.commit()


def store_raw(con, run_ts, endpoint, status, body):
    con.execute(
        "INSERT INTO raw_responses (run_ts, endpoint, status, body) VALUES (?,?,?,?)",
        (run_ts, endpoint, status, json.dumps(body)),
    )


# --------------------------------------------------------------------------
# Extract: call the API and keep the raw responses
# --------------------------------------------------------------------------
def fetch(endpoint):
    """GET one endpoint with simple retries. Returns (status, parsed_json_or_None)."""
    url = f"{BASE}/{endpoint}"
    for attempt in range(3):
        try:
            r = _session.get(url, timeout=20)
            if r.status_code in (200, 404):
                return r.status_code, r.json()
            log.warning("%s: HTTP %s (attempt %d)", endpoint, r.status_code, attempt + 1)
        except (requests.RequestException, ValueError) as exc:
            log.warning("%s: %s (attempt %d)", endpoint, exc, attempt + 1)
        time.sleep(2 ** attempt)
    return 0, None


def get_genes(con, run_ts, override=None):
    """Return the genes to pull: the override list, or everything GET /genes lists."""
    if override:
        return list(override)
    status, body = fetch("genes")
    store_raw(con, run_ts, "genes", status, body)
    genes = parse_gene_list(body) if status == 200 else []
    if not genes:
        con.commit()
        raise RuntimeError("Could not read a gene list from GET /genes. "
                           "Run `python pipeline.py inspect genes` to see the response, "
                           "or pass --genes.")
    if isinstance(body, dict) and body.get("count") not in (None, len(genes)):
        log.warning("API says count=%s but parsed %d genes (the list may be paginated)",
                    body.get("count"), len(genes))
    log.info("API lists %d genes", len(genes))
    return genes


def extract(con, run_ts, genes_override=None):
    """Store the gene list, then call /pgx/{gene} for every gene. Returns the list of genes."""
    genes = get_genes(con, run_ts, genes_override)
    for gene in genes:
        endpoint = f"pgx/{quote(gene)}"
        status, body = fetch(endpoint)
        store_raw(con, run_ts, endpoint, status, body)
        log.info("%s: %s", endpoint, status)
        time.sleep(DELAY_SECONDS)
    con.commit()
    return genes


def extract_variants(con, run_ts, genes):
    """Fetch GET /genes/{gene} for each gene, then GET /variants/{rsID} for every rsID found.
    Only stores raw responses; transform_variants() builds the tables from them."""
    rsids = set()
    for gene in genes:
        endpoint = f"genes/{quote(gene)}"
        status, body = fetch(endpoint)
        store_raw(con, run_ts, endpoint, status, body)
        time.sleep(DELAY_SECONDS)
        if status == 200:
            rsids.update(find_rsids(body))
    log.info("Found %d unique rsIDs across %d genes", len(rsids), len(genes))

    for rsid in sorted(rsids):
        endpoint = f"variants/{rsid}"
        status, body = fetch(endpoint)
        store_raw(con, run_ts, endpoint, status, body)
        time.sleep(DELAY_SECONDS)
    con.commit()
    return len(rsids)


# --------------------------------------------------------------------------
# Transform: flatten the raw responses into clean tables
# --------------------------------------------------------------------------
def transform(con, run_ts, genes_listed=None):
    """Build drug_gene_pairs, phenotypes, pair_sources, and run_summary from raw data."""
    rows = con.execute(
        "SELECT endpoint, status, body FROM raw_responses "
        "WHERE run_ts = ? AND endpoint LIKE 'pgx/%'", (run_ts,)
    ).fetchall()
    found = missing = pairs = 0

    for endpoint, status, body in rows:
        gene = unquote(endpoint.split("/", 1)[1])
        if status != 200:
            missing += 1
            continue
        found += 1
        data = json.loads(body) or {}
        items = data.get(FIELD_MAP["drug_list_key"]) or []

        # data-quality check: the API reports how many pairs it returned
        if data.get("count") != len(items):
            log.warning("%s: count=%s but parsed %d pairs", gene, data.get("count"), len(items))

        for item in items:
            if not isinstance(item, dict):
                continue
            drug = item.get(FIELD_MAP["drug_name_key"])
            con.execute(
                "INSERT INTO drug_gene_pairs (run_ts, gene, drug, guideline, level, summary) "
                "VALUES (?,?,?,?,?,?)",
                (run_ts, gene, drug,
                 item.get(FIELD_MAP["guideline_key"]),
                 item.get(FIELD_MAP["level"]),
                 item.get(FIELD_MAP["summary"])),
            )
            for ph in as_list(item.get(FIELD_MAP["phenotypes"])):
                ph_name, ph_code, ph_desc = parse_phenotype(ph)
                con.execute(
                    "INSERT INTO phenotypes "
                    "(run_ts, gene, drug, phenotype, phenotype_name, code, description) "
                    "VALUES (?,?,?,?,?,?,?)",
                    (run_ts, gene, drug, as_text(ph), ph_name, ph_code, ph_desc),
                )
            for src in as_list(item.get(FIELD_MAP["sources"])):
                label, url = parse_source(src)
                con.execute(
                    "INSERT INTO pair_sources (run_ts, gene, drug, source, url) VALUES (?,?,?,?,?)",
                    (run_ts, gene, drug, label, url),
                )
            pairs += 1

    con.execute(
        "INSERT INTO run_summary (run_ts, genes_found, genes_missing, pairs, genes_listed) "
        "VALUES (?,?,?,?,?)",
        (run_ts, found, missing, pairs, genes_listed),
    )
    con.commit()
    log.info("genes listed=%s, with /pgx data=%d, without=%d, drug-gene pairs=%d",
             genes_listed, found, missing, pairs)


def transform_variants(con, run_ts):
    """Build the variants and variant_sources tables from the stored raw responses."""
    rows = con.execute(
        "SELECT endpoint, status, body FROM raw_responses WHERE run_ts = ? "
        "AND (endpoint LIKE 'genes/%' OR endpoint LIKE 'variants/%')", (run_ts,)
    ).fetchall()

    rsid_genes, records = {}, {}   # rsID -> genes it was found under; rsID -> its record
    for endpoint, status, body in rows:
        if status != 200:
            continue
        kind, _, rest = endpoint.partition("/")
        body = json.loads(body)
        if kind == "genes":
            for rsid in find_rsids(body):
                rsid_genes.setdefault(rsid, set()).add(unquote(rest))
        elif isinstance(body, dict):
            # some responses wrap the record in a "variant" key; handle both layouts
            records[rest] = body.get("variant") if isinstance(body.get("variant"), dict) else body

    if not rsid_genes:
        log.info("No raw variant data for this run.")
        return 0

    stored = 0
    for rsid, gene_set in sorted(rsid_genes.items()):
        data = records.get(rsid)
        if data is None:
            log.warning("variants/%s: no data (skipped)", rsid)
            continue
        for gene in sorted(gene_set):
            con.execute(
                "INSERT INTO variants (run_ts, rsid, gene, name, alleles, effect, frequency) "
                "VALUES (?,?,?,?,?,?,?)",
                (run_ts, rsid, gene,
                 as_text(data.get(VARIANT_MAP["name"])),
                 as_text(data.get(VARIANT_MAP["alleles"])),
                 as_text(data.get(VARIANT_MAP["effect"])),
                 as_text(data.get(VARIANT_MAP["frequency"]))),
            )
            for src in as_list(data.get(VARIANT_MAP["sources"])):
                label, url = parse_source(src)
                con.execute(
                    "INSERT INTO variant_sources (run_ts, rsid, gene, source, url) VALUES (?,?,?,?,?)",
                    (run_ts, rsid, gene, label, url),
                )
            stored += 1

    con.execute("UPDATE run_summary SET variants = ? WHERE run_ts = ?", (stored, run_ts))
    con.commit()
    log.info("variants stored=%d", stored)
    return stored


def reprocess(con, run_ts):
    """Rebuild every derived table for one run from its stored raw responses (no API calls).
    Use this after changing FIELD_MAP, VARIANT_MAP, or the table layout.
    Returns False (and changes nothing) if the run's raw responses were pruned."""
    has_raw_pgx = con.execute(
        "SELECT COUNT(*) FROM raw_responses WHERE run_ts = ? AND endpoint LIKE 'pgx/%'",
        (run_ts,)).fetchone()[0] > 0
    if not has_raw_pgx:
        log.warning("%s: no raw responses stored (pruned?), skipping", run_ts)
        return False

    prev = con.execute("SELECT genes_listed FROM run_summary WHERE run_ts = ?", (run_ts,)).fetchone()
    has_raw_variants = con.execute(
        "SELECT COUNT(*) FROM raw_responses WHERE run_ts = ? AND endpoint LIKE 'variants/%'",
        (run_ts,)).fetchone()[0] > 0

    tables = ["drug_gene_pairs", "phenotypes", "pair_sources", "run_summary"]
    if has_raw_variants:
        tables += ["variants", "variant_sources"]
    for table in tables:
        con.execute(f"DELETE FROM {table} WHERE run_ts = ?", (run_ts,))

    transform(con, run_ts, prev[0] if prev else None)
    if has_raw_variants:
        transform_variants(con, run_ts)
    return True


# --------------------------------------------------------------------------
# Quality checks
# --------------------------------------------------------------------------
def quality_checks(con, run_ts, compare=True):
    """Look for empty data and, if compare is True, for drops since the previous run.
    Returns a list of warnings."""
    warnings = []
    cur = con.execute("SELECT genes_listed, pairs, variants FROM run_summary WHERE run_ts = ?",
                      (run_ts,)).fetchone()
    if cur is None:
        return ["the run has no summary row"]
    if not cur[1]:
        warnings.append("no drug-gene pairs were stored")

    prev = con.execute(
        "SELECT genes_listed, pairs, variants FROM run_summary "
        "WHERE run_ts < ? AND pairs > 0 ORDER BY run_ts DESC LIMIT 1", (run_ts,)).fetchone()
    if prev and compare:
        for name, new, old in zip(("genes listed", "drug-gene pairs", "variants"), cur, prev):
            if new is not None and old is not None and new < old:
                warnings.append(f"{name} dropped from {old} to {new} since the previous run")

    total, no_effect = con.execute(
        "SELECT COUNT(*), COALESCE(SUM(effect IS NULL), 0) FROM variants WHERE run_ts = ?",
        (run_ts,)).fetchone()
    if total and no_effect == total:
        warnings.append("every variant is missing its effect; check VARIANT_MAP")
    return warnings


# --------------------------------------------------------------------------
# Commands
# --------------------------------------------------------------------------
def run_pipeline(db_path, genes=None, skip_variants=False):
    """Run extract + transform (+ variants) and return (run_ts, quality warnings)."""
    con = sqlite3.connect(db_path)
    run_ts = now()
    try:
        init_db(con)
        # timestamps have one-second resolution, so never reuse one that already has a summary
        while con.execute("SELECT 1 FROM run_summary WHERE run_ts = ?", (run_ts,)).fetchone():
            time.sleep(1)
            run_ts = now()
        gene_list = extract(con, run_ts, genes)
        transform(con, run_ts, len(gene_list))
        if not skip_variants:
            extract_variants(con, run_ts, gene_list)
            transform_variants(con, run_ts)
        # a hand-picked --genes list is a partial run, so don't compare it with full runs
        warnings = quality_checks(con, run_ts, compare=genes is None)
    except BaseException:
        con.rollback()
        # drop raw rows from a run that never finished, so they don't pile up
        con.execute("DELETE FROM raw_responses WHERE run_ts = ? "
                    "AND run_ts NOT IN (SELECT run_ts FROM run_summary)", (run_ts,))
        con.commit()
        raise
    finally:
        con.close()
    return run_ts, warnings


def show_summary(db_path):
    con = sqlite3.connect(db_path)
    cur = con.execute("SELECT * FROM run_summary ORDER BY run_ts DESC")
    header = [c[0] for c in cur.description]
    rows = [["" if v is None else str(v) for v in r] for r in cur.fetchall()]
    con.close()
    if not rows:
        print("No runs yet. Run: python pipeline.py run")
        return
    widths = [max(len(h), *(len(r[i]) for r in rows)) for i, h in enumerate(header)]
    print("  ".join(h.ljust(w) for h, w in zip(header, widths)))
    for r in rows:
        print("  ".join(v.ljust(w) for v, w in zip(r, widths)))


def inspect(target, name=None, max_chars=3000):
    """Print a raw API response so you can check the field names."""
    endpoint = {
        "genes": "genes",
        "gene": f"genes/{quote(name or 'CYP2C19')}",
        "variant": f"variants/{name or 'rs4244285'}",
        "pgx": f"pgx/{quote(name or 'CYP2C19')}",
    }[target]
    status, body = fetch(endpoint)
    print(f"GET /{endpoint} -> {status}")
    print(json.dumps(body, indent=2)[:max_chars])
    if target == "genes":
        print("\nParsed symbols:", parse_gene_list(body))
    elif target == "gene":
        print("\nrsIDs found:", find_rsids(body))


def reprocess_runs(db_path, run_ts=None, every=False):
    con = sqlite3.connect(db_path)
    init_db(con)
    if every:
        runs = [r[0] for r in con.execute("SELECT run_ts FROM run_summary ORDER BY run_ts")]
    else:
        runs = [run_ts or con.execute("SELECT MAX(run_ts) FROM run_summary").fetchone()[0]]
    for ts in filter(None, runs):
        log.info("Reprocessing %s", ts)
        reprocess(con, ts)
    con.close()


def prune(db_path, keep, keep_raw=None):
    """Keep only the newest `keep` runs. With keep_raw, also drop the bulky raw API responses
    for all but the newest `keep_raw` runs while keeping those runs' clean tables."""
    con = sqlite3.connect(db_path)
    kept = [r[0] for r in con.execute(
        "SELECT run_ts FROM run_summary ORDER BY run_ts DESC LIMIT ?", (keep,))]
    if not kept:
        log.info("Nothing to prune.")
        con.close()
        return
    marks = ",".join("?" * len(kept))
    for table in TABLES:
        cur = con.execute(f"DELETE FROM {table} WHERE run_ts NOT IN ({marks})", kept)
        log.info("%s: removed %d rows", table, cur.rowcount)

    if keep_raw is not None:
        raw_kept = kept[:max(keep_raw, 0)]  # kept is newest-first
        if raw_kept:
            cur = con.execute("DELETE FROM raw_responses WHERE run_ts NOT IN (%s)"
                              % ",".join("?" * len(raw_kept)), raw_kept)
        else:
            cur = con.execute("DELETE FROM raw_responses")
        log.info("raw_responses: removed %d more rows (keeping raw data for the newest %d runs)",
                 cur.rowcount, len(raw_kept))
    con.commit()
    con.execute("VACUUM")
    con.close()


# --------------------------------------------------------------------------
# Command line
# --------------------------------------------------------------------------
def build_parser():
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--db", default=DEFAULT_DB, help="path to the SQLite database (default: %(default)s)")
    common.add_argument("--delay", type=float, default=DELAY_SECONDS,
                        help="seconds to wait between API calls (default: %(default)s)")

    parser = argparse.ArgumentParser(description="DeepDNA drug-gene pipeline.")
    sub = parser.add_subparsers(dest="command")

    run = sub.add_parser("run", parents=[common], help="pull data and build the tables")
    run.add_argument("--genes", help="comma-separated gene symbols (default: every gene the API lists)")
    run.add_argument("--skip-variants", action="store_true", help="skip the variants step")
    run.add_argument("--strict", action="store_true", help="exit with code 2 if a quality check warns")

    sub.add_parser("summary", parents=[common], help="show the run history")

    ins = sub.add_parser("inspect", parents=[common], help="print a raw API response")
    ins.add_argument("target", choices=["genes", "gene", "variant", "pgx"])
    ins.add_argument("name", nargs="?", help="gene symbol or rsID (not needed for 'genes')")

    rep = sub.add_parser("reprocess", parents=[common], help="rebuild tables from stored raw data")
    rep.add_argument("--run", help="run timestamp to rebuild (default: the newest)")
    rep.add_argument("--all", action="store_true", help="rebuild every run")

    pr = sub.add_parser("prune", parents=[common], help="delete all but the newest runs")
    pr.add_argument("--keep", type=int, default=10, help="how many runs to keep (default: %(default)s)")
    pr.add_argument("--keep-raw", type=int, default=None,
                    help="also delete raw API responses for all but the newest N runs "
                         "(keeps the clean tables, so history stays small)")
    return parser


def main(argv=None):
    global DELAY_SECONDS
    argv = sys.argv[1:] if argv is None else argv
    args = build_parser().parse_args(argv or ["run"])
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    DELAY_SECONDS = args.delay

    try:
        if args.command in (None, "run"):
            genes = [g.strip() for g in args.genes.split(",") if g.strip()] if getattr(args, "genes", None) else None
            run_ts, warnings = run_pipeline(args.db, genes, getattr(args, "skip_variants", False))
            for w in warnings:
                log.warning("Quality check: %s", w)
            log.info("Run %s finished with %d warning(s)", run_ts, len(warnings))
            return 2 if warnings and getattr(args, "strict", False) else 0
        if args.command == "summary":
            show_summary(args.db)
        elif args.command == "inspect":
            inspect(args.target, args.name)
        elif args.command == "reprocess":
            reprocess_runs(args.db, args.run, args.all)
        elif args.command == "prune":
            prune(args.db, args.keep, args.keep_raw)
        return 0
    except Exception as exc:  # noqa: BLE001 - report any failure with a non-zero exit code
        log.error("%s: %s", type(exc).__name__, exc)
        return 1


if __name__ == "__main__":
    sys.exit(main())
