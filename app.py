"""DeepDNA drug-gene explorer (Streamlit).

Run:   streamlit run app.py
Reads deepdna.db, the SQLite file built by the pipeline notebook. Set the
DEEPDNA_DB environment variable to point at a database somewhere else.
Educational information only, not medical advice.
"""
import json
import os
import sqlite3
from pathlib import Path
from urllib.parse import urlparse

import altair as alt
import pandas as pd
import streamlit as st

DB_PATH = Path(os.environ.get("DEEPDNA_DB", Path(__file__).parent / "deepdna.db"))

# Cool, low-key palette: teal first, amber for contrast, then muted tones.
PALETTE = ["#1F6F8B", "#C9822B", "#7A5C99", "#3E8E6B", "#8A97A3"]
DISCLAIMER = ("Educational information only, not medical advice or diagnosis. "
              "Data comes from the DeepDNA API, curated from sources such as CPIC and DPWG.")

st.set_page_config(page_title="DeepDNA drug-gene explorer", layout="wide")


# --------------------------------------------------------------------------
# Data access
# --------------------------------------------------------------------------
if not DB_PATH.exists():
    st.error(f"Database not found at `{DB_PATH}`.")
    st.write("Run the pipeline notebook first so it creates `deepdna.db`, then put the file "
             "next to `app.py`, or set the `DEEPDNA_DB` environment variable to its path.")
    st.stop()

DB_MTIME = DB_PATH.stat().st_mtime  # part of the cache key, so new runs show up on refresh


@st.cache_data(show_spinner=False)
def _query(sql, params, mtime):
    """Run a read-only query. A missing table or column gives an empty table."""
    try:
        con = sqlite3.connect(f"file:{DB_PATH}?mode=ro", uri=True)
        try:
            return pd.read_sql(sql, con, params=params)
        finally:
            con.close()
    except (sqlite3.OperationalError, pd.errors.DatabaseError):
        return pd.DataFrame()


def q(sql, params=()):
    return _query(sql, tuple(params), DB_MTIME)


def ensure_columns(df, columns):
    """Make sure an empty or older table still has the columns the app reads."""
    for col in columns:
        if col not in df.columns:
            df[col] = pd.Series(dtype="object")
    return df


@st.cache_data(show_spinner=False)
def load_catalog(run_ts, mtime):
    """Gene names, categories, and summaries from the raw /genes response."""
    raw = _query("SELECT body FROM raw_responses WHERE run_ts = ? AND endpoint = 'genes' LIMIT 1",
                 (run_ts,), mtime)
    empty = pd.DataFrame(columns=["symbol", "name", "category", "summary"])
    if raw.empty:
        return empty
    body = json.loads(raw.body.iloc[0]) or {}
    items = body.get("genes") if isinstance(body, dict) else body
    rows = [{"symbol": g.get("symbol"),
             "name": g.get("name"),
             "category": str(g.get("category") or "uncategorized").replace("_", " ").capitalize(),
             "summary": g.get("summary")}
            for g in (items or []) if isinstance(g, dict) and g.get("symbol")]
    return pd.DataFrame(rows) if rows else empty


# --------------------------------------------------------------------------
# Charts
# --------------------------------------------------------------------------
def hbar(df, value, label, x_title, color=None):
    """Horizontal bar chart with whole-number ticks, largest bar on top."""
    encode = dict(
        x=alt.X(f"{value}:Q", title=x_title, axis=alt.Axis(tickMinStep=1, format="d")),
        y=alt.Y(f"{label}:N", sort="-x", title=None,
                axis=alt.Axis(labelLimit=240, labelOverlap=False)),
        tooltip=list(df.columns),
    )
    if color:
        cats = sorted(df[color].dropna().unique())
        encode["color"] = alt.Color(f"{color}:N", title=None,
                                    scale=alt.Scale(domain=cats, range=PALETTE[:len(cats)] or PALETTE),
                                    legend=alt.Legend(orient="top"))
    else:
        encode["color"] = alt.value(PALETTE[0])
    return (alt.Chart(df).mark_bar().encode(**encode)
            .properties(width="container", height=alt.Step(28)))


def domain_counts(df):
    """Count citations per website."""
    if df.empty:
        return pd.DataFrame(columns=["site", "citations"])
    site = df["url"].map(lambda u: urlparse(u).netloc if isinstance(u, str) and u else None)
    counts = site.fillna("no link").value_counts().head(10)
    return counts.rename_axis("site").reset_index(name="citations")


def source_links(df):
    """Render sources as a bulleted list of links."""
    if df.empty:
        st.caption("No sources listed.")
        return
    for s in df.drop_duplicates(["source", "url"]).itertuples():
        st.markdown(f"- [{s.source}]({s.url})" if isinstance(s.url, str) and s.url else f"- {s.source}")


# --------------------------------------------------------------------------
# Sidebar: pick a run and filter genes
# --------------------------------------------------------------------------
runs = q("SELECT * FROM run_summary WHERE pairs > 0 ORDER BY run_ts DESC")
if runs.empty:
    st.warning("No completed runs in the database yet. Run the pipeline notebook, then refresh.")
    st.stop()
runs = ensure_columns(runs, ["genes_listed", "variants"])

with st.sidebar:
    st.header("Filters")
    run_ts = st.selectbox(
        "Snapshot", list(runs.run_ts),
        format_func=lambda t: pd.to_datetime(t).strftime("%d %b %Y, %H:%M UTC"),
        help="Each pipeline run saves a snapshot of the API. The newest is selected.")
    st.caption(DISCLAIMER)

pairs_all = ensure_columns(q("SELECT * FROM drug_gene_pairs WHERE run_ts = ?", (run_ts,)),
                           ["gene", "drug", "guideline", "level", "summary"])
phen = ensure_columns(q("SELECT * FROM phenotypes WHERE run_ts = ?", (run_ts,)),
                      ["gene", "drug", "phenotype_name", "code", "description"])
pair_sources = ensure_columns(q("SELECT * FROM pair_sources WHERE run_ts = ?", (run_ts,)),
                              ["gene", "drug", "source", "url"])
variants_all = ensure_columns(q("SELECT * FROM variants WHERE run_ts = ?", (run_ts,)),
                              ["rsid", "gene", "name", "alleles", "effect", "frequency"])
variant_sources = ensure_columns(q("SELECT * FROM variant_sources WHERE run_ts = ?", (run_ts,)),
                                 ["rsid", "gene", "source", "url"])
catalog = load_catalog(run_ts, DB_MTIME)

all_genes = sorted(set(pairs_all.gene) | set(variants_all.gene) | set(catalog.symbol))
with st.sidebar:
    chosen = st.multiselect("Genes", all_genes, placeholder="All genes")
genes_in_view = chosen or all_genes

pairs = pairs_all[pairs_all.gene.isin(genes_in_view)].copy()
variants = variants_all[variants_all.gene.isin(genes_in_view)].copy()
pair_sources_v = pair_sources[pair_sources.gene.isin(genes_in_view)]
variant_sources_v = variant_sources[variant_sources.gene.isin(genes_in_view)]

category = catalog.set_index("symbol").category.to_dict()
pairs["category"] = pairs.gene.map(category).fillna("All genes")

# --------------------------------------------------------------------------
# Header and key numbers
# --------------------------------------------------------------------------
st.title("Drug-gene explorer")
st.caption(f"Snapshot from {pd.to_datetime(run_ts).strftime('%d %b %Y, %H:%M UTC')}. {DISCLAIMER}")

links = pd.concat([pair_sources_v.url, variant_sources_v.url]).dropna()
k1, k2, k3, k4, k5 = st.columns(5)
k1.metric("Genes", len(genes_in_view))
k2.metric("Genes with drug pairs", pairs.gene.nunique())
k3.metric("Drug-gene pairs", len(pairs))
k4.metric("Key variants", len(variants))
k5.metric("Unique source links", links.nunique())

tab_overview, tab_pairs, tab_variants, tab_evidence, tab_history = st.tabs(
    ["Overview", "Drug-gene pairs", "Variants", "Evidence", "Run history"])

# --------------------------------------------------------------------------
# Overview
# --------------------------------------------------------------------------
with tab_overview:
    left, right = st.columns(2)
    with left:
        st.subheader("Drug pairs per gene")
        per_gene = pairs.groupby(["gene", "category"]).size().reset_index(name="pairs")
        if per_gene.empty:
            st.info("No drug pairs for the selected genes.")
        else:
            st.altair_chart(hbar(per_gene, "pairs", "gene", "Drug pairs", color="category"))
    with right:
        st.subheader("Evidence level")
        levels = pairs.level.fillna("Not provided").value_counts().rename_axis("level").reset_index(name="pairs")
        if levels.empty:
            st.info("No evidence levels to show.")
        else:
            st.altair_chart(hbar(levels, "pairs", "level", "Pairs"))

    st.subheader("Where the guidance comes from")
    n = len(pairs)
    if n:
        flags = pd.DataFrame({
            "source": ["CPIC", "DPWG", "FDA label"],
            "pairs": [pairs.guideline.str.contains(key, na=False).sum() for key in ("CPIC", "DPWG", "FDA")],
        })
        flags["share of pairs"] = (flags.pairs / n).map("{:.0%}".format)
        c1, c2 = st.columns([2, 1])
        c1.altair_chart(hbar(flags, "pairs", "source", "Pairs"))
        c2.dataframe(flags, hide_index=True)
        missing = int(pairs.guideline.isna().sum())
        st.caption(f"Shares use all {n} pairs and can add up to more than 100%, because a pair can cite "
                   f"more than one source. {missing} {'pair lists' if missing == 1 else 'pairs list'} no guideline.")

# --------------------------------------------------------------------------
# Drug-gene pairs
# --------------------------------------------------------------------------
with tab_pairs:
    search = st.text_input("Search", placeholder="A drug or gene, such as clopidogrel or CYP2D6")
    view = pairs
    if search:
        view = pairs[pairs.drug.str.contains(search, case=False, na=False)
                     | pairs.gene.str.contains(search, case=False, na=False)]

    table = view[["gene", "drug", "guideline", "level"]].sort_values(["gene", "drug"])
    st.dataframe(table, hide_index=True)
    st.download_button("Download table as CSV", table.to_csv(index=False),
                       file_name="drug_gene_pairs.csv", mime="text/csv")

    if view.empty:
        st.info("No pairs match that search.")
    else:
        st.subheader("Pair details")
        options = sorted((view.gene + " / " + view.drug).tolist())
        choice = st.selectbox("Choose a pair", options)
        gene, drug = choice.split(" / ", 1)
        row = view[(view.gene == gene) & (view.drug == drug)].iloc[0]

        st.markdown(f"**{gene} and {drug}**")
        st.caption(f"Guideline: {row.guideline if pd.notna(row.guideline) else 'not listed'}. "
                   f"Evidence level: {row.level if pd.notna(row.level) else 'not listed'}.")
        if pd.notna(row.summary):
            st.write(row.summary)

        ph = phen[(phen.gene == gene) & (phen.drug == drug)]
        if not ph.empty:
            st.markdown("**Phenotype groups**")
            st.dataframe(ph[["phenotype_name", "code", "description"]]
                         .rename(columns={"phenotype_name": "phenotype"}), hide_index=True)
        st.markdown("**Sources**")
        source_links(pair_sources[(pair_sources.gene == gene) & (pair_sources.drug == drug)])

    with st.expander("Drugs that depend on more than one gene"):
        multi = (pairs_all.groupby("drug")
                 .agg(genes=("gene", lambda s: ", ".join(sorted(set(s)))), gene_count=("gene", "nunique"))
                 .query("gene_count > 1")
                 .sort_values(["gene_count", "drug"], ascending=[False, True])
                 .reset_index())
        if multi.empty:
            st.write("No drug in this snapshot is linked to more than one gene.")
        else:
            st.dataframe(multi, hide_index=True)

# --------------------------------------------------------------------------
# Variants
# --------------------------------------------------------------------------
with tab_variants:
    if variants.empty:
        st.info("No variants stored for this snapshot. Run the variant steps in the pipeline notebook.")
    else:
        st.subheader("Key variants per gene")
        per_gene_v = variants.groupby("gene").size().reset_index(name="variants")
        st.altair_chart(hbar(per_gene_v, "variants", "gene", "Variants"))

        st.subheader("Look up a gene")
        gene = st.selectbox("Gene", sorted(variants.gene.unique()))
        info = catalog[catalog.symbol == gene]
        if not info.empty:
            st.caption(f"{info.name.iloc[0]}. Category: {info.category.iloc[0].lower()}.")
            if pd.notna(info.summary.iloc[0]):
                st.write(info.summary.iloc[0])

        gv = variants[variants.gene == gene].sort_values(["name", "rsid"])
        st.dataframe(gv[["rsid", "name", "effect", "frequency"]]
                     .rename(columns={"name": "variant"}), hide_index=True)
        if gv.effect.isna().all():
            st.caption("No functional effects were stored for this gene. If that looks wrong, check "
                       "`VARIANT_MAP` in the notebook and rebuild with `reprocess`.")

        for v in gv.itertuples():
            with st.expander(f"{v.rsid}  {v.name if pd.notna(v.name) else ''}"):
                st.markdown("**Alleles**")
                st.write(v.alleles if pd.notna(v.alleles) else "Not listed.")
                st.markdown("**Functional effect**")
                st.write(v.effect if pd.notna(v.effect) else "Not listed.")
                st.markdown("**Frequency**")
                st.write(v.frequency if pd.notna(v.frequency) else "Not listed.")
                st.markdown("**Sources**")
                source_links(variant_sources[(variant_sources.rsid == v.rsid) & (variant_sources.gene == v.gene)])

# --------------------------------------------------------------------------
# Evidence
# --------------------------------------------------------------------------
with tab_evidence:
    st.write("Which websites the citations point to. Each citation counts once per pair or variant it supports.")
    e1, e2 = st.columns(2)
    with e1:
        st.subheader("Drug-gene pairs")
        pc = domain_counts(pair_sources_v)
        if pc.empty:
            st.info("No citations stored for pairs.")
        else:
            st.altair_chart(hbar(pc, "citations", "site", "Citations"))
    with e2:
        st.subheader("Variants")
        vc = domain_counts(variant_sources_v)
        if vc.empty:
            st.info("No citations stored for variants.")
        else:
            st.altair_chart(hbar(vc, "citations", "site", "Citations"))

    st.subheader("Metabolizer codes")
    codes = phen[phen.gene.isin(genes_in_view)].code.dropna()
    if codes.empty:
        st.info("No phenotype codes stored for the selected genes.")
    else:
        code_counts = codes.value_counts().head(12).rename_axis("code").reset_index(name="rows")
        st.altair_chart(hbar(code_counts, "rows", "code", "Phenotype rows"))
        st.caption("Codes come from text like \"Poor metabolizer (PM): ...\". Genes that describe "
                   "genotypes or activity scores show those labels instead.")

# --------------------------------------------------------------------------
# Run history
# --------------------------------------------------------------------------
with tab_history:
    history = q("SELECT * FROM run_summary ORDER BY run_ts")
    history = ensure_columns(history, ["genes_listed", "variants"])
    st.write("One row per pipeline run. The data is curated and changes slowly, so expect flat lines "
             "until the API adds genes or pairs.")
    st.dataframe(history.sort_values("run_ts", ascending=False), hide_index=True)

    long = (history.assign(run=pd.to_datetime(history.run_ts))
            .melt(id_vars="run", value_vars=["genes_listed", "pairs", "variants"],
                  var_name="measure", value_name="count").dropna())
    if long.run.nunique() > 1:
        chart = (alt.Chart(long).mark_line(point=True)
                 .encode(x=alt.X("run:T", title=None, axis=alt.Axis(format="%d %b, %H:%M")),
                         y=alt.Y("count:Q", title="Count", axis=alt.Axis(tickMinStep=1, format="d")),
                         color=alt.Color("measure:N", title=None,
                                         scale=alt.Scale(range=PALETTE[:3])),
                         tooltip=["run:T", "measure", "count"])
                 .properties(width="container", height=280))
        st.altair_chart(chart)
    else:
        st.caption("Run the pipeline again to see how coverage changes over time.")
