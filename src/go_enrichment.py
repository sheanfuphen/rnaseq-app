"""
src/go_enrichment.py
────────────────────
Gene Ontology enrichment analysis using the gseapy library (Enrichr API).

For a given gene list (e.g. significant DEGs from one comparison),
queries Enrichr for enriched GO terms and returns a bar/dot plot.

Requires: gseapy >= 1.0
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st


# ── Supported gene set libraries ─────────────────────────────────────────────
GO_LIBRARIES = {
    "GO Biological Process": "GO_Biological_Process_2023",
    "GO Molecular Function": "GO_Molecular_Function_2023",
    "GO Cellular Component": "GO_Cellular_Component_2023",
    "KEGG Pathways":         "KEGG_2021_Human",
    "Reactome":              "Reactome_2022",
    "MSigDB Hallmarks":      "MSigDB_Hallmark_2020",
}

GO_BP_LIBRARY = GO_LIBRARIES["GO Biological Process"]


def load_go_biological_process_genesets() -> Dict[str, List[str]]:
    """
    Download (or reuse the local Enrichr cache of) GO Biological Process
    gene sets. Returns {term_name: [gene, ...]}.
    """
    import gseapy as gp

    return gp.get_library(
        name=GO_BP_LIBRARY,
        organism="Human",
        min_size=0,
        max_size=100_000,
    )


def run_enrichment(
    gene_list: List[str],
    library: str = "GO_Biological_Process_2023",
    organism: str = "human",
    cutoff: float = 0.05,
    top_n: int = 20,
) -> Optional[pd.DataFrame]:
    """
    Run Enrichr-based GO enrichment via gseapy.

    Parameters
    ----------
    gene_list : list of gene symbols (HGNC)
    library   : Enrichr gene set library name
    organism  : 'human' or 'mouse'
    cutoff    : adjusted p-value cutoff
    top_n     : maximum number of terms to return

    Returns
    -------
    DataFrame with columns:
        Term, Overlap, P-value, Adjusted P-value, Odds Ratio,
        Combined Score, Genes, -log10(padj), Gene Count
    or None if gseapy is not installed / no results
    """
    try:
        import gseapy as gp
    except ImportError:
        st.error(
            "gseapy is not installed. Run `pip install gseapy` and restart the app."
        )
        return None

    if not gene_list:
        st.warning("No genes provided for enrichment analysis.")
        return None

    try:
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=library,
            organism=organism,
            cutoff=cutoff,
            outdir=None,
            verbose=False,
        )
    except Exception as e:
        st.error(f"Enrichr query failed: {e}")
        return None

    if enr is None or enr.results is None or enr.results.empty:
        st.warning("No enriched terms found at the selected cutoff.")
        return None

    res = enr.results.copy()

    # Standardise column names across gseapy versions
    res.columns = [c.strip() for c in res.columns]
    rename = {
        "Term":               "Term",
        "Overlap":            "Overlap",
        "P-value":            "P-value",
        "Adjusted P-value":   "Adjusted P-value",
        "Odds Ratio":         "Odds Ratio",
        "Combined Score":     "Combined Score",
        "Genes":              "Genes",
    }
    res = res.rename(columns={k: v for k, v in rename.items() if k in res.columns})

    res = res[res["Adjusted P-value"] <= cutoff].copy()
    if res.empty:
        st.warning("No enriched terms found at the selected cutoff.")
        return None

    res["-log10(padj)"] = -np.log10(res["Adjusted P-value"].clip(lower=1e-300))

    # Parse gene count from Overlap string e.g. "5/200"
    if "Overlap" in res.columns:
        res["Gene Count"] = res["Overlap"].str.split("/").str[0].astype(int)
    else:
        res["Gene Count"] = 1

    res = res.sort_values("Adjusted P-value").head(top_n)
    return res


def plot_go_bars(
    res: pd.DataFrame,
    title: str = "GO Enrichment",
    color: str = "#4361ee",
    width: int = 800,
    height: int = 600,
) -> go.Figure:
    """
    Horizontal bar chart of top enriched GO terms.
    Bars are coloured by -log10(adjusted p-value).
    """
    df = res.sort_values("-log10(padj)", ascending=True).copy()

    # Truncate long term names
    df["Term_short"] = df["Term"].str.replace(r"\(GO:\d+\)", "", regex=True).str.strip()
    df["Term_short"] = df["Term_short"].apply(lambda x: x[:55] + "…" if len(x) > 55 else x)

    hover = (
        df["Term"] + "<br>"
        + "Adjusted p-value: " + df["Adjusted P-value"].map("{:.2e}".format) + "<br>"
        + "-log10(padj): "     + df["-log10(padj)"].round(2).astype(str) + "<br>"
        + "Gene Count: "       + df["Gene Count"].astype(str) + "<br>"
        + "Genes: "            + df["Genes"].astype(str)
    )

    fig = go.Figure(go.Bar(
        x=df["-log10(padj)"],
        y=df["Term_short"],
        orientation="h",
        marker=dict(
            color=df["-log10(padj)"],
            colorscale=[[0, _lighten(color, 0.6)], [1, color]],
            showscale=True,
            colorbar=dict(
                title=dict(text="-log₁₀(padj)", side="right"),
                thickness=12,
            ),
        ),
        hovertext=hover,
        hoverinfo="text",
        text=df["Gene Count"].astype(str) + " genes",
        textposition="outside",
        textfont=dict(size=11, color="#111111"),
    ))

    fig.add_vline(x=-np.log10(0.05), line_dash="dot",
                  line_color="#888", line_width=1,
                  annotation_text="padj=0.05",
                  annotation_font=dict(size=11, color="#555"))

    fig.update_layout(
        title=dict(text=title, font=dict(size=16, color="#111111")),
        xaxis_title="-log₁₀(adjusted p-value)",
        xaxis=dict(tickfont=dict(size=12, color="#111111"),
                   title_font=dict(size=13, color="#111111")),
        yaxis=dict(tickfont=dict(size=11, color="#111111"), automargin=True),
        plot_bgcolor="#f5f5f5",
        paper_bgcolor="#ffffff",
        font=dict(family="DM Sans, sans-serif", size=13, color="#111111"),
        margin=dict(l=300, r=80, t=60, b=60),
        width=width,
        height=max(height, len(df) * 32 + 120),
    )
    return fig


def plot_go_dots(
    res: pd.DataFrame,
    title: str = "GO Enrichment",
    color: str = "#4361ee",
    width: int = 800,
    height: int = 600,
) -> go.Figure:
    """
    Dot plot: X = -log10(padj), dot size = gene count, colour intensity = odds ratio.
    """
    df = res.sort_values("-log10(padj)", ascending=True).copy()
    df["Term_short"] = df["Term"].str.replace(r"\(GO:\d+\)", "", regex=True).str.strip()
    df["Term_short"] = df["Term_short"].apply(lambda x: x[:55] + "…" if len(x) > 55 else x)

    hover = (
        df["Term"] + "<br>"
        + "Adjusted p-value: " + df["Adjusted P-value"].map("{:.2e}".format) + "<br>"
        + "Gene Count: "       + df["Gene Count"].astype(str) + "<br>"
        + "Odds Ratio: "       + df["Odds Ratio"].round(2).astype(str) + "<br>"
        + "Genes: "            + df["Genes"].astype(str)
    )

    # Scale dot sizes between 8 and 28
    gc      = df["Gene Count"]
    gc_norm = (gc - gc.min()) / (gc.max() - gc.min() + 1e-9)
    sizes   = (gc_norm * 20 + 8).tolist()

    fig = go.Figure(go.Scatter(
        x=df["-log10(padj)"],
        y=df["Term_short"],
        mode="markers",
        marker=dict(
            size=sizes,
            color=df["-log10(padj)"],
            colorscale=[[0, _lighten(color, 0.6)], [1, color]],
            showscale=True,
            colorbar=dict(
                title=dict(text="-log₁₀(padj)", side="right"),
                thickness=12,
            ),
            line=dict(width=0.5, color="#ffffff"),
        ),
        hovertext=hover,
        hoverinfo="text",
    ))

    fig.add_vline(x=-np.log10(0.05), line_dash="dot",
                  line_color="#888", line_width=1,
                  annotation_text="padj=0.05",
                  annotation_font=dict(size=11, color="#555"))

    fig.update_layout(
        title=dict(text=title, font=dict(size=16, color="#111111")),
        xaxis_title="-log₁₀(adjusted p-value)",
        xaxis=dict(tickfont=dict(size=12, color="#111111"),
                   title_font=dict(size=13, color="#111111")),
        yaxis=dict(tickfont=dict(size=11, color="#111111"), automargin=True),
        plot_bgcolor="#f5f5f5",
        paper_bgcolor="#ffffff",
        font=dict(family="DM Sans, sans-serif", size=13, color="#111111"),
        margin=dict(l=300, r=80, t=60, b=60),
        width=width,
        height=max(height, len(df) * 32 + 120),
    )
    return fig


def _lighten(hex_color: str, factor: float) -> str:
    """Lighten a hex colour by blending with white."""
    hex_color = hex_color.lstrip("#")
    r, g, b = (int(hex_color[i:i+2], 16) for i in (0, 2, 4))
    r = int(r + (255 - r) * factor)
    g = int(g + (255 - g) * factor)
    b = int(b + (255 - b) * factor)
    return f"#{r:02x}{g:02x}{b:02x}"
