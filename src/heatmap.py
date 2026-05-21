"""
src/heatmap.py
──────────────
Interactive heatmap built on log2(TPM + 1) Z-scores.

Layout:
  Y axis  — genes (labels draggable via Plotly pan mode)
  X axis  — individual samples (one column per TSV)
  Below X — coloured group annotation bar matching PCA colours
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def order_genes_by_group_contrast(
    log_tpm: pd.DataFrame,
    groups: Dict[str, List[str]],
    group1: str,
    group2: str,
    genes: List[str],
) -> List[str]:
    """
    Order genes for heatmap rows: high in group1 (top) → similar (middle) → high in group2 (bottom).

    Uses mean log₂(TPM + 1) per group; sort key is group1_mean − group2_mean.
    """
    samples_g1 = [s for s in groups.get(group1, []) if s in log_tpm.columns]
    samples_g2 = [s for s in groups.get(group2, []) if s in log_tpm.columns]
    if not samples_g1 or not samples_g2:
        return genes

    present = [g for g in genes if g in log_tpm.index]
    if not present:
        return genes

    mean_g1 = log_tpm.loc[present, samples_g1].mean(axis=1)
    mean_g2 = log_tpm.loc[present, samples_g2].mean(axis=1)
    score = mean_g1 - mean_g2
    ordered = score.sort_values(ascending=False).index.tolist()

    missing = [g for g in genes if g not in log_tpm.index]
    return ordered + missing


def order_gene_list_for_heatmap(
    tpm_df: pd.DataFrame,
    groups: Dict[str, List[str]],
    gene_list: List[str],
    contrast_groups: Optional[Tuple[str, str]] = None,
) -> List[str]:
    """Reorder a user gene list using group contrast (or return unchanged if no contrast pair)."""
    idx_upper = {g.upper(): g for g in tpm_df.index}
    matched = [idx_upper[g.upper()] for g in gene_list if g.upper() in idx_upper]
    if not matched:
        return gene_list

    log_tpm = np.log2(tpm_df.loc[matched].astype(float) + 1)
    if contrast_groups:
        g1, g2 = contrast_groups
        matched = order_genes_by_group_contrast(log_tpm, groups, g1, g2, matched)

    canon_to_input = {idx_upper[g.upper()]: g for g in gene_list if g.upper() in idx_upper}
    ordered = [canon_to_input[m] for m in matched]
    not_found = [g for g in gene_list if g.upper() not in idx_upper]
    return ordered + not_found


def build_heatmap(
    tpm_df: pd.DataFrame,
    sample_meta: pd.DataFrame,
    gene_list: List[str],
    group_colors: Dict[str, str],
    groups: Dict[str, List[str]],
) -> go.Figure:
    """
    Parameters
    ----------
    tpm_df      : genes × samples TPM matrix
    sample_meta : indexed by sample_key, column 'group'
    gene_list   : genes to display (case-insensitive matching against tpm_df.index)
    group_colors: {group_name: hex_color}  — must match PCA colours
    groups      : {group_name: [sample_key, ...]}  — used for column ordering

    Returns
    -------
    Plotly Figure (heatmap rows + thin colour annotation bar)
    """

    # ── Case-insensitive gene matching ────────────────────────────────────────
    idx_upper = {g.upper(): g for g in tpm_df.index}
    matched   = [idx_upper[g.upper()] for g in gene_list if g.upper() in idx_upper]
    missing   = [g for g in gene_list if g.upper() not in idx_upper]

    if not matched:
        raise ValueError(
            f"None of the requested genes were found in the TPM matrix. "
            f"Sample of missing: {missing[:10]}"
        )

    # ── Order samples group-by-group ──────────────────────────────────────────
    ordered_samples: List[str] = []
    for gname, samples in groups.items():
        for s in samples:
            if s in sample_meta.index:
                ordered_samples.append(s)

    sub = tpm_df.loc[matched, ordered_samples].astype(float)

    # ── log2(TPM + 1) then Z-score per gene ──────────────────────────────────
    log_tpm  = np.log2(sub + 1)
    row_mean = log_tpm.mean(axis=1)
    row_std  = log_tpm.std(axis=1).replace(0, 1)
    z        = log_tpm.subtract(row_mean, axis=0).divide(row_std, axis=0)

    n_genes   = len(matched)
    n_samples = len(ordered_samples)

    # ── Subplot: heatmap (row 1) + annotation bar (row 2) ────────────────────
    fig = make_subplots(
        rows=2, cols=1,
        row_heights=[0.93, 0.07],
        vertical_spacing=0.01,
        shared_xaxes=True,
    )

    # ── Hover text ────────────────────────────────────────────────────────────
    hover = []
    for gene in matched:
        row = []
        for samp in ordered_samples:
            grp   = sample_meta.loc[samp, "group"] if samp in sample_meta.index else "?"
            raw   = round(float(sub.loc[gene, samp]), 2)
            zval  = round(float(z.loc[gene, samp]), 3)
            row.append(
                f"<b>{gene}</b><br>"
                f"Sample: {samp}<br>"
                f"Group: {grp}<br>"
                f"TPM: {raw}<br>"
                f"Z-score: {zval}"
            )
        hover.append(row)

    # ── Heatmap trace ─────────────────────────────────────────────────────────
    fig.add_trace(
        go.Heatmap(
            z=z.values,
            x=ordered_samples,
            y=matched,
            text=hover,
            hoverinfo="text",
            colorscale="RdBu_r",
            zmid=0,
            colorbar=dict(
                title=dict(text="Z-score<br>(log₂ TPM)", side="right"),
                thickness=12,
                len=0.9,
                y=0.5,
                yanchor="middle",
            ),
        ),
        row=1, col=1,
    )

    # ── Group annotation bar ──────────────────────────────────────────────────
    # Track which groups have been added to legend already
    legend_groups_added: set = set()

    for samp in ordered_samples:
        grp   = sample_meta.loc[samp, "group"] if samp in sample_meta.index else "?"
        color = group_colors.get(grp, "#aaaaaa")
        show_legend = grp not in legend_groups_added
        legend_groups_added.add(grp)

        fig.add_trace(
            go.Bar(
                x=[samp],
                y=[1],
                marker_color=color,
                name=grp,
                showlegend=show_legend,
                legendgroup=grp,
                hovertext=grp,
                hoverinfo="text",
            ),
            row=2, col=1,
        )

    # ── Layout ────────────────────────────────────────────────────────────────
    gene_px    = max(16, min(34, int(650 / max(n_genes, 1))))
    fig_height = max(400, gene_px * n_genes + 100)

    fig.update_layout(
        height=fig_height,
        barmode="stack",
        plot_bgcolor="white",
        paper_bgcolor="white",
        font=dict(family="DM Sans, sans-serif", size=11, color="#111111"),
        legend=dict(
            title=dict(
                text="Group",
                font=dict(color="#111111", size=12),
                side="top",
            ),
            orientation="v",
            x=1.08, y=0.97,
            font=dict(color="#111111"),
        ),
        margin=dict(l=130, r=120, t=40, b=20),
        dragmode="pan",
    )

    # Heatmap axes
    fig.update_xaxes(
        tickangle=-40, tickfont=dict(size=10, color="#111111"),
        showgrid=False, row=1, col=1,
    )
    fig.update_yaxes(
        title_text="Gene",
        tickfont=dict(size=10, color="#111111"),
        title_font=dict(color="#111111"),
        autorange="reversed",
        fixedrange=False,
        showgrid=False, row=1, col=1,
    )

    # Annotation bar axes
    fig.update_xaxes(showticklabels=False, showgrid=False, row=2, col=1)
    fig.update_yaxes(
        showticklabels=False, showgrid=False,
        range=[0, 1], fixedrange=True, row=2, col=1,
    )

    # ── Warning for unmatched genes ───────────────────────────────────────────
    if missing:
        shown = ", ".join(missing[:8]) + ("…" if len(missing) > 8 else "")
        fig.add_annotation(
            text=f"⚠️ {len(missing)} gene(s) not found in TPM matrix: {shown}",
            xref="paper", yref="paper",
            x=0, y=-0.03,
            showarrow=False,
            font=dict(size=10, color="#cc0033"),
            xanchor="left",
        )

    return fig
