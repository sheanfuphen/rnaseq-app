"""
src/heatmap.py
──────────────
Interactive heatmap of log₂(TPM+1) Z-scores or raw TPM.

Layout:
  Y axis  — genes (labels draggable via Plotly pan mode)
  X axis  — individual samples (one column per TSV)
  Below X — coloured group annotation bar matching PCA colours
"""

from __future__ import annotations

from typing import Dict, List, Literal, Optional, Tuple

ValueMode = Literal["zscore", "tpm"]

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.colors import sample_colorscale
from plotly.subplots import make_subplots

_FONT = dict(family="Arial, sans-serif", size=12, color="#111111")


def _value_to_fillcolor(
    value: float,
    vmin: float,
    vmax: float,
    colorscale: str = "RdBu_r",
) -> str:
    """Map a numeric value to a fill colour along a colorscale."""
    if vmax <= vmin:
        t = 0.5
    else:
        t = (float(value) - vmin) / (vmax - vmin)
    t = max(0.0, min(1.0, t))
    return sample_colorscale(colorscale, [t])[0]


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
    value_mode: ValueMode = "zscore",
) -> go.Figure:
    """
    Parameters
    ----------
    tpm_df      : genes × samples TPM matrix
    sample_meta : indexed by sample_key, column 'group'
    gene_list   : genes to display (case-insensitive matching against tpm_df.index)
    group_colors: {group_name: hex_color}  — must match PCA colours
    groups      : {group_name: [sample_key, ...]}  — used for column ordering
    value_mode  : "zscore" for per-gene Z-scored log₂(TPM+1), or "tpm" for raw TPM

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

    if value_mode == "tpm":
        vmin = float(sub.values.min()) if sub.size else 0.0
        vmax = float(sub.values.max()) if sub.size else 1.0
        if vmax <= vmin:
            vmax = vmin + 1.0
        colorscale = "Viridis"
        colorbar_title = "TPM"
    else:
        z_abs = float(np.nanmax(np.abs(z.values))) if z.size else 1.0
        z_abs = max(z_abs, 0.01)
        vmin, vmax = -z_abs, z_abs
        colorscale = "RdBu_r"
        colorbar_title = "Z-score<br>(log₂ TPM)"

    # ── One SVG rectangle per cell (editable in Illustrator) ─────────────────
    # Plotly Heatmap exports as a single embedded image in SVG; shapes export
    # as separate <path> elements under individual shape-group nodes.
    hover_x, hover_y, hover_text = [], [], []
    for gi, gene in enumerate(matched):
        for si, samp in enumerate(ordered_samples):
            raw_tpm = float(sub.loc[gene, samp])
            log_val = float(log_tpm.loc[gene, samp])
            zval    = float(z.loc[gene, samp])
            grp     = sample_meta.loc[samp, "group"] if samp in sample_meta.index else "?"
            cell_val = raw_tpm if value_mode == "tpm" else zval
            fig.add_shape(
                type="rect",
                x0=si - 0.5,
                x1=si + 0.5,
                y0=gi - 0.5,
                y1=gi + 0.5,
                fillcolor=_value_to_fillcolor(cell_val, vmin, vmax, colorscale),
                line=dict(width=0),
                layer="below",
                row=1,
                col=1,
            )
            hover_x.append(si)
            hover_y.append(gi)
            hover_text.append(
                f"<b>{gene}</b><br>"
                f"Sample: {samp}<br>"
                f"Group: {grp}<br>"
                f"TPM: {raw_tpm:.2f}<br>"
                f"log₂(TPM+1): {log_val:.3f}<br>"
                f"Z-score: {zval:.3f}"
            )

    # Invisible markers for tooltips only (no line segments in SVG export)
    fig.add_trace(
        go.Scatter(
            x=hover_x,
            y=hover_y,
            mode="markers",
            marker=dict(size=14, opacity=0, symbol="square"),
            line=dict(width=0, color="rgba(0,0,0,0)"),
            text=hover_text,
            hoverinfo="text",
            showlegend=False,
        ),
        row=1,
        col=1,
    )

    # Invisible scale trace for the colour bar only
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            marker=dict(
                colorscale=colorscale,
                cmin=vmin,
                cmax=vmax,
                color=[0],
                showscale=True,
                colorbar=dict(
                    title=dict(
                        text=colorbar_title,
                        side="right",
                        font=_FONT,
                    ),
                    tickfont=_FONT,
                    thickness=12,
                    len=0.9,
                    y=0.5,
                    yanchor="middle",
                    outlinewidth=0,
                ),
            ),
            hoverinfo="skip",
            showlegend=False,
        ),
        row=1,
        col=1,
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
                marker=dict(color=color, line=dict(width=0)),
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
        font=_FONT,
        legend=dict(
            title=dict(text="Group", font=_FONT, side="top"),
            orientation="v",
            x=1.08,
            y=0.97,
            font=_FONT,
        ),
        margin=dict(l=130, r=120, t=40, b=20),
        dragmode="pan",
    )

    # Heatmap axes (numeric coordinates so each cell is a unit square)
    fig.update_xaxes(
        range=[-0.5, n_samples - 0.5],
        tickvals=list(range(n_samples)),
        ticktext=ordered_samples,
        tickangle=-40,
        tickfont=_FONT,
        showgrid=False,
        zeroline=False,
        showline=False,
        mirror=False,
        row=1,
        col=1,
    )
    fig.update_yaxes(
        range=[-0.5, n_genes - 0.5],
        tickvals=list(range(n_genes)),
        ticktext=matched,
        title_text="Gene",
        tickfont=_FONT,
        title_font=_FONT,
        autorange="reversed",
        fixedrange=False,
        showgrid=False,
        zeroline=False,
        showline=False,
        mirror=False,
        row=1,
        col=1,
    )

    # Annotation bar axes
    fig.update_xaxes(
        showticklabels=False,
        showgrid=False,
        zeroline=False,
        showline=False,
        mirror=False,
        row=2,
        col=1,
    )
    fig.update_yaxes(
        showticklabels=False,
        showgrid=False,
        zeroline=False,
        showline=False,
        mirror=False,
        range=[0, 1],
        fixedrange=True,
        row=2,
        col=1,
    )

    # Disable axis baselines on every subplot (avoids lines at x=0 / y=0 in SVG)
    fig.update_xaxes(zeroline=False, showline=False, mirror=False, showgrid=False)
    fig.update_yaxes(zeroline=False, showline=False, mirror=False, showgrid=False)

    # ── Warning for unmatched genes ───────────────────────────────────────────
    if missing:
        shown = ", ".join(missing[:8]) + ("…" if len(missing) > 8 else "")
        fig.add_annotation(
            text=f"⚠️ {len(missing)} gene(s) not found in TPM matrix: {shown}",
            xref="paper", yref="paper",
            x=0, y=-0.03,
            showarrow=False,
            font=dict(family="Arial, sans-serif", size=12, color="#cc0033"),
            xanchor="left",
        )

    return fig
