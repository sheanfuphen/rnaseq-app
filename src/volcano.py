"""
src/volcano.py
──────────────
Volcano plot for one pairwise DGE comparison.

X axis: log2FoldChange  (positive → higher in group1 / right side)
Y axis: -log10(padj)

Group name labels appear at the top-left and top-right corners.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go

_FONT = dict(family="Arial, sans-serif", size=12, color="#111111")


def plot_volcano(
    df: pd.DataFrame,
    group1: str,
    group2: str,
    fc_thresh: float = 1.0,
    neg_log10_padj_thresh: float = 1.301,   # ≈ padj < 0.05
    color_up: str = "#4361ee",
    color_down: str = "#f72585",
) -> go.Figure:
    """
    Parameters
    ----------
    df                    : DGE results (index=gene, cols include
                            log2FoldChange, padj, optionally baseMean)
    group1                : name for the right (positive FC) side
    group2                : name for the left (negative FC) side
    fc_thresh             : absolute log2FC cut-off for colouring
    neg_log10_padj_thresh : significance cut-off on Y axis
    color_up              : colour for genes up in group1
    color_down            : colour for genes up in group2
    """
    plot_df = df.copy()
    plot_df = plot_df[
        np.isfinite(plot_df["log2FoldChange"]) & (plot_df["padj"] > 0)
    ].copy()
    plot_df["-log10padj"] = -np.log10(plot_df["padj"])

    def _classify(row):
        sig = row["-log10padj"] >= neg_log10_padj_thresh
        up  = row["log2FoldChange"] >=  fc_thresh
        dn  = row["log2FoldChange"] <= -fc_thresh
        if sig and up:   return "up"
        if sig and dn:   return "down"
        if sig:          return "sig_only"
        return "ns"

    plot_df["class"] = plot_df.apply(_classify, axis=1)

    palette  = {"up": color_up, "down": color_down, "sig_only": "#aaa", "ns": "#ddd"}
    sizes    = {"up": 7,        "down": 7,           "sig_only": 5,      "ns": 4}
    opacs    = {"up": 0.85,     "down": 0.85,        "sig_only": 0.6,    "ns": 0.35}

    fig = go.Figure()

    for cls, label in [
        ("ns",       "Not significant"),
        ("sig_only", "Significant (low |FC|)"),
        ("down",     f"Up in {group2}"),
        ("up",       f"Up in {group1}"),
    ]:
        sub = plot_df[plot_df["class"] == cls]
        if sub.empty:
            continue

        hover = (
            sub.index
            + "<br>log2FC: " + sub["log2FoldChange"].round(3).astype(str)
            + "<br>padj: "   + sub["padj"].map("{:.2e}".format)
        )
        if "baseMean" in sub.columns:
            hover += "<br>baseMean: " + sub["baseMean"].round(1).astype(str)

        # Label genes in the up / down groups (user-coloured points)
        show_labels = cls in ("up", "down")
        mode = "markers+text" if show_labels else "markers"
        text = sub.index.astype(str) if show_labels else None
        textfont = dict(family="Arial, sans-serif", size=12, color=palette[cls]) if show_labels else None

        fig.add_trace(go.Scatter(
            x=sub["log2FoldChange"],
            y=sub["-log10padj"],
            mode=mode,
            name=label,
            text=text,
            textposition="top center",
            textfont=textfont,
            hovertext=hover,
            hoverinfo="text",
            marker=dict(color=palette[cls], size=sizes[cls],
                        opacity=opacs[cls], line=dict(width=0)),
        ))

    x_range = plot_df["log2FoldChange"].abs().max() * 1.1
    y_max   = plot_df["-log10padj"].max() * 1.1

    # Threshold reference lines
    fig.add_hline(y=neg_log10_padj_thresh, line_dash="dot",
                  line_color="#999", line_width=1)
    for xv in [-fc_thresh, fc_thresh]:
        fig.add_vline(x=xv, line_dash="dot", line_color="#999", line_width=1)

    # ── Group labels at top corners ───────────────────────────────────────────
    label_y = y_max * 0.97
    fig.add_annotation(
        x=-x_range * 0.95, y=label_y,
        text=f"← {group2}",
        showarrow=False,
        font=dict(size=12, color=color_down, family="Arial, sans-serif"),
        xanchor="left",
    )
    fig.add_annotation(
        x=x_range * 0.95, y=label_y,
        text=f"{group1} →",
        showarrow=False,
        font=dict(size=12, color=color_up, family="Arial, sans-serif"),
        xanchor="right",
    )

    up_n = (plot_df["class"] == "up").sum()
    dn_n = (plot_df["class"] == "down").sum()

    fig.update_layout(
        title=dict(
            text=(f"<b>{group1}</b> vs <b>{group2}</b>"
                  f"  —  {up_n} up in {group1} · {dn_n} up in {group2}"),
            font=_FONT,
        ),
        xaxis_title="log₂ Fold Change",
        yaxis_title="-log₁₀(adjusted p-value)",
        legend=dict(
            title=dict(text="<b><u>Legend</u></b>", font=_FONT, side="top center"),
            x=1.02,
            y=1.0,
            xanchor="left",
            yanchor="top",
            bgcolor="rgba(255,255,255,0.9)",
            bordercolor="#cccccc",
            borderwidth=1,
            font=_FONT,
        ),
        plot_bgcolor="#f5f5f5",
        paper_bgcolor="#ffffff",
        font=_FONT,
        width=850,
        height=700,
        margin=dict(l=70, r=180, t=80, b=70),
        hovermode="closest",
        xaxis=dict(range=[-x_range, x_range], showgrid=True,
                   gridcolor="#dddddd", zeroline=True,
                   zerolinecolor="#aaaaaa", zerolinewidth=1,
                   tickfont=_FONT,
                   title_font=_FONT),
        yaxis=dict(range=[0, y_max], showgrid=True, gridcolor="#dddddd",
                   tickfont=_FONT,
                   title_font=_FONT),
    )
    return fig
