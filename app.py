"""
app.py
──────
RNAseq Analysis Suite — Streamlit entry point.

Single combined TSV input (all samples in one file):
    target_id | length | eff_length | est_counts | tpm | gene_name | srr_id

Groups and sample assignments are stored entirely in session_state so that
widget interactions (typing, clicking) never reset other widgets.
"""

import streamlit as st
import pandas as pd
import io
import re
from pathlib import Path

st.set_page_config(
    page_title="RNAseq Analysis Suite",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── Custom CSS ────────────────────────────────────────────────────────────────
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=DM+Mono:wght@400;500&family=DM+Sans:wght@300;400;500;600&display=swap');

/* ── Base font size boost (Streamlit default is ~14px, we want ~18px) ── */
html { font-size: 18px !important; }

/* ── Dark backgrounds ── */
.stApp { background-color: #0f1117; }
section[data-testid="stSidebar"] { background-color: #1a1d27; }
section[data-testid="stSidebar"] > div { background-color: #1a1d27; }

/* ── White text globally ── */
.stApp, .stApp * {
    color: #ffffff;
    font-family: 'DM Sans', sans-serif;
}

/* ── Sidebar text ── */
section[data-testid="stSidebar"] * { color: #ffffff; }

/* ── Inputs dark ── */
input, textarea, [data-baseweb="input"] input {
    background-color: #2a2d3a !important;
    color: #ffffff !important;
}

/* ── Multiselect dropdown dark ── */
[data-baseweb="select"] > div { background-color: #2a2d3a !important; }
[data-baseweb="menu"] { background-color: #2a2d3a !important; }
[data-baseweb="option"] { background-color: #2a2d3a !important; color: #ffffff !important; }
[data-baseweb="tag"] { background-color: #4361ee !important; color: #ffffff !important; }

/* ── Buttons ── */
.stButton > button {
    background-color: #2a2d3a;
    color: #ffffff;
    border: 1px solid #444;
    font-size: 1rem;
}

/* ── Tab text ── */
button[data-baseweb="tab"] { color: #ffffff !important; }

.main-title {
    font-family: 'DM Mono', monospace;
    font-size: 2.2rem;
    font-weight: 500;
    color: #ffffff;
    letter-spacing: -0.02em;
    margin-bottom: 0;
}
.subtitle { color: #aaaaaa; font-size: 1rem; margin-top: 0.1rem; margin-bottom: 0; }
.schema-box {
    background: #1e2130;
    border-left: 3px solid #4361ee;
    border-radius: 6px;
    padding: 0.6rem 1rem;
    font-family: 'DM Mono', monospace;
    font-size: 0.9rem;
    color: #dddddd;
    margin: 0.6rem 0 1rem 0;
    line-height: 1.8;
}
/* ── Fix duplicate upload button text ── */
[data-testid="stIconMaterial"] {
    display: none !important;
}

.group-card {
    background: #1e2130;
    border: 1.5px solid #333a55;
    border-radius: 10px;
    padding: 0.85rem 1rem;
    margin-bottom: 0.75rem;
}
</style>
""", unsafe_allow_html=True)

# ── Src imports ───────────────────────────────────────────────────────────────
from src.data_loader  import load_combined_tsv, validate_counts
from src.dge          import run_dge_all_pairs
from src.pca          import compute_pca, plot_pca_2d, plot_pca_3d
from src.volcano      import plot_volcano
from src.heatmap      import build_heatmap, order_gene_list_for_heatmap
from src.go_enrichment import run_enrichment, plot_go_bars, plot_go_dots, GO_LIBRARIES


def export_buttons(fig, filename_stem: str, include_pdf: bool = False):
    """Render PNG and SVG download buttons for a Plotly figure."""
    try:
        png_bytes = fig.to_image(format="png", scale=3)
        st.download_button(
            "⬇ Download PNG",
            data=png_bytes,
            file_name=f"{filename_stem}.png",
            mime="image/png",
            key=f"dl_png_{filename_stem}",
        )
    except Exception:
        st.caption("⚠️ PNG export unavailable — run `pip install kaleido` to enable.")

    try:
        svg_bytes = fig.to_image(format="svg")
        st.download_button(
            "⬇ Download SVG",
            data=svg_bytes,
            file_name=f"{filename_stem}.svg",
            mime="image/svg+xml",
            key=f"dl_svg_{filename_stem}",
        )
    except Exception:
        pass

    if include_pdf:
        try:
            pdf_bytes = fig.to_image(format="pdf", scale=3)
            st.download_button(
                "⬇ Download PDF",
                data=pdf_bytes,
                file_name=f"{filename_stem}.pdf",
                mime="application/pdf",
                key=f"dl_pdf_{filename_stem}",
            )
        except Exception:
            st.caption("⚠️ PDF export unavailable — run `pip install kaleido` to enable.")

# ── Session state defaults ────────────────────────────────────────────────────
def _init(key, val):
    if key not in st.session_state:
        st.session_state[key] = val

_init("gene_sample_df", None)
_init("all_srr_ids",    [])
_init("group_ids",      [])          # ordered list of stable group ids (ints)
_init("group_names",    {})          # {gid: str}
_init("group_colors",   {})          # {gid: hex str}
_init("group_samples",  {})          # {gid: [srr_id, ...]}
_init("next_gid",       0)
_init("dge_results",    {})

DEFAULT_COLORS = [
    "#4361ee", "#f72585", "#4cc9f0", "#7209b7",
    "#06d6a0", "#f77f00", "#ef233c", "#3a86ff",
]

# ── Helper: add / remove groups ───────────────────────────────────────────────
def add_group():
    gid = st.session_state.next_gid
    st.session_state.next_gid += 1
    idx = len(st.session_state.group_ids)
    st.session_state.group_ids.append(gid)
    st.session_state.group_names[gid]   = f"Group_{idx + 1}"
    st.session_state.group_colors[gid]  = DEFAULT_COLORS[idx % len(DEFAULT_COLORS)]
    st.session_state.group_samples[gid] = []

def remove_group(gid):
    st.session_state.group_ids.remove(gid)
    for d in (st.session_state.group_names,
              st.session_state.group_colors,
              st.session_state.group_samples):
        d.pop(gid, None)

# ── Callbacks that write directly into session_state ─────────────────────────
def on_name_change(gid):
    st.session_state.group_names[gid] = st.session_state[f"name_{gid}"]

def on_color_change(gid):
    st.session_state.group_colors[gid] = st.session_state[f"color_{gid}"]

def on_samples_change(gid):
    st.session_state.group_samples[gid] = st.session_state[f"samples_{gid}"]

# ── Header ────────────────────────────────────────────────────────────────────
st.markdown('<p class="main-title">🧬 RNAseq Analysis Suite</p>', unsafe_allow_html=True)
st.markdown(
    '<p class="subtitle">Upload your combined TSV · Assign srr_ids to groups · '
    'Explore differential expression</p>',
    unsafe_allow_html=True,
)
st.markdown("""
<div class="schema-box">
Required TSV columns (all other columns are ignored):<br>
<b>est_counts</b> &nbsp;·&nbsp; <b>tpm</b>
&nbsp;·&nbsp; <b>gene_name</b> &nbsp;·&nbsp; <b>srr_id</b>
</div>
""", unsafe_allow_html=True)
st.divider()

# ═════════════════════════════════════════════════════════════════════════════
# SIDEBAR
# ═════════════════════════════════════════════════════════════════════════════
with st.sidebar:

    # ── Step 1: Upload ────────────────────────────────────────────────────────
    st.markdown("### 📂 Step 1 — Upload Combined TSV")
    st.caption("One file containing all samples. The srr_id column identifies each sample.")

    uploaded = st.file_uploader(
        "Upload combined TSV",
        type=["tsv", "txt"],
        accept_multiple_files=False,
        label_visibility="collapsed",
    )

    if uploaded is not None and st.session_state.gene_sample_df is None:
        # Only parse once — if user re-uploads a new file they must refresh
        try:
            gene_sample_df, srr_ids = load_combined_tsv(uploaded)
            st.session_state.gene_sample_df = gene_sample_df
            st.session_state.all_srr_ids    = srr_ids
            st.success(f"✅ {len(srr_ids)} sample(s) detected.")
        except Exception as e:
            st.error(f"Failed to load file: {e}")

    if st.session_state.gene_sample_df is not None:
        st.caption(f"Samples: {', '.join(st.session_state.all_srr_ids)}")

        # Button to clear and re-upload
        if st.button("↩ Clear & upload new file", use_container_width=True):
            for k in ["gene_sample_df", "all_srr_ids", "dge_results",
                      "_counts_df", "_tpm_df", "_sample_meta"]:
                st.session_state.pop(k, None)
            st.session_state.gene_sample_df = None
            st.session_state.all_srr_ids    = []
            st.session_state.group_ids      = []
            st.session_state.group_names    = {}
            st.session_state.group_colors   = {}
            st.session_state.group_samples  = {}
            st.session_state.next_gid       = 0
            st.session_state.dge_results    = {}
            st.rerun()

    st.divider()

    # ── Step 2: Group Builder ─────────────────────────────────────────────────
    st.markdown("### 🗂 Step 2 — Define Groups")
    st.caption("Create groups, name them, pick a colour, then assign srr_ids.")

    all_srr_ids = st.session_state.all_srr_ids

    if not all_srr_ids:
        st.caption("⬆️ Upload a TSV file first.")
    else:
        st.button("➕ Add Group", on_click=add_group, use_container_width=True)

        for gid in list(st.session_state.group_ids):
            st.markdown('<div class="group-card">', unsafe_allow_html=True)

            col_name, col_color = st.columns([3, 1])
            with col_name:
                st.text_input(
                    "Group name",
                    value=st.session_state.group_names[gid],
                    key=f"name_{gid}",
                    label_visibility="collapsed",
                    on_change=on_name_change,
                    args=(gid,),
                )
            with col_color:
                st.color_picker(
                    "Color",
                    value=st.session_state.group_colors[gid],
                    key=f"color_{gid}",
                    label_visibility="collapsed",
                    on_change=on_color_change,
                    args=(gid,),
                )

            st.multiselect(
                "Assign srr_ids",
                options=all_srr_ids,
                default=st.session_state.group_samples[gid],
                key=f"samples_{gid}",
                label_visibility="collapsed",
                placeholder="Select srr_ids for this group…",
                on_change=on_samples_change,
                args=(gid,),
            )

            st.button(
                "🗑 Remove group",
                key=f"del_{gid}",
                on_click=remove_group,
                args=(gid,),
                use_container_width=True,
            )

            st.markdown("</div>", unsafe_allow_html=True)

    st.divider()

    # ── Step 3: Run Analysis ──────────────────────────────────────────────────
    st.markdown("### 🚀 Step 3 — Run Analysis")

    # Build groups dict {name: [srr_ids]} for downstream modules
    groups_dict = {
        st.session_state.group_names[gid]: st.session_state.group_samples[gid]
        for gid in st.session_state.group_ids
    }
    colors_dict = {
        st.session_state.group_names[gid]: st.session_state.group_colors[gid]
        for gid in st.session_state.group_ids
    }
    valid_groups = {k: v for k, v in groups_dict.items() if v}
    can_run = st.session_state.gene_sample_df is not None and len(valid_groups) >= 2

    if st.button("Run DGE Analysis", type="primary",
                 use_container_width=True, disabled=not can_run):
        with st.spinner("Building matrices and running DGE…"):
            try:
                counts_df, tpm_df, sample_meta = validate_counts(
                    st.session_state.gene_sample_df, groups_dict,
                )
                results = run_dge_all_pairs(counts_df, sample_meta)
                st.session_state.dge_results        = results
                st.session_state["_counts_df"]      = counts_df
                st.session_state["_tpm_df"]         = tpm_df
                st.session_state["_sample_meta"]    = sample_meta
                st.session_state["_groups_dict"]    = groups_dict
                st.session_state["_colors_dict"]    = colors_dict
                st.success(f"Done — {len(results)} comparison(s) ready.")
            except Exception as e:
                st.error(f"Analysis error: {e}")

    if not can_run:
        if st.session_state.gene_sample_df is None:
            st.caption("⬆️ Upload a TSV file to get started.")
        else:
            st.caption("⚠️ Define ≥ 2 groups with srr_ids assigned.")

# ═════════════════════════════════════════════════════════════════════════════
# MAIN AREA — Tabs
# ═════════════════════════════════════════════════════════════════════════════
tab_samples, tab_dge, tab_volcano, tab_pca, tab_heatmap, tab_go = st.tabs([
    "🧪 Samples", "📊 DGE Results", "🌋 Volcano Plots", "🔵 PCA", "🔥 Heatmap", "🧬 GO Enrichment",
])

# ── Samples overview ──────────────────────────────────────────────────────────
with tab_samples:
    if st.session_state.gene_sample_df is None:
        st.info("Upload a combined TSV file using the sidebar to get started.")
    else:
        st.markdown("#### Samples detected in uploaded file")
        group_map = {
            s: st.session_state.group_names[gid]
            for gid in st.session_state.group_ids
            for s in st.session_state.group_samples[gid]
        }
        summary = pd.DataFrame({
            "srr_id": st.session_state.all_srr_ids,
            "assigned_group": [
                group_map.get(s, "— unassigned —")
                for s in st.session_state.all_srr_ids
            ],
        })
        st.dataframe(summary, use_container_width=True, hide_index=True)
        unassigned = [s for s in st.session_state.all_srr_ids if s not in group_map]
        if unassigned:
            st.warning(f"{len(unassigned)} sample(s) not yet assigned: {', '.join(unassigned)}")
        else:
            st.success("All samples are assigned to a group. ✅")

# ── DGE Results ───────────────────────────────────────────────────────────────
with tab_dge:
    if not st.session_state.dge_results:
        st.info("Run DGE analysis using the sidebar to see results here.")
    else:
        for (g1, g2), df in st.session_state.dge_results.items():
            with st.expander(f"**{g1}** vs **{g2}**", expanded=True):
                c1, c2, c3 = st.columns(3)
                c1.metric("Genes tested", len(df))
                c2.metric("Significant (padj < 0.05)", int((df["padj"] < 0.05).sum()))
                c3.metric(f"Up in {g1} (log2FC > 1)",
                          int(((df["padj"] < 0.05) & (df["log2FoldChange"] > 1)).sum()))
                st.dataframe(df.sort_values("padj").round(5),
                             use_container_width=True, height=320)
                st.download_button(
                    f"⬇ Download {g1}_vs_{g2}.csv",
                    data=df.to_csv(index=True).encode(),
                    file_name=f"DGE_{g1}_vs_{g2}.csv",
                    mime="text/csv",
                    key=f"dl_{g1}_{g2}",
                )

# ── Volcano Plots ─────────────────────────────────────────────────────────────
with tab_volcano:
    if not st.session_state.dge_results:
        st.info("Run DGE analysis first.")
    else:
        colors_dict = st.session_state.get("_colors_dict", {})

        col_fc, col_pv = st.columns(2)
        fc_thresh = col_fc.slider("log₂FC threshold", 0.5, 4.0, 1.0, 0.25)
        pv_thresh = col_pv.slider("-log₁₀(padj) threshold", 1.0, 10.0, 1.301, 0.1,
                                  help="1.301 ≈ padj < 0.05")

        st.markdown("**Dot colours** — defaults to your group colours; override with a hex code below:")

        for (g1, g2), df in st.session_state.dge_results.items():
            st.markdown(f"#### {g1} vs {g2}")

            cc1, cc2 = st.columns(2)
            default_up   = colors_dict.get(g1, "#4361ee")
            default_down = colors_dict.get(g2, "#f72585")
            color_up   = cc1.text_input(f"Color for {g1} (hex)", value=default_up,
                                        key=f"vcol_up_{g1}_{g2}")
            color_down = cc2.text_input(f"Color for {g2} (hex)", value=default_down,
                                        key=f"vcol_dn_{g1}_{g2}")

            # Validate hex — fall back to default if invalid
            if not re.match(r"^#[0-9a-fA-F]{6}$", color_up):
                color_up = default_up
                cc1.caption("⚠️ Invalid hex, using default.")
            if not re.match(r"^#[0-9a-fA-F]{6}$", color_down):
                color_down = default_down
                cc2.caption("⚠️ Invalid hex, using default.")

            fig = plot_volcano(df, g1, g2,
                               fc_thresh=fc_thresh,
                               neg_log10_padj_thresh=pv_thresh,
                               color_up=color_up,
                               color_down=color_down)
            st.plotly_chart(fig, use_container_width=False)
            export_buttons(fig, f"volcano_{g1}_vs_{g2}")

# ── PCA ───────────────────────────────────────────────────────────────────────
with tab_pca:
    tpm_ready = "_tpm_df" in st.session_state
    if not tpm_ready:
        valid_for_pca = {
            st.session_state.group_names[gid]: st.session_state.group_samples[gid]
            for gid in st.session_state.group_ids
            if st.session_state.group_samples[gid]
        }
        if st.session_state.gene_sample_df is not None and len(valid_for_pca) >= 2:
            if st.button("▶ Compute PCA"):
                try:
                    _, tpm_df, sample_meta = validate_counts(
                        st.session_state.gene_sample_df, valid_for_pca,
                    )
                    st.session_state["_tpm_df"]      = tpm_df
                    st.session_state["_sample_meta"] = sample_meta
                    st.session_state["_colors_dict"] = {
                        st.session_state.group_names[gid]: st.session_state.group_colors[gid]
                        for gid in st.session_state.group_ids
                    }
                    st.rerun()
                except Exception as e:
                    st.error(f"PCA error: {e}")
        else:
            st.info("Upload a file and define ≥ 2 groups to compute PCA.")

    if tpm_ready:
        tpm_df      = st.session_state["_tpm_df"]
        sample_meta = st.session_state["_sample_meta"]
        colors_dict = st.session_state.get("_colors_dict", {}).copy()

        col_w, col_h = st.columns(2)
        pca_w = col_w.slider("Plot width (px)",  400, 1400, 850, 50)
        pca_h = col_h.slider("Plot height (px)", 300, 900,  550, 50)

        st.markdown("**Group colours** — defaults to sidebar pickers; override with a hex code:")
        color_cols = st.columns(max(1, len(colors_dict)))
        pca_colors = {}
        for i, (gname, default_color) in enumerate(colors_dict.items()):
            with color_cols[i % len(color_cols)]:
                hex_input = st.text_input(f"{gname} color (hex)", value=default_color,
                                          key=f"pca_color_{gname}")
                if re.match(r"^#[0-9a-fA-F]{6}$", hex_input):
                    pca_colors[gname] = hex_input
                else:
                    pca_colors[gname] = default_color
                    if hex_input != default_color:
                        st.caption("⚠️ Invalid hex, using default.")

        try:
            coords_2d, coords_3d, explained = compute_pca(tpm_df)
            st.markdown("#### 2D PCA")
            fig2d = plot_pca_2d(coords_2d, sample_meta, explained, pca_colors,
                                width=pca_w, height=pca_h)
            st.plotly_chart(fig2d, use_container_width=False)
            export_buttons(fig2d, "PCA_2D")

            st.markdown("#### 3D PCA")
            fig3d = plot_pca_3d(coords_3d, sample_meta, explained, pca_colors,
                                width=pca_w, height=pca_h)
            st.plotly_chart(fig3d, use_container_width=False)
            export_buttons(fig3d, "PCA_3D")
        except Exception as e:
            st.error(f"PCA failed: {e}")

# ── Heatmap ───────────────────────────────────────────────────────────────────
with tab_heatmap:
    if "_tpm_df" not in st.session_state:
        st.info("Run DGE analysis (or compute PCA) first so TPM data is available.")
    else:
        tpm_df      = st.session_state["_tpm_df"]
        sample_meta = st.session_state["_sample_meta"]
        colors_dict = st.session_state.get("_colors_dict", {})
        groups_dict = st.session_state.get("_groups_dict", {})

        # ── Gene list CSVs ────────────────────────────────────────────────────
        gene_list_dir = Path("gene_lists")
        gene_list_dir.mkdir(exist_ok=True)
        csv_files = sorted(gene_list_dir.glob("*.csv"))
        selected_genes: set = set()

        if csv_files:
            chosen_csvs = st.multiselect(
                "📋 Select gene list CSV(s)",
                options=[f.stem for f in csv_files],
                help="Place CSV files in the gene_lists/ folder. First column = gene names.",
            )
            for csv_name in chosen_csvs:
                try:
                    gdf = pd.read_csv(gene_list_dir / f"{csv_name}.csv")
                    col = gdf.columns[0]
                    selected_genes.update(
                        gdf[col].dropna().astype(str).str.strip().str.upper().tolist()
                    )
                except Exception as e:
                    st.warning(f"Could not read {csv_name}.csv: {e}")
        else:
            st.caption("No CSVs found in gene_lists/. Add gene list CSVs there to use the dropdown.")

        manual = st.text_input(
            "✏️ Additional genes (comma-separated, case-insensitive)",
            placeholder="e.g.  ACTB, GAPDH, TP53",
        )
        if manual:
            for g in manual.split(","):
                g = g.strip().upper()
                if g:
                    selected_genes.add(g)

        if selected_genes:
            st.caption(f"**{len(selected_genes)}** unique gene(s) selected.")

        # ── Gene row ordering by group contrast ───────────────────────────────
        group_names = [k for k, v in groups_dict.items() if v]
        contrast_pair = None
        if len(group_names) >= 2:
            st.markdown("**📐 Gene row order** — cluster by expression contrast between two groups:")
            c1, c2 = st.columns(2)
            with c1:
                hm_group1 = st.selectbox(
                    "Higher at top of heatmap",
                    options=group_names,
                    index=0,
                    key="hm_contrast_g1",
                )
            with c2:
                other = [g for g in group_names if g != hm_group1]
                hm_group2 = st.selectbox(
                    "Higher at bottom of heatmap",
                    options=other if other else group_names,
                    index=0,
                    key="hm_contrast_g2",
                )
            if hm_group1 != hm_group2:
                contrast_pair = (hm_group1, hm_group2)
                st.caption(
                    f"Genes are sorted by mean log₂(TPM+1): enriched in **{hm_group1}** at the top, "
                    f"enriched in **{hm_group2}** at the bottom, and similar between groups in the middle."
                )
            else:
                st.warning("Pick two different groups for contrast ordering.")

        contrast_key = contrast_pair
        gene_selection_key = (tuple(sorted(selected_genes)), contrast_key)
        if st.session_state.get("_heatmap_gene_key") != gene_selection_key:
            st.session_state["_heatmap_gene_key"] = gene_selection_key
            gene_list = list(selected_genes)
            if contrast_pair:
                st.session_state["_heatmap_gene_order"] = order_gene_list_for_heatmap(
                    tpm_df, groups_dict, gene_list, contrast_pair,
                )
            else:
                st.session_state["_heatmap_gene_order"] = gene_list

        # ── Gene reordering UI ────────────────────────────────────────────────
        if selected_genes:
            st.markdown("**🔀 Manual reorder** (optional) — select a gene and move it up or down:")

            gene_order = st.session_state["_heatmap_gene_order"]

            # Show current order as a small preview — above the buttons
            with st.expander("📋 Current gene order (top → bottom on heatmap)", expanded=False):
                for i, g in enumerate(st.session_state["_heatmap_gene_order"]):
                    st.write(f"{i+1}. {g}")

            col_list, col_btns = st.columns([3, 1])

            with col_list:
                selected_gene = st.selectbox(
                    "Select gene to move",
                    options=gene_order,
                    key="hm_gene_select",
                    label_visibility="collapsed",
                )

            with col_btns:
                b1, b2, b3 = st.columns(3)
                move_top  = b1.button("⏫", key="hm_top",  help="Move to top")
                move_up   = b2.button("⬆",  key="hm_up",   help="Move up one")
                move_down = b3.button("⬇",  key="hm_down", help="Move down one")

            if selected_gene and selected_gene in gene_order:
                idx = gene_order.index(selected_gene)
                if move_top and idx > 0:
                    gene_order.insert(0, gene_order.pop(idx))
                    st.session_state["_heatmap_gene_order"] = gene_order
                    st.rerun()
                if move_up and idx > 0:
                    gene_order[idx], gene_order[idx-1] = gene_order[idx-1], gene_order[idx]
                    st.session_state["_heatmap_gene_order"] = gene_order
                    st.rerun()
                if move_down and idx < len(gene_order) - 1:
                    gene_order[idx], gene_order[idx+1] = gene_order[idx+1], gene_order[idx]
                    st.session_state["_heatmap_gene_order"] = gene_order
                    st.rerun()

        hm_use_tpm = st.toggle(
            "Show raw TPM (instead of per-gene Z-scores on log₂ TPM)",
            key="hm_value_toggle",
            help="Z-score mode scales each gene row; TPM mode colours cells by raw TPM.",
        )
        hm_value_mode = "tpm" if hm_use_tpm else "zscore"

        if st.button("🔥 Generate Heatmap", type="primary",
                     disabled=len(selected_genes) == 0):
            try:
                gene_order = st.session_state.get(
                    "_heatmap_gene_order", list(selected_genes),
                )
                fig = build_heatmap(
                    tpm_df, sample_meta, gene_order,
                    colors_dict, groups_dict,
                    value_mode=hm_value_mode,
                )
                st.session_state["_heatmap_fig"] = fig
                st.session_state["_heatmap_value_mode"] = hm_value_mode
            except Exception as e:
                st.error(f"Heatmap error: {e}")

        # Keep the figure visible after reordering without re-clicking generate
        if "_heatmap_fig" in st.session_state and selected_genes:
            if st.session_state.get("_heatmap_value_mode") != hm_value_mode:
                st.info("Value mode changed — click **Generate Heatmap** to update the plot.")
            fig = st.session_state["_heatmap_fig"]
            st.plotly_chart(fig, use_container_width=True)
            export_buttons(fig, "heatmap", include_pdf=True)


# ── GO Enrichment ─────────────────────────────────────────────────────────────
with tab_go:
    if not st.session_state.dge_results:
        st.info("Run DGE analysis first to enable GO enrichment.")
    else:
        st.markdown(
            "Select a comparison and a gene set to query Enrichr for enriched GO terms. "
            "Genes are filtered by your significance thresholds before submission."
        )

        # ── Settings ──────────────────────────────────────────────────────────
        col_a, col_b, col_c = st.columns(3)

        comparison_options = [f"{g1} vs {g2}" for (g1, g2) in st.session_state.dge_results]
        chosen_comp = col_a.selectbox("Comparison", comparison_options, key="go_comp")

        library_label = col_b.selectbox("Gene set library", list(GO_LIBRARIES.keys()),
                                        key="go_library")
        library = GO_LIBRARIES[library_label]

        direction = col_c.selectbox(
            "Gene direction",
            ["Up in group 1", "Up in group 2", "All significant"],
            key="go_direction",
        )

        col_d, col_e, col_f = st.columns(3)
        go_fc_thresh  = col_d.slider("log₂FC threshold", 0.5, 4.0, 1.0, 0.25, key="go_fc")
        go_pv_thresh  = col_e.slider("padj threshold", 0.001, 0.1, 0.05, 0.001,
                                     format="%.3f", key="go_pv")
        go_top_n      = col_f.slider("Top N terms", 5, 50, 20, 5, key="go_topn")

        plot_type = st.radio("Plot type", ["Bar chart", "Dot plot"],
                             horizontal=True, key="go_plottype")

        # Hex color override
        go_color_input = st.text_input("Plot color (hex)", value="#4361ee", key="go_color")
        go_color = go_color_input if re.match(r"^#[0-9a-fA-F]{6}$", go_color_input) else "#4361ee"
        if go_color_input != go_color:
            st.caption("⚠️ Invalid hex, using default #4361ee.")

        if st.button("🧬 Run GO Enrichment", type="primary", key="go_run"):
            # Parse chosen comparison
            g1, g2 = chosen_comp.split(" vs ", 1)
            dge_df = st.session_state.dge_results.get((g1, g2))
            if dge_df is None:
                st.error(f"Could not find results for {chosen_comp}.")
            else:
                sig = dge_df[dge_df["padj"] <= go_pv_thresh].copy()

                if direction == "Up in group 1":
                    sig = sig[sig["log2FoldChange"] >= go_fc_thresh]
                elif direction == "Up in group 2":
                    sig = sig[sig["log2FoldChange"] <= -go_fc_thresh]
                else:
                    sig = sig[sig["log2FoldChange"].abs() >= go_fc_thresh]

                gene_list = sig.index.tolist()
                st.caption(f"Submitting **{len(gene_list)}** genes to Enrichr…")

                if len(gene_list) == 0:
                    st.warning("No genes pass the selected thresholds. "
                               "Try relaxing the FC or padj cutoffs.")
                else:
                    with st.spinner("Querying Enrichr…"):
                        results_df = run_enrichment(
                            gene_list,
                            library=library,
                            cutoff=go_pv_thresh,
                            top_n=go_top_n,
                        )

                    if results_df is not None and not results_df.empty:
                        title = (f"{library_label} — {direction}<br>"
                                 f"{g1} vs {g2}")

                        if plot_type == "Bar chart":
                            fig = plot_go_bars(results_df, title=title,
                                               color=go_color, width=900)
                        else:
                            fig = plot_go_dots(results_df, title=title,
                                               color=go_color, width=900)

                        st.plotly_chart(fig, use_container_width=False)
                        export_buttons(fig, f"GO_{g1}_vs_{g2}_{library_label.replace(' ','_')}")

                        st.download_button(
                            "⬇ Download GO results CSV",
                            data=results_df.to_csv(index=False).encode(),
                            file_name=f"GO_{g1}_vs_{g2}_{library_label.replace(' ','_')}.csv",
                            mime="text/csv",
                            key="go_dl",
                        )
