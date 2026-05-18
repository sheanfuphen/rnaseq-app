# RNAseq Analysis & Visualization App

A user-friendly Streamlit application for RNA sequencing analysis and publication-ready plot generation — no coding required.

---

## Overview

This app is designed for bench scientists and bioinformaticians who want to perform differential gene expression analysis and generate high-quality visualizations from quantified RNA-seq count data without writing a single line of code.

**Key Features:**
- Upload a single combined TSV containing all samples — the app reads the `srr_id` column to identify samples automatically
- Assign samples to named, color-coded groups through a point-and-click interface
- Differential Gene Expression (DGE) analysis between every pair of groups using DESeq2
- Interactive Volcano Plots with customizable colors and square dimensions
- 2D and 3D PCA plots with per-group hex color customization
- Gene list-driven Heatmaps with step-by-step gene reordering
- GO Term Enrichment analysis via the Enrichr API (bar chart and dot plot)
- Export any plot as a high-resolution PNG or SVG

---

## Input Data Format

The app accepts a **single combined TSV file** containing all samples. Each row is one transcript or gene entry. The file can have any number of columns — only these four are required:

| Column | Description |
|---|---|
| `est_counts` | Estimated raw counts (used for DGE) |
| `tpm` | Transcripts per million (used for PCA and heatmap) |
| `gene_name` | Gene symbol (e.g. `ACTB`, `IL6`) |
| `srr_id` | Sample identifier — distinguishes samples within the file |

All other columns (e.g. `target_id`, `length`, `eff_length`) are ignored automatically.

If multiple rows share the same `gene_name` and `srr_id` (i.e. multiple transcripts per gene), `est_counts` and `tpm` are summed to produce gene-level values.

This format is compatible with kallisto output that has been post-processed to add a `gene_name` and `srr_id` column, as well as any other quantification tool that can produce this schema.

---

## Getting Started

### 1. Installation

**Prerequisites:** Python 3.9+

```bash
git clone https://github.com/am15443/rnaseq-app.git
cd rnaseq-app
pip install -r requirements.txt
```

### 2. Configure Upload Limit

For large TSV files (>200MB), increase Streamlit's default upload limit by creating `.streamlit/config.toml`:

```toml
[server]
maxUploadSize = 2000
```

This is already included in the repository.

### 3. Run the App

```bash
streamlit run app.py
```

The app will open in your browser at `http://localhost:8501`.

---

## Usage Guide

### Step 1 — Upload Combined TSV

Upload your single combined TSV file using the file uploader in the sidebar. The app will automatically detect all unique `srr_id` values and display them as available samples.

### Step 2 — Create Sample Groups

- Click **"Add Group"** to create a new group
- Give each group a descriptive name (e.g. `Control`, `Treatment`)
- Pick a color using the color picker — this color is shared across PCA, heatmap annotation bar, and serves as the default for volcano plots
- Use the dropdown to assign `srr_id` values to each group
- Create as many groups as needed — there is no limit
- The **Samples tab** shows a table of all detected samples and their current group assignments

### Step 3 — Run DGE Analysis

Click **"Run DGE Analysis"** in the sidebar. DGE is computed between every pair of groups automatically. For example, groups A, B, and C will produce comparisons A vs B, A vs C, and B vs C.

Once complete, all analysis tabs become available.

---

## Analysis Tabs

### 📊 DGE Results

- One results table per pairwise comparison
- Columns: gene, baseMean, log2FoldChange, pvalue, padj
- Summary metrics shown: total genes tested, significant genes (padj < 0.05), upregulated count
- Download results as CSV per comparison

### 🌋 Volcano Plots

- One square volcano plot per group pair
- X axis: log₂ fold change — positive = higher in group 1
- Y axis: −log₁₀(adjusted p-value)
- Group names labeled at top-left and top-right corners
- Adjustable log₂FC and padj significance thresholds
- **Per-comparison hex color inputs** — override the default group colors with any hex code
- Export each plot as PNG or SVG

### 🔵 PCA

- One point per sample, colored by group
- **Per-group hex color inputs** above the plots — override sidebar colors with any hex code
- Both **2D and 3D PCA** generated from log₂(TPM + 1) values
- Adjustable plot width and height
- Export each plot as PNG or SVG

### 🔥 Heatmap

- Select genes from one or more **gene list CSVs** placed in the `gene_lists/` folder
- Add additional genes manually by typing comma-separated names (case-insensitive)
- **Gene reordering UI** — select any gene from the dropdown and move it up, down, or to the top; the current order is shown in a numbered list before generating
- Heatmap layout:
  - Y axis: genes (in your specified order)
  - X axis: individual samples (one column per srr_id), ordered by group
  - Colored annotation bar below X axis — one color per group, matching PCA colors
- Values shown as Z-scored log₂(TPM + 1)
- Hover tooltip shows gene, sample, group, raw TPM, and Z-score
- Export as PNG or SVG

### 🧬 GO Enrichment

- Select a pairwise comparison and gene direction (up in group 1, up in group 2, or all significant)
- Adjust log₂FC and padj thresholds independently
- Choose a gene set library:
  - GO Biological Process
  - GO Molecular Function
  - GO Cellular Component
  - KEGG Pathways
  - Reactome
  - MSigDB Hallmarks
- Choose **Bar chart** or **Dot plot** visualization
- Override plot color with any hex code
- Results queried live from the **Enrichr API** via `gseapy` — no local annotation files needed
- Download enrichment results as CSV
- Export plot as PNG or SVG

---

## Gene List CSVs

Place CSV files in the `gene_lists/` directory at the root of the project. Each CSV should have at least one column containing gene names or symbols. The column header can be anything — the app reads the first column.

```
gene_lists/
├── cytokines.csv
├── cell_cycle_genes.csv
└── custom_pathway.csv
```

Example format:
```
gene_name
IL6
TNF
CXCL10
```

The app detects all CSVs in this folder automatically and lists them by filename in the heatmap dropdown. If the same gene appears in multiple selected CSVs it is deduplicated.

---

## Project Structure

```
rnaseq-app/
├── app.py                   # Main Streamlit application entry point
├── requirements.txt         # Python dependencies
├── .streamlit/
│   └── config.toml          # Upload size and server config
├── gene_lists/              # Place gene list CSVs here
│   └── .gitkeep
└── src/
    ├── __init__.py
    ├── data_loader.py       # TSV parsing, column validation, matrix building
    ├── dge.py               # DESeq2-based DGE (pydeseq2) with t-test fallback
    ├── pca.py               # PCA computation and 2D/3D Plotly figures
    ├── volcano.py           # Volcano plot generation
    ├── heatmap.py           # Heatmap with group annotation bar
    └── go_enrichment.py     # GO enrichment via Enrichr API (gseapy)
```

---

## Dependencies

| Package | Purpose |
|---|---|
| `streamlit` | Web app framework |
| `pandas` | Data loading and manipulation |
| `numpy` | Numerical computation |
| `pydeseq2` | DESeq2-based differential gene expression |
| `plotly` | Interactive volcano, PCA, heatmap, and GO plots |
| `scikit-learn` | PCA computation |
| `scipy` | Statistical testing (fallback DGE) |
| `statsmodels` | Benjamini-Hochberg FDR correction |
| `gseapy` | GO enrichment via Enrichr API |
| `kaleido` | PNG and SVG export for Plotly figures |

---

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/your-feature`)
3. Commit your changes (`git commit -m 'Add your feature'`)
4. Push to the branch (`git push origin feature/your-feature`)
5. Open a Pull Request
