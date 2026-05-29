"""
report.py — Automated HTML report generation for RNA-seq analysis results.

Generates a self-contained HTML report with:
  - Summary statistics
  - Interactive DEA results table
  - Embedded plots (volcano, MA, PCA, heatmap)
  - Enrichment results
"""

import os
import base64
import datetime
import pandas as pd
import numpy as np
from typing import Optional


def _img_to_base64(path: str) -> str:
    """Convert image file to base64 string for embedding in HTML."""
    if not os.path.exists(path):
        return ""
    with open(path, "rb") as f:
        data = base64.b64encode(f.read()).decode("utf-8")
    ext = os.path.splitext(path)[1].lstrip(".")
    return f"data:image/{ext};base64,{data}"


def _df_to_html_table(df: pd.DataFrame, max_rows: int = 200) -> str:
    """Convert DataFrame to styled HTML table."""
    df_show = df.head(max_rows).copy()
    # Round floats
    for col in df_show.select_dtypes(include=[float]).columns:
        df_show[col] = df_show[col].apply(
            lambda x: f"{x:.4e}" if abs(x) < 0.001 and x != 0 else f"{x:.4f}"
            if pd.notna(x) else "NA"
        )
    return df_show.to_html(
        classes="table table-striped table-hover table-sm",
        border=0,
        escape=True,
    )


def generate_html_report(
    results: pd.DataFrame,
    norm_counts: Optional[pd.DataFrame],
    metadata: pd.DataFrame,
    output_dir: str,
    design: str,
    title: str = "RNA-seq Analysis Report",
    tool_version: str = "0.1.0",
) -> str:
    """
    Generate a self-contained HTML report.

    Parameters
    ----------
    results : pd.DataFrame
        DEA results table.
    norm_counts : pd.DataFrame or None
        Normalized count matrix.
    metadata : pd.DataFrame
        Sample metadata.
    output_dir : str
        Directory containing plots/ subdirectory and where report will be saved.
    design : str
        Design formula used in analysis.
    title : str
        Report title.
    tool_version : str
        Package version string.

    Returns
    -------
    str
        Path to the generated HTML report.
    """
    plots_dir = os.path.join(output_dir, "plots")
    report_path = os.path.join(output_dir, "report.html")

    # Compute summary stats
    n_genes_total = len(results)
    sig_up = ((results["padj"] < 0.05) & (results["log2FoldChange"] > 1)).sum()
    sig_down = ((results["padj"] < 0.05) & (results["log2FoldChange"] < -1)).sum()
    n_samples = len(metadata)

    # Load plots as base64
    def _img_tag(name, alt):
        path = os.path.join(plots_dir, name)
        b64 = _img_to_base64(path)
        if not b64:
            return f'<p class="text-muted">Plot not available: {name}</p>'
        return (f'<figure class="figure">'
                f'<img src="{b64}" class="img-fluid figure-img" alt="{alt}">'
                f'<figcaption class="figure-caption text-center">{alt}</figcaption>'
                f'</figure>')

    # Top 20 DEGs table
    top_degs = results.dropna(subset=["padj"]).nsmallest(20, "padj")

    # GO enrichment table
    go_path = os.path.join(output_dir, "enrichment", "go", "go_enrichment.csv")
    go_html = ""
    if os.path.exists(go_path):
        go_df = pd.read_csv(go_path).head(20)
        go_html = _df_to_html_table(go_df)

    # KEGG enrichment table
    kegg_path = os.path.join(output_dir, "enrichment", "kegg", "kegg_enrichment.csv")
    kegg_html = ""
    if os.path.exists(kegg_path):
        kegg_df = pd.read_csv(kegg_path).head(20)
        kegg_html = _df_to_html_table(kegg_df)

    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>{title}</title>
  <link rel="stylesheet"
    href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css">
  <style>
    body {{ font-family: 'Segoe UI', sans-serif; background: #f8f9fa; }}
    .hero {{ background: linear-gradient(135deg, #1a3a5c 0%, #2d7dd2 100%);
             color: white; padding: 3rem 2rem; margin-bottom: 2rem; }}
    .stat-card {{ background: white; border-radius: 8px; padding: 1.5rem;
                  text-align: center; box-shadow: 0 2px 8px rgba(0,0,0,0.08); }}
    .stat-number {{ font-size: 2.5rem; font-weight: bold; }}
    .stat-up {{ color: #d62728; }}
    .stat-down {{ color: #1f77b4; }}
    .stat-total {{ color: #2ca02c; }}
    .section {{ background: white; border-radius: 8px; padding: 2rem;
                margin-bottom: 1.5rem; box-shadow: 0 2px 8px rgba(0,0,0,0.06); }}
    .table {{ font-size: 0.85rem; }}
    .badge-sig {{ background: #d62728; color: white; padding: 2px 6px;
                  border-radius: 4px; font-size: 0.75rem; }}
    footer {{ text-align: center; padding: 2rem; color: #666; font-size: 0.85rem; }}
  </style>
</head>
<body>

<div class="hero">
  <div class="container">
    <h1 class="display-5 fw-bold">&#x1F9EC; {title}</h1>
    <p class="lead">Generated by <strong>rnaseq-toolkit v{tool_version}</strong>
      &nbsp;|&nbsp; {timestamp}</p>
    <p class="mb-0">Design: <code>{design}</code></p>
  </div>
</div>

<div class="container">

  <!-- Summary stats -->
  <div class="section">
    <h2 class="h4 mb-4">&#x1F4CA; Analysis Summary</h2>
    <div class="row g-3">
      <div class="col-md-3">
        <div class="stat-card">
          <div class="stat-number stat-total">{n_samples}</div>
          <div class="text-muted">Samples</div>
        </div>
      </div>
      <div class="col-md-3">
        <div class="stat-card">
          <div class="stat-number stat-total">{n_genes_total:,}</div>
          <div class="text-muted">Genes Tested</div>
        </div>
      </div>
      <div class="col-md-3">
        <div class="stat-card">
          <div class="stat-number stat-up">{int(sig_up):,}</div>
          <div class="text-muted">Upregulated (padj&lt;0.05, |LFC|&gt;1)</div>
        </div>
      </div>
      <div class="col-md-3">
        <div class="stat-card">
          <div class="stat-number stat-down">{int(sig_down):,}</div>
          <div class="text-muted">Downregulated (padj&lt;0.05, |LFC|&gt;1)</div>
        </div>
      </div>
    </div>
  </div>

  <!-- Sample metadata -->
  <div class="section">
    <h2 class="h4 mb-3">&#x1F9EA; Sample Metadata</h2>
    {_df_to_html_table(metadata)}
  </div>

  <!-- Plots -->
  <div class="section">
    <h2 class="h4 mb-4">&#x1F4C8; Visualizations</h2>
    <div class="row g-4">
      <div class="col-md-6">{_img_tag("volcano.png", "Volcano Plot")}</div>
      <div class="col-md-6">{_img_tag("ma_plot.png", "MA Plot")}</div>
      <div class="col-md-6">{_img_tag("pca.png", "PCA of Normalized Counts")}</div>
      <div class="col-md-6">{_img_tag("heatmap.png", "Top DEGs Heatmap")}</div>
    </div>
  </div>

  <!-- GSEA plots -->
  <div class="section">
    <h2 class="h4 mb-4">&#x1F9EC; Pathway Enrichment Plots</h2>
    <div class="row g-4">
      <div class="col-md-6">
        {_img_tag("gsea_KEGG_2021_Human.png", "GSEA — KEGG 2021 Human")}
      </div>
      <div class="col-md-6">
        {_img_tag("gsea_GO_Biological_Process_2021.png",
                  "GSEA — GO Biological Process 2021")}
      </div>
    </div>
  </div>

  <!-- Top DEGs table -->
  <div class="section">
    <h2 class="h4 mb-3">&#x1F3AF; Top 20 Differentially Expressed Genes</h2>
    {_df_to_html_table(top_degs)}
  </div>

  <!-- GO enrichment -->
  {"<div class='section'><h2 class='h4 mb-3'>&#x1F9EC; GO Enrichment (top 20)</h2>" + go_html + "</div>" if go_html else ""}

  <!-- KEGG enrichment -->
  {"<div class='section'><h2 class='h4 mb-3'>&#x1F5FA;&#xFE0F; KEGG Pathway Enrichment (top 20)</h2>" + kegg_html + "</div>" if kegg_html else ""}

  <!-- Methods -->
  <div class="section">
    <h2 class="h4 mb-3">&#x1F4D6; Methods</h2>
    <p>Differential expression analysis was performed using
    <strong>rnaseq-toolkit v{tool_version}</strong>, which provides a unified
    interface to PyDESeq2 (Love <em>et al.</em>, 2014) and an edgeR-like
    quasi-likelihood pipeline. Low-count genes were removed prior to analysis
    (minimum 10 counts across all samples, present in at least 2 samples).
    Normalization was performed using the median-of-ratios method implemented
    in PyDESeq2. Multiple testing correction was applied using the
    Benjamini–Hochberg procedure. Pathway enrichment was assessed using
    pre-ranked GSEA (Subramanian <em>et al.</em>, 2005) via gseapy, and
    over-representation analysis was performed using Enrichr gene set libraries.
    All analyses were conducted in Python 3.12 on Ubuntu 24.04 LTS.</p>
  </div>

</div>

<footer>
  Generated by <strong>rnaseq-toolkit v{tool_version}</strong> &mdash;
  <a href="https://github.com/rnaseq-toolkit/rnaseq-toolkit">GitHub</a>
  &mdash; {timestamp}
</footer>

<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/js/bootstrap.bundle.min.js"></script>
</body>
</html>"""

    with open(report_path, "w", encoding="utf-8") as f:
        f.write(html)

    print(f"  HTML report saved to {report_path}")
    return report_path
