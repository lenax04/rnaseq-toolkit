"""
visualization.py — Publication-quality plots for RNA-seq analysis.

All functions accept a results DataFrame and return a matplotlib Figure,
optionally saving to disk.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from typing import Optional, List, Tuple
import os


# ── Shared style ───────────────────────────────────────────────────────────────

def _set_style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 11,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.dpi": 150,
    })


# ── Volcano plot ───────────────────────────────────────────────────────────────

def plot_volcano(
    results: pd.DataFrame,
    lfc_col: str = "log2FoldChange",
    pval_col: str = "padj",
    lfc_threshold: float = 1.0,
    pval_threshold: float = 0.05,
    top_n_labels: int = 10,
    title: str = "Volcano Plot",
    output_path: Optional[str] = None,
    figsize: Tuple[int, int] = (9, 7),
) -> plt.Figure:
    """
    Generate a publication-quality volcano plot.

    Parameters
    ----------
    results : pd.DataFrame
        DEA results with log2FoldChange and padj columns. Index = gene names.
    lfc_col : str
        Column name for log2 fold-change values.
    pval_col : str
        Column name for adjusted p-values.
    lfc_threshold : float
        Absolute log2FC threshold for significance (default 1.0).
    pval_threshold : float
        Adjusted p-value threshold (default 0.05).
    top_n_labels : int
        Number of top significant genes to label.
    title : str
        Plot title.
    output_path : str, optional
        If provided, save figure to this path.
    figsize : tuple
        Figure size in inches.

    Returns
    -------
    matplotlib.figure.Figure
    """
    _set_style()
    res = results.dropna(subset=[lfc_col, pval_col]).copy()
    res["neg_log10_p"] = -np.log10(res[pval_col].clip(lower=1e-300))

    up = (res[pval_col] < pval_threshold) & (res[lfc_col] > lfc_threshold)
    down = (res[pval_col] < pval_threshold) & (res[lfc_col] < -lfc_threshold)
    ns = ~(up | down)

    fig, ax = plt.subplots(figsize=figsize)

    ax.scatter(res.loc[ns, lfc_col], res.loc[ns, "neg_log10_p"],
               c="#AAAAAA", alpha=0.4, s=12, linewidths=0, label="NS")
    ax.scatter(res.loc[up, lfc_col], res.loc[up, "neg_log10_p"],
               c="#D62728", alpha=0.7, s=18, linewidths=0,
               label=f"Up ({up.sum()})")
    ax.scatter(res.loc[down, lfc_col], res.loc[down, "neg_log10_p"],
               c="#1F77B4", alpha=0.7, s=18, linewidths=0,
               label=f"Down ({down.sum()})")

    # Threshold lines
    ax.axhline(-np.log10(pval_threshold), color="black", linestyle="--",
               linewidth=0.8, alpha=0.6)
    ax.axvline(lfc_threshold, color="black", linestyle="--",
               linewidth=0.8, alpha=0.6)
    ax.axvline(-lfc_threshold, color="black", linestyle="--",
               linewidth=0.8, alpha=0.6)

    # Label top genes
    sig_genes = res[up | down].copy()
    sig_genes = sig_genes.nsmallest(top_n_labels, pval_col)
    for gene, row in sig_genes.iterrows():
        ax.annotate(
            str(gene),
            xy=(row[lfc_col], row["neg_log10_p"]),
            xytext=(5, 3), textcoords="offset points",
            fontsize=7, color="black",
            arrowprops=dict(arrowstyle="-", color="gray", lw=0.5),
        )

    ax.set_xlabel(r"$\log_2$ Fold Change", fontsize=12)
    ax.set_ylabel(r"$-\log_{10}$(adjusted p-value)", fontsize=12)
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.legend(frameon=False, fontsize=10)

    plt.tight_layout()
    if output_path:
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    return fig


# ── MA plot ────────────────────────────────────────────────────────────────────

def plot_ma(
    results: pd.DataFrame,
    lfc_col: str = "log2FoldChange",
    mean_col: str = "baseMean",
    pval_col: str = "padj",
    pval_threshold: float = 0.05,
    title: str = "MA Plot",
    output_path: Optional[str] = None,
    figsize: Tuple[int, int] = (9, 6),
) -> plt.Figure:
    """MA plot (mean expression vs. log2 fold-change)."""
    _set_style()
    res = results.dropna(subset=[lfc_col, mean_col, pval_col]).copy()
    sig = res[pval_col] < pval_threshold

    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(np.log2(res.loc[~sig, mean_col] + 1), res.loc[~sig, lfc_col],
               c="#AAAAAA", alpha=0.4, s=10, linewidths=0, label="NS")
    ax.scatter(np.log2(res.loc[sig, mean_col] + 1), res.loc[sig, lfc_col],
               c="#D62728", alpha=0.7, s=14, linewidths=0,
               label=f"Significant ({sig.sum()})")
    ax.axhline(0, color="black", linewidth=1)
    ax.set_xlabel(r"$\log_2$ Mean Expression", fontsize=12)
    ax.set_ylabel(r"$\log_2$ Fold Change", fontsize=12)
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.legend(frameon=False)
    plt.tight_layout()
    if output_path:
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    return fig


# ── PCA plot ───────────────────────────────────────────────────────────────────

def plot_pca(
    norm_counts: pd.DataFrame,
    metadata: pd.DataFrame,
    color_by: str,
    n_top_genes: int = 500,
    title: str = "PCA of Normalized Counts",
    output_path: Optional[str] = None,
    figsize: Tuple[int, int] = (8, 6),
) -> plt.Figure:
    """
    PCA plot of normalized count data, coloured by a metadata variable.

    Parameters
    ----------
    norm_counts : pd.DataFrame
        Normalized count matrix (genes × samples).
    metadata : pd.DataFrame
        Sample metadata. Index must match columns of norm_counts.
    color_by : str
        Column in metadata used for colouring points.
    n_top_genes : int
        Number of most-variable genes to use for PCA (default 500).
    """
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    _set_style()

    # Select most variable genes
    common = norm_counts.columns.intersection(metadata.index)
    mat = norm_counts[common].T  # samples × genes
    variances = mat.var(axis=0)
    top_genes = variances.nlargest(min(n_top_genes, len(variances))).index
    mat_top = mat[top_genes]

    # Scale and PCA
    scaled = StandardScaler().fit_transform(mat_top)
    pca = PCA(n_components=2)
    coords = pca.fit_transform(scaled)
    var_exp = pca.explained_variance_ratio_ * 100

    groups = metadata.loc[common, color_by].astype(str)
    unique_groups = groups.unique()
    palette = plt.cm.get_cmap("tab10", len(unique_groups))
    color_map = {g: palette(i) for i, g in enumerate(unique_groups)}
    colors = [color_map[g] for g in groups]

    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(coords[:, 0], coords[:, 1], c=colors, s=60, alpha=0.85,
               edgecolors="white", linewidths=0.5)

    # Legend
    patches = [mpatches.Patch(color=color_map[g], label=g)
               for g in unique_groups]
    ax.legend(handles=patches, frameon=False, fontsize=9,
              title=color_by, title_fontsize=10)

    ax.set_xlabel(f"PC1 ({var_exp[0]:.1f}%)", fontsize=12)
    ax.set_ylabel(f"PC2 ({var_exp[1]:.1f}%)", fontsize=12)
    ax.set_title(title, fontsize=13, fontweight="bold")
    plt.tight_layout()
    if output_path:
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    return fig


# ── Heatmap ────────────────────────────────────────────────────────────────────

def plot_heatmap(
    norm_counts: pd.DataFrame,
    metadata: pd.DataFrame,
    results: pd.DataFrame,
    color_by: str,
    n_top_genes: int = 50,
    pval_threshold: float = 0.05,
    lfc_threshold: float = 1.0,
    title: str = "Top DEGs Heatmap",
    output_path: Optional[str] = None,
    figsize: Tuple[int, int] = (12, 10),
) -> plt.Figure:
    """
    Clustered heatmap of top differentially expressed genes.
    """
    _set_style()

    # Select top DEGs
    sig = results.dropna(subset=["padj", "log2FoldChange"])
    sig = sig[(sig["padj"] < pval_threshold) &
              (sig["log2FoldChange"].abs() > lfc_threshold)]
    sig = sig.nsmallest(n_top_genes, "padj")

    common = norm_counts.columns.intersection(metadata.index)
    mat = norm_counts.loc[sig.index.intersection(norm_counts.index), common]

    if mat.empty:
        print("[Heatmap] No significant genes found with current thresholds.")
        return plt.figure()

    # Z-score per gene
    mat_z = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1) + 1e-9, axis=0)

    # Annotation bar
    groups = metadata.loc[common, color_by].astype(str)
    unique_groups = groups.unique()
    palette = sns.color_palette("tab10", len(unique_groups))
    col_colors = pd.Series(
        [dict(zip(unique_groups, palette))[g] for g in groups],
        index=common
    )

    g = sns.clustermap(
        mat_z,
        col_colors=col_colors,
        cmap="RdBu_r",
        center=0,
        vmin=-3, vmax=3,
        yticklabels=True,
        xticklabels=True,
        figsize=figsize,
        dendrogram_ratio=(0.1, 0.05),
        cbar_pos=(0.02, 0.8, 0.03, 0.15),
    )
    g.ax_heatmap.set_title(title, fontsize=13, fontweight="bold", pad=20)
    g.ax_heatmap.set_xlabel("Sample", fontsize=10)
    g.ax_heatmap.set_ylabel("Gene", fontsize=10)
    g.ax_heatmap.tick_params(axis="y", labelsize=7)
    g.ax_heatmap.tick_params(axis="x", labelsize=8, rotation=45)

    if output_path:
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        g.savefig(output_path, dpi=300, bbox_inches="tight")
    return g.fig


# ── GSEA dotplot ───────────────────────────────────────────────────────────────

def plot_gsea_dotplot(
    gsea_results,
    n_top: int = 20,
    title: str = "GSEA Enrichment",
    output_path: Optional[str] = None,
    figsize: Tuple[int, int] = (10, 8),
) -> plt.Figure:
    """
    Dotplot of top GSEA enriched pathways.

    Parameters
    ----------
    gsea_results : gseapy PreRank result or pd.DataFrame
        Result from gseapy.prerank() or a DataFrame with columns:
        Term, NES, FDR q-val, Tag %.
    """
    _set_style()

    if hasattr(gsea_results, "res2d"):
        df = gsea_results.res2d.copy()
    else:
        df = gsea_results.copy()

    # Standardise column names
    col_map = {}
    for c in df.columns:
        cl = c.lower().replace(" ", "_").replace("-", "_")
        if "nes" in cl:
            col_map[c] = "NES"
        elif "fdr" in cl or "q_val" in cl or "qval" in cl:
            col_map[c] = "FDR"
        elif "tag" in cl or "gene_percent" in cl:
            col_map[c] = "GeneRatio"
        elif "term" in cl or "pathway" in cl:
            col_map[c] = "Term"
    df = df.rename(columns=col_map)

    needed = {"NES", "FDR", "Term"}
    missing = needed - set(df.columns)
    if missing:
        print(f"[GSEA dotplot] Missing columns: {missing}. Skipping plot.")
        return plt.figure()

    df = df[df["FDR"] < 0.25].copy()
    df = df.reindex(df["NES"].abs().nlargest(n_top).index)
    df = df.sort_values("NES")

    colors = ["#1F77B4" if n < 0 else "#D62728" for n in df["NES"]]
    sizes = df.get("GeneRatio", pd.Series([10] * len(df))).fillna(10)
    if isinstance(sizes, pd.Series) and sizes.dtype == object:
        sizes = sizes.str.replace("%", "").astype(float)

    fig, ax = plt.subplots(figsize=figsize)
    sc = ax.scatter(df["NES"], range(len(df)), c=colors,
                    s=sizes * 3, alpha=0.8, edgecolors="white", linewidths=0.5)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["Term"].str[:60], fontsize=8)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Normalized Enrichment Score (NES)", fontsize=11)
    ax.set_title(title, fontsize=13, fontweight="bold")

    # Legend for direction
    up_patch = mpatches.Patch(color="#D62728", label="Upregulated")
    dn_patch = mpatches.Patch(color="#1F77B4", label="Downregulated")
    ax.legend(handles=[up_patch, dn_patch], frameon=False, fontsize=9)

    plt.tight_layout()
    if output_path:
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    return fig
