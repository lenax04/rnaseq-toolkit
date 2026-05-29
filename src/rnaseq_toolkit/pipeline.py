"""
pipeline.py — High-level RNAseqPipeline class.

Provides a single-entry-point API that orchestrates:
  1. Data loading and quality control
  2. Normalization (switchable via one argument)
  3. Differential expression analysis (DESeq2 or edgeR-like)
  4. Pathway enrichment (GSEA, GO, KEGG)
  5. Visualization (volcano, MA, PCA, heatmap, GSEA dotplot)
  6. HTML report generation

Example
-------
>>> from rnaseq_toolkit import RNAseqPipeline
>>> pipe = RNAseqPipeline(
...     counts_path="data/counts.csv",
...     metadata_path="data/metadata.csv",
...     design="~condition",
...     output_dir="results/",
... )
>>> pipe.run(
...     norm_method="deseq2",
...     dea_method="deseq2",
...     contrast=["condition", "treated", "control"],
...     gene_sets=["KEGG_2021_Human", "GO_Biological_Process_2021"],
... )
"""

import os
import pandas as pd
import numpy as np
from typing import Optional, List, Literal

from .normalization import normalize_counts
from .dea import run_deseq2, run_edger_like
from .visualization import (
    plot_volcano, plot_ma, plot_pca, plot_heatmap, plot_gsea_dotplot
)
from .enrichment import run_gsea, run_go_enrichment, run_kegg_enrichment
from .report import generate_html_report


class RNAseqPipeline:
    """
    Streamlined RNA-seq analysis pipeline.

    Attributes
    ----------
    counts : pd.DataFrame
        Raw count matrix (genes × samples).
    metadata : pd.DataFrame
        Sample annotation table.
    design : str
        Design formula (e.g. '~condition').
    output_dir : str
        Root directory for all output files.
    results : pd.DataFrame or None
        DEA results after calling run() or run_dea().
    norm_counts : pd.DataFrame or None
        Normalized counts after calling run() or run_normalization().
    """

    def __init__(
        self,
        counts_path: str,
        metadata_path: str,
        design: str,
        output_dir: str = "results",
        min_count: int = 10,
        min_samples: int = 2,
    ):
        """
        Initialise the pipeline by loading and pre-filtering data.

        Parameters
        ----------
        counts_path : str
            Path to raw count matrix (CSV or TSV). Genes as rows, samples as columns.
        metadata_path : str
            Path to sample metadata (CSV or TSV). Samples as rows.
        design : str
            Design formula, e.g. '~condition'.
        output_dir : str
            Directory for output files. Created if it does not exist.
        min_count : int
            Minimum total count per gene across all samples for filtering (default 10).
        min_samples : int
            Minimum number of samples with count > 0 per gene (default 2).
        """
        self.design = design
        self.output_dir = output_dir
        self.results = None
        self.norm_counts = None
        self._gsea_result = None

        os.makedirs(output_dir, exist_ok=True)

        self.counts = self._load_table(counts_path)
        self.metadata = self._load_table(metadata_path)

        # Align samples
        common = self.counts.columns.intersection(self.metadata.index)
        if len(common) == 0:
            raise ValueError(
                "No common samples found between counts columns and metadata index. "
                "Check that sample names match exactly."
            )
        self.counts = self.counts[common]
        self.metadata = self.metadata.loc[common]

        # Pre-filter low-count genes
        keep = (
            (self.counts.sum(axis=1) >= min_count) &
            ((self.counts > 0).sum(axis=1) >= min_samples)
        )
        n_before = len(self.counts)
        self.counts = self.counts[keep]
        print(f"[Pipeline] Loaded {n_before} genes, "
              f"kept {len(self.counts)} after low-count filtering.")
        print(f"[Pipeline] Samples: {len(common)}")

    # ── Public API ─────────────────────────────────────────────────────────────

    def run(
        self,
        norm_method: str = "deseq2",
        dea_method: Literal["deseq2", "edger"] = "deseq2",
        contrast: Optional[List[str]] = None,
        gene_sets: Optional[List[str]] = None,
        lfc_threshold: float = 1.0,
        pval_threshold: float = 0.05,
        top_n_labels: int = 15,
        generate_report: bool = True,
    ) -> "RNAseqPipeline":
        """
        Run the complete pipeline in one call.

        Parameters
        ----------
        norm_method : str
            Normalization method. One of: 'deseq2', 'tmm', 'cpm', 'vst', 'rlog'.
        dea_method : str
            DEA method: 'deseq2' or 'edger'.
        contrast : list of str, optional
            Contrast for DEA, e.g. ['condition', 'treated', 'control'].
        gene_sets : list of str, optional
            Gene set libraries for enrichment analysis.
            Default: ['KEGG_2021_Human', 'GO_Biological_Process_2021'].
        lfc_threshold : float
            |log2FC| threshold for significance (default 1.0).
        pval_threshold : float
            Adjusted p-value threshold (default 0.05).
        top_n_labels : int
            Number of top genes to label in volcano plot.
        generate_report : bool
            Whether to generate an HTML report (default True).

        Returns
        -------
        self : RNAseqPipeline
            Returns self for method chaining.
        """
        if gene_sets is None:
            gene_sets = ["KEGG_2021_Human", "GO_Biological_Process_2021"]

        print("\n" + "=" * 60)
        print("  rnaseq-toolkit pipeline")
        print("=" * 60)

        # Step 1: Normalization
        print("\n[Step 1/5] Normalization...")
        self.run_normalization(method=norm_method)

        # Step 2: DEA
        print("\n[Step 2/5] Differential Expression Analysis...")
        self.run_dea(method=dea_method, contrast=contrast,
                     lfc_threshold=lfc_threshold, pval_threshold=pval_threshold)

        # Step 3: Visualization
        print("\n[Step 3/5] Generating plots...")
        self.run_visualization(
            lfc_threshold=lfc_threshold,
            pval_threshold=pval_threshold,
            top_n_labels=top_n_labels,
        )

        # Step 4: Enrichment
        print("\n[Step 4/5] Pathway enrichment analysis...")
        self.run_enrichment(
            gene_sets=gene_sets,
            lfc_threshold=lfc_threshold,
            pval_threshold=pval_threshold,
        )

        # Step 5: Report
        if generate_report:
            print("\n[Step 5/5] Generating HTML report...")
            self.run_report()

        print("\n[Pipeline] Complete. Results in:", self.output_dir)
        return self

    def run_normalization(self, method: str = "deseq2") -> pd.DataFrame:
        """Normalize counts and store in self.norm_counts."""
        self.norm_counts = normalize_counts(
            counts=self.counts,
            method=method,
            metadata=self.metadata,
            design=self.design,
        )
        out_path = os.path.join(self.output_dir, "normalized_counts.csv")
        self.norm_counts.to_csv(out_path)
        print(f"  Normalized counts saved to {out_path}")
        return self.norm_counts

    def run_dea(
        self,
        method: Literal["deseq2", "edger"] = "deseq2",
        contrast: Optional[List[str]] = None,
        lfc_threshold: float = 0.0,
        pval_threshold: float = 0.05,
    ) -> pd.DataFrame:
        """Run DEA and store results in self.results."""
        if method == "deseq2":
            self.results = run_deseq2(
                counts=self.counts,
                metadata=self.metadata,
                design=self.design,
                contrast=contrast,
                lfc_threshold=lfc_threshold,
                alpha=pval_threshold,
            )
        elif method == "edger":
            design_col = self.design.replace("~", "").strip()
            contrast_pair = None
            if contrast and len(contrast) >= 3:
                contrast_pair = (contrast[1], contrast[2])
            self.results = run_edger_like(
                counts=self.counts,
                metadata=self.metadata,
                design=design_col,
                contrast=contrast_pair,
                alpha=pval_threshold,
            )
        else:
            raise ValueError(f"Unknown DEA method: '{method}'. Use 'deseq2' or 'edger'.")

        out_path = os.path.join(self.output_dir, f"{method}_results.csv")
        self.results.to_csv(out_path)
        print(f"  DEA results saved to {out_path}")
        return self.results

    def run_visualization(
        self,
        lfc_threshold: float = 1.0,
        pval_threshold: float = 0.05,
        top_n_labels: int = 15,
    ):
        """Generate all standard plots."""
        if self.results is None:
            raise RuntimeError("Run DEA first (run_dea or run).")

        plots_dir = os.path.join(self.output_dir, "plots")
        os.makedirs(plots_dir, exist_ok=True)

        plot_volcano(
            self.results,
            lfc_threshold=lfc_threshold,
            pval_threshold=pval_threshold,
            top_n_labels=top_n_labels,
            output_path=os.path.join(plots_dir, "volcano.png"),
        )
        plot_ma(
            self.results,
            pval_threshold=pval_threshold,
            output_path=os.path.join(plots_dir, "ma_plot.png"),
        )

        if self.norm_counts is not None:
            design_col = self.design.replace("~", "").strip()
            plot_pca(
                self.norm_counts,
                self.metadata,
                color_by=design_col,
                output_path=os.path.join(plots_dir, "pca.png"),
            )
            plot_heatmap(
                self.norm_counts,
                self.metadata,
                self.results,
                color_by=design_col,
                lfc_threshold=lfc_threshold,
                pval_threshold=pval_threshold,
                output_path=os.path.join(plots_dir, "heatmap.png"),
            )
        print(f"  Plots saved to {plots_dir}/")

    def run_enrichment(
        self,
        gene_sets: Optional[List[str]] = None,
        lfc_threshold: float = 1.0,
        pval_threshold: float = 0.05,
    ):
        """Run GSEA, GO, and KEGG enrichment analyses."""
        if self.results is None:
            raise RuntimeError("Run DEA first.")

        if gene_sets is None:
            gene_sets = ["KEGG_2021_Human", "GO_Biological_Process_2021"]

        enrich_dir = os.path.join(self.output_dir, "enrichment")

        # GSEA (pre-ranked)
        for gs in gene_sets:
            try:
                gsea_res = run_gsea(
                    self.results,
                    gene_sets=gs,
                    output_dir=os.path.join(enrich_dir, "gsea", gs.replace(" ", "_")),
                )
                self._gsea_result = gsea_res
                # Dotplot
                plots_dir = os.path.join(self.output_dir, "plots")
                plot_gsea_dotplot(
                    gsea_res,
                    title=f"GSEA — {gs}",
                    output_path=os.path.join(
                        plots_dir, f"gsea_{gs.replace(' ', '_')}.png"
                    ),
                )
            except Exception as e:
                print(f"  [Warning] GSEA with {gs} failed: {e}")

        # GO and KEGG over-representation on significant genes
        sig_genes = self._get_sig_genes(lfc_threshold, pval_threshold)
        if sig_genes:
            run_go_enrichment(
                sig_genes,
                output_dir=os.path.join(enrich_dir, "go"),
            )
            run_kegg_enrichment(
                sig_genes,
                output_dir=os.path.join(enrich_dir, "kegg"),
            )

    def run_report(self):
        """Generate HTML report."""
        generate_html_report(
            results=self.results,
            norm_counts=self.norm_counts,
            metadata=self.metadata,
            output_dir=self.output_dir,
            design=self.design,
        )

    # ── Private helpers ────────────────────────────────────────────────────────

    def _load_table(self, path: str) -> pd.DataFrame:
        sep = "," if path.endswith(".csv") else "\t"
        df = pd.read_csv(path, sep=sep, index_col=0)
        return df

    def _get_sig_genes(
        self, lfc_threshold: float, pval_threshold: float
    ) -> List[str]:
        if self.results is None:
            return []
        sig = self.results.dropna(subset=["padj", "log2FoldChange"])
        sig = sig[
            (sig["padj"] < pval_threshold) &
            (sig["log2FoldChange"].abs() > lfc_threshold)
        ]
        return sig.index.tolist()
