"""
rnaseq_toolkit — Streamlined differential expression analysis and
pathway enrichment visualization from RNA-seq data.

Version: 0.1.0
License: MIT
"""

__version__ = "0.1.0"
__author__ = "RNAseqKit Contributors"
__email__ = "contact@rnaseqkit.org"

from .pipeline import RNAseqPipeline
from .normalization import normalize_counts
from .dea import run_deseq2, run_edger_like
from .visualization import (
    plot_volcano,
    plot_heatmap,
    plot_pca,
    plot_ma,
    plot_gsea_dotplot,
)
from .enrichment import run_gsea, run_go_enrichment, run_kegg_enrichment
from .report import generate_html_report

__all__ = [
    "RNAseqPipeline",
    "normalize_counts",
    "run_deseq2",
    "run_edger_like",
    "plot_volcano",
    "plot_heatmap",
    "plot_pca",
    "plot_ma",
    "plot_gsea_dotplot",
    "run_gsea",
    "run_go_enrichment",
    "run_kegg_enrichment",
    "generate_html_report",
]
