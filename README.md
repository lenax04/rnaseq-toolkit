# rnaseq-toolkit

[![CI](https://github.com/rnaseq-toolkit/rnaseq-toolkit/actions/workflows/ci.yml/badge.svg)](https://github.com/rnaseq-toolkit/rnaseq-toolkit/actions)
[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-pending-lightgrey)](CITATION.cff)

**rnaseq-toolkit** is an open-source Python package that provides a unified, switchable interface for RNA-seq differential expression analysis (DEA) and pathway enrichment visualization. It automates the complete pipeline from raw count matrices through normalization, DEA, and pathway analysis to self-contained HTML reports — all with a single command or a few lines of Python.

## Key Features

- **Unified normalization interface**: switch between DESeq2 median-of-ratios, TMM, CPM, TPM, RPKM, VST, and rlog with one argument (`--norm-method`)
- **Switchable DEA methods**: DESeq2 (via PyDESeq2) and edgeR-like quasi-likelihood (via statsmodels NB GLM) with identical output format
- **Pathway enrichment**: pre-ranked GSEA, GO over-representation, and KEGG pathway analysis via gseapy/Enrichr
- **Publication-quality plots**: volcano, MA, PCA, clustered heatmap, GSEA dotplot — all at 300 DPI
- **Automated HTML reports**: self-contained reports with embedded plots and interactive tables
- **Snakemake workflow**: reproducible, scalable pipeline management
- **Docker support**: fully containerized for reproducibility
- **Benchmarked**: validated on three public GEO datasets (Batten disease, COVID-19, breast cancer)

## Installation

```bash
pip install rnaseq_toolkit
```

Or install from source:

```bash
git clone https://github.com/rnaseq-toolkit/rnaseq-toolkit.git
cd rnaseq-toolkit
pip install -e ".[dev]"
```

## Quick Start

### Command-line interface

```bash
# Basic DESeq2 analysis
rnaseq-toolkit \
    --counts data/counts.csv \
    --metadata data/metadata.csv \
    --design "~condition" \
    --contrast condition treated control \
    --output results/

# Switch to edgeR-like with TMM normalization — one argument change
rnaseq-toolkit \
    --counts data/counts.csv \
    --metadata data/metadata.csv \
    --design "~condition" \
    --norm-method tmm \
    --dea-method edger \
    --output results_edger/
```

### Python API

```python
from rnaseq_toolkit import RNAseqPipeline

# Initialize pipeline
pipe = RNAseqPipeline(
    counts_path="data/counts.csv",
    metadata_path="data/metadata.csv",
    design="~condition",
    output_dir="results/",
)

# Run complete pipeline with DESeq2
pipe.run(
    norm_method="deseq2",
    dea_method="deseq2",
    contrast=["condition", "treated", "control"],
    gene_sets=["KEGG_2021_Human", "GO_Biological_Process_2021"],
)

# Switch to edgeR-like with TMM — one line change
pipe.run(norm_method="tmm", dea_method="edger")
```

### Individual modules

```python
from rnaseq_toolkit import normalize_counts, run_deseq2, plot_volcano
import pandas as pd

counts = pd.read_csv("data/counts.csv", index_col=0)
meta   = pd.read_csv("data/metadata.csv", index_col=0)

# Normalize
norm = normalize_counts(counts, method="deseq2", metadata=meta, design="~condition")

# DEA
results = run_deseq2(counts, meta, "~condition",
                     contrast=["condition", "treated", "control"])

# Visualize
plot_volcano(results, output_path="results/volcano.png")
```

## Input Format

**Count matrix** (`counts.csv`): genes as rows, samples as columns.

```
gene_id,Sample_01,Sample_02,Sample_03,Sample_04,Sample_05,Sample_06
ENSG00000000003,1234,1456,1123,2345,2567,2234
ENSG00000000005,456,512,489,234,198,267
...
```

**Metadata** (`metadata.csv`): samples as rows, variables as columns.

```
sample_id,condition,batch
Sample_01,control,A
Sample_02,control,A
Sample_03,control,B
Sample_04,treated,A
Sample_05,treated,A
Sample_06,treated,B
```

## Normalization Methods

| Method   | Description                                   | Use case                         |
|----------|-----------------------------------------------|----------------------------------|
| `deseq2` | Median-of-ratios (DESeq2)                     | Default; DEA with PyDESeq2       |
| `tmm`    | Trimmed Mean of M-values (edgeR)              | DEA with edgeR-like              |
| `cpm`    | Counts Per Million                            | Quick exploration                |
| `tpm`    | Transcripts Per Million (needs gene lengths)  | Cross-sample comparison          |
| `rpkm`   | RPKM (needs gene lengths)                     | Legacy compatibility             |
| `vst`    | Variance Stabilizing Transformation           | PCA, heatmaps, clustering        |
| `rlog`   | Regularized log (approximate)                 | Small sample sizes               |

## Snakemake Workflow

```bash
# Edit workflow/config.yaml, then run:
snakemake --cores 4 --configfile workflow/config.yaml
```

## Docker

```bash
docker build -t rnaseq-toolkit:0.1.0 .

docker run --rm \
    -v $(pwd)/data:/workspace/data \
    -v $(pwd)/results:/workspace/results \
    rnaseq-toolkit:0.1.0 \
    --counts data/counts.csv \
    --metadata data/metadata.csv \
    --design "~condition" \
    --contrast condition treated control \
    --output results/
```

## Benchmark

rnaseq-toolkit was benchmarked against standalone DESeq2 and edgeR-like implementations on three public GEO datasets:

| Dataset   | Disease / Condition | GEO Accession | Samples |
|-----------|--------------------|--------------:|--------:|
| Batten    | CLN2 (rare disease) | GSE210143    | 12      |
| COVID-19  | SARS-CoV-2 infection | GSE157103   | 126     |
| Breast Ca | TNBC vs. luminal   | GSE183947    | 60      |

Results demonstrate high concordance (Pearson r > 0.95 for log2FC) between rnaseq-toolkit and standalone tools, with significantly reduced setup time.

## Citation

If you use rnaseq-toolkit in your research, please cite:

```bibtex
@software{rnaseq_toolkit_2026,
  title   = {rnaseq-toolkit: Streamlined Differential Expression Analysis
             and Pathway Enrichment Visualization from RNA-seq Data},
  version = {0.1.0},
  year    = {2026},
  url     = {https://github.com/rnaseq-toolkit/rnaseq-toolkit},
  license = {MIT}
}
```

See also [CITATION.cff](CITATION.cff) for full citation metadata.

## Dependencies

- [PyDESeq2](https://github.com/owkin/PyDESeq2) — DESeq2 in Python
- [gseapy](https://github.com/zqfang/GSEApy) — GSEA and Enrichr
- [statsmodels](https://www.statsmodels.org) — NB GLM for edgeR-like
- [scikit-learn](https://scikit-learn.org) — PCA
- [seaborn](https://seaborn.pydata.org) / [matplotlib](https://matplotlib.org) — visualization

## License

MIT License — see [LICENSE](LICENSE).
