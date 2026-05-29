"""
cli.py — Command-line interface for rnaseq-toolkit.

Usage
-----
    rnaseq-toolkit \\
        --counts data/counts.csv \\
        --metadata data/metadata.csv \\
        --design "~condition" \\
        --contrast condition treated control \\
        --norm-method deseq2 \\
        --dea-method deseq2 \\
        --gene-sets KEGG_2021_Human GO_Biological_Process_2021 \\
        --output results/

    rnaseq-toolkit --help
"""

import argparse
import sys
import os


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="rnaseq-toolkit",
        description=(
            "rnaseq-toolkit: Streamlined differential expression analysis "
            "and pathway enrichment visualization from RNA-seq data."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic DESeq2 analysis
  rnaseq-toolkit --counts counts.csv --metadata meta.csv \\
                 --design "~condition" --contrast condition treated control

  # Switch to edgeR-like method with TMM normalization
  rnaseq-toolkit --counts counts.csv --metadata meta.csv \\
                 --design "~condition" --norm-method tmm --dea-method edger

  # Custom gene sets and output directory
  rnaseq-toolkit --counts counts.csv --metadata meta.csv \\
                 --design "~condition" --output my_results/ \\
                 --gene-sets KEGG_2021_Human Reactome_2022
        """,
    )

    # Required
    req = parser.add_argument_group("Required arguments")
    req.add_argument("--counts", required=True,
                     help="Path to raw count matrix (CSV or TSV). "
                          "Genes as rows, samples as columns.")
    req.add_argument("--metadata", required=True,
                     help="Path to sample metadata (CSV or TSV). "
                          "Samples as rows, variables as columns.")
    req.add_argument("--design", required=True,
                     help="Design formula, e.g. '~condition'.")

    # Optional
    opt = parser.add_argument_group("Optional arguments")
    opt.add_argument("--contrast", nargs=3,
                     metavar=("FACTOR", "NUMERATOR", "DENOMINATOR"),
                     help="Contrast for DEA, e.g. --contrast condition treated control")
    opt.add_argument("--norm-method", default="deseq2",
                     choices=["deseq2", "tmm", "cpm", "tpm", "rpkm", "vst", "rlog"],
                     help="Normalization method (default: deseq2).")
    opt.add_argument("--dea-method", default="deseq2",
                     choices=["deseq2", "edger"],
                     help="Differential expression method (default: deseq2).")
    opt.add_argument("--gene-sets", nargs="+",
                     default=["KEGG_2021_Human", "GO_Biological_Process_2021"],
                     help="Gene set libraries for enrichment analysis.")
    opt.add_argument("--lfc-threshold", type=float, default=1.0,
                     help="Log2 fold-change threshold for significance (default: 1.0).")
    opt.add_argument("--pval-threshold", type=float, default=0.05,
                     help="Adjusted p-value threshold (default: 0.05).")
    opt.add_argument("--min-count", type=int, default=10,
                     help="Minimum total count per gene for filtering (default: 10).")
    opt.add_argument("--min-samples", type=int, default=2,
                     help="Minimum samples with count > 0 per gene (default: 2).")
    opt.add_argument("--output", default="results",
                     help="Output directory (default: results/).")
    opt.add_argument("--no-report", action="store_true",
                     help="Skip HTML report generation.")
    opt.add_argument("--cpus", type=int, default=2,
                     help="Number of CPUs for parallel computation (default: 2).")
    opt.add_argument("--version", action="version", version="rnaseq-toolkit 0.1.0")

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    # Import here to allow --help without loading heavy dependencies
    try:
        from .pipeline import RNAseqPipeline
    except ImportError as e:
        print(f"[Error] Failed to import rnaseq_toolkit: {e}", file=sys.stderr)
        print("Install dependencies: pip install rnaseq_toolkit", file=sys.stderr)
        sys.exit(1)

    # Validate inputs
    for path_arg, name in [(args.counts, "--counts"), (args.metadata, "--metadata")]:
        if not os.path.exists(path_arg):
            print(f"[Error] File not found: {path_arg} ({name})", file=sys.stderr)
            sys.exit(1)

    print(f"\nrnaseq-toolkit v0.1.0")
    print(f"Counts:      {args.counts}")
    print(f"Metadata:    {args.metadata}")
    print(f"Design:      {args.design}")
    print(f"Contrast:    {args.contrast}")
    print(f"Norm method: {args.norm_method}")
    print(f"DEA method:  {args.dea_method}")
    print(f"Gene sets:   {', '.join(args.gene_sets)}")
    print(f"Output:      {args.output}\n")

    pipe = RNAseqPipeline(
        counts_path=args.counts,
        metadata_path=args.metadata,
        design=args.design,
        output_dir=args.output,
        min_count=args.min_count,
        min_samples=args.min_samples,
    )

    pipe.run(
        norm_method=args.norm_method,
        dea_method=args.dea_method,
        contrast=args.contrast,
        gene_sets=args.gene_sets,
        lfc_threshold=args.lfc_threshold,
        pval_threshold=args.pval_threshold,
        generate_report=not args.no_report,
    )


if __name__ == "__main__":
    main()
