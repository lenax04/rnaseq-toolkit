"""
enrichment.py — Pathway and gene-set enrichment analysis.

Provides unified wrappers for:
  - GSEA (pre-ranked, via gseapy)
  - GO over-representation analysis (via gseapy Enrichr)
  - KEGG pathway enrichment (via gseapy Enrichr)
"""

import numpy as np
import pandas as pd
import os
from typing import Optional, List, Literal


OrgDb = Literal["Human", "Mouse", "Rat"]


def run_gsea(
    results: pd.DataFrame,
    gene_sets: str = "KEGG_2021_Human",
    lfc_col: str = "log2FoldChange",
    pval_col: str = "pvalue",
    output_dir: str = "results/gsea",
    permutation_num: int = 1000,
    min_size: int = 5,
    max_size: int = 1000,
    threads: int = 2,
    seed: int = 42,
    quiet: bool = False,
):
    """
    Run pre-ranked GSEA using gseapy.

    The ranking metric is: sign(log2FC) × -log10(pvalue).

    Parameters
    ----------
    results : pd.DataFrame
        DEA results. Index = gene names.
    gene_sets : str or list
        Gene set library name(s) from Enrichr or a GMT file path.
        Common options: 'KEGG_2021_Human', 'GO_Biological_Process_2021',
        'Reactome_2022', 'MSigDB_Hallmark_2020'.
    lfc_col : str
        Column name for log2 fold-change.
    pval_col : str
        Column name for raw p-values.
    output_dir : str
        Directory to save GSEA output files.
    permutation_num : int
        Number of permutations for p-value estimation (default 1000).
    min_size : int
        Minimum gene set size (default 5).
    max_size : int
        Maximum gene set size (default 1000).
    threads : int
        Number of parallel threads.
    seed : int
        Random seed for reproducibility.
    quiet : bool
        Suppress output.

    Returns
    -------
    gseapy PreRank result object.
    """
    import gseapy as gp

    os.makedirs(output_dir, exist_ok=True)

    res = results.dropna(subset=[lfc_col, pval_col]).copy()
    res["rank_metric"] = (
        np.sign(res[lfc_col]) * (-np.log10(res[pval_col].clip(lower=1e-300)))
    )
    rnk = res[["rank_metric"]].sort_values("rank_metric", ascending=False)

    if not quiet:
        print(f"[GSEA] Running pre-ranked GSEA with gene set: {gene_sets}")
        print(f"[GSEA] Ranked genes: {len(rnk)}")

    pre_res = gp.prerank(
        rnk=rnk,
        gene_sets=gene_sets,
        threads=threads,
        min_size=min_size,
        max_size=max_size,
        permutation_num=permutation_num,
        outdir=output_dir,
        seed=seed,
        verbose=not quiet,
    )

    n_sig = (pre_res.res2d["FDR q-val"] < 0.25).sum() if hasattr(pre_res, "res2d") else "?"
    if not quiet:
        print(f"[GSEA] Significant pathways (FDR < 0.25): {n_sig}")

    return pre_res


def run_go_enrichment(
    gene_list: List[str],
    organism: OrgDb = "Human",
    go_categories: Optional[List[str]] = None,
    pval_threshold: float = 0.05,
    output_dir: str = "results/go",
    quiet: bool = False,
) -> pd.DataFrame:
    """
    Gene Ontology over-representation analysis via Enrichr (gseapy).

    Parameters
    ----------
    gene_list : list of str
        List of significant gene symbols.
    organism : str
        Organism: 'Human', 'Mouse', or 'Rat'.
    go_categories : list of str, optional
        GO categories to test. Default: all three (BP, MF, CC).
    pval_threshold : float
        Adjusted p-value threshold for filtering results.
    output_dir : str
        Directory to save results.
    quiet : bool
        Suppress output.

    Returns
    -------
    pd.DataFrame
        Combined GO enrichment results.
    """
    import gseapy as gp

    os.makedirs(output_dir, exist_ok=True)

    if go_categories is None:
        go_categories = [
            f"GO_Biological_Process_2021",
            f"GO_Molecular_Function_2021",
            f"GO_Cellular_Component_2021",
        ]

    all_results = []
    for cat in go_categories:
        try:
            enr = gp.enrichr(
                gene_list=gene_list,
                gene_sets=cat,
                organism=organism,
                outdir=os.path.join(output_dir, cat.replace(" ", "_")),
                no_plot=True,
                verbose=not quiet,
            )
            df = enr.results.copy()
            df["Category"] = cat
            all_results.append(df)
        except Exception as e:
            if not quiet:
                print(f"[GO] Warning: {cat} failed — {e}")

    if not all_results:
        return pd.DataFrame()

    combined = pd.concat(all_results, ignore_index=True)
    combined = combined[combined["Adjusted P-value"] < pval_threshold]
    combined = combined.sort_values("Adjusted P-value")

    out_path = os.path.join(output_dir, "go_enrichment.csv")
    combined.to_csv(out_path, index=False)
    if not quiet:
        print(f"[GO] Results saved to {out_path} ({len(combined)} terms)")

    return combined


def run_kegg_enrichment(
    gene_list: List[str],
    organism: OrgDb = "Human",
    pval_threshold: float = 0.05,
    output_dir: str = "results/kegg",
    quiet: bool = False,
) -> pd.DataFrame:
    """
    KEGG pathway over-representation analysis via Enrichr (gseapy).

    Parameters
    ----------
    gene_list : list of str
        List of significant gene symbols.
    organism : str
        Organism: 'Human', 'Mouse', or 'Rat'.
    pval_threshold : float
        Adjusted p-value threshold.
    output_dir : str
        Directory to save results.
    quiet : bool
        Suppress output.

    Returns
    -------
    pd.DataFrame
        KEGG enrichment results.
    """
    import gseapy as gp

    os.makedirs(output_dir, exist_ok=True)

    kegg_lib = "KEGG_2021_Human" if organism == "Human" else f"KEGG_2019_{organism}"

    try:
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=kegg_lib,
            organism=organism,
            outdir=output_dir,
            no_plot=True,
            verbose=not quiet,
        )
        df = enr.results.copy()
    except Exception as e:
        if not quiet:
            print(f"[KEGG] Error: {e}")
        return pd.DataFrame()

    df = df[df["Adjusted P-value"] < pval_threshold]
    df = df.sort_values("Adjusted P-value")

    out_path = os.path.join(output_dir, "kegg_enrichment.csv")
    df.to_csv(out_path, index=False)
    if not quiet:
        print(f"[KEGG] Results saved to {out_path} ({len(df)} pathways)")

    return df
