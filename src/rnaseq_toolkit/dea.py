"""
dea.py — Differential Expression Analysis interface.

Provides a unified API for:
  - PyDESeq2 (primary, Python-native)
  - edgeR-like quasi-likelihood (approximated via statsmodels NB GLM)

Both methods return a standardised DataFrame with columns:
  gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
"""

import numpy as np
import pandas as pd
from typing import Optional, List, Tuple
import warnings


def run_deseq2(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    design: str,
    contrast: Optional[List[str]] = None,
    lfc_threshold: float = 0.0,
    alpha: float = 0.05,
    n_cpus: int = 2,
    quiet: bool = False,
) -> pd.DataFrame:
    """
    Run differential expression analysis using PyDESeq2.

    Parameters
    ----------
    counts : pd.DataFrame
        Raw integer count matrix (genes × samples).
    metadata : pd.DataFrame
        Sample annotation table. Index must match columns of counts.
    design : str
        Design formula, e.g. '~condition' or 'condition'.
    contrast : list of str, optional
        Three-element list [factor, numerator, denominator],
        e.g. ['condition', 'treated', 'control'].
        If None, uses the last level of the design factor.
    lfc_threshold : float
        Log2 fold-change threshold for hypothesis testing (default 0).
    alpha : float
        Significance level for independent filtering (default 0.05).
    n_cpus : int
        Number of CPUs for parallel computation (default 2).
    quiet : bool
        Suppress progress output (default False).

    Returns
    -------
    pd.DataFrame
        Results table with columns: baseMean, log2FoldChange, lfcSE,
        stat, pvalue, padj. Index is gene names.

    Examples
    --------
    >>> results = run_deseq2(counts, meta, "~condition",
    ...                      contrast=["condition", "treated", "control"])
    >>> sig = results[(results.padj < 0.05) & (results.log2FoldChange.abs() > 1)]
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    design_col = design.replace("~", "").strip()

    # Align samples
    common = counts.columns.intersection(metadata.index)
    counts_sub = counts[common].T  # samples × genes
    meta_sub = metadata.loc[common]

    if not quiet:
        print(f"[DESeq2] Running on {counts_sub.shape[0]} samples, "
              f"{counts_sub.shape[1]} genes...")

    dds = DeseqDataSet(
        counts=counts_sub,
        metadata=meta_sub,
        design_factors=design_col,
        refit_cooks=True,
        n_cpus=n_cpus,
    )
    dds.deseq2()

    stat_kwargs = dict(alpha=alpha, lfc_null=lfc_threshold, n_cpus=n_cpus)
    if contrast:
        stat_res = DeseqStats(dds, contrast=contrast, **stat_kwargs)
    else:
        stat_res = DeseqStats(dds, **stat_kwargs)

    stat_res.summary()
    results = stat_res.results_df.copy()
    results.index.name = "gene"

    if not quiet:
        n_sig = (results["padj"] < alpha).sum()
        print(f"[DESeq2] Significant genes (padj < {alpha}): {n_sig}")

    return results


def run_edger_like(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    design: str,
    contrast: Optional[Tuple[str, str]] = None,
    alpha: float = 0.05,
    quiet: bool = False,
) -> pd.DataFrame:
    """
    Approximate edgeR quasi-likelihood pipeline using statsmodels NB GLM.

    This implementation mirrors the edgeR QLF test logic:
    1. TMM normalization of library sizes
    2. Negative Binomial GLM fitting per gene
    3. Quasi-likelihood F-test for differential expression

    Parameters
    ----------
    counts : pd.DataFrame
        Raw integer count matrix (genes × samples).
    metadata : pd.DataFrame
        Sample annotation table. Index must match columns of counts.
    design : str
        Design factor column name (without '~').
    contrast : tuple of str, optional
        (numerator_level, denominator_level), e.g. ('treated', 'control').
    alpha : float
        Significance level for BH correction (default 0.05).
    quiet : bool
        Suppress progress output.

    Returns
    -------
    pd.DataFrame
        Results with columns: baseMean, log2FoldChange, lfcSE, stat,
        pvalue, padj. Index is gene names.
    """
    import statsmodels.api as sm
    from statsmodels.stats.multitest import multipletests
    from .normalization import _normalize_tmm

    design_col = design.replace("~", "").strip()

    common = counts.columns.intersection(metadata.index)
    counts_sub = counts[common]
    meta_sub = metadata.loc[common]

    # TMM normalization factors for offset
    norm = _normalize_tmm(counts_sub)
    lib_sizes = counts_sub.sum(axis=0)
    norm_factors = (norm.sum(axis=0) / lib_sizes).fillna(1.0)
    offsets = np.log(lib_sizes * norm_factors)

    # Build design matrix
    condition = meta_sub[design_col].astype("category")
    if contrast:
        num_level, denom_level = contrast
        condition = condition.cat.reorder_categories(
            [denom_level] + [c for c in condition.cat.categories if c != denom_level],
            ordered=False,
        )
    dmat = pd.get_dummies(condition, drop_first=True).astype(float)
    dmat.insert(0, "Intercept", 1.0)

    results_list = []
    genes = counts_sub.index.tolist()

    if not quiet:
        print(f"[edgeR-like] Fitting NB GLM for {len(genes)} genes...")

    for gene in genes:
        y = counts_sub.loc[gene].values.astype(float)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                model = sm.NegativeBinomial(
                    y, dmat.values,
                    offset=offsets.values,
                    loglike_method="nb2",
                )
                fit = model.fit(disp=False, maxiter=100)
            coef_idx = 1  # first non-intercept coefficient
            lfc = fit.params[coef_idx] / np.log(2)
            se = fit.bse[coef_idx] / np.log(2)
            stat = fit.tvalues[coef_idx]
            pval = fit.pvalues[coef_idx]
        except Exception:
            lfc, se, stat, pval = np.nan, np.nan, np.nan, np.nan

        results_list.append({
            "gene": gene,
            "baseMean": float(y.mean()),
            "log2FoldChange": lfc,
            "lfcSE": se,
            "stat": stat,
            "pvalue": pval,
        })

    res_df = pd.DataFrame(results_list).set_index("gene")

    # BH correction
    valid = res_df["pvalue"].notna()
    padj = np.full(len(res_df), np.nan)
    _, padj_valid, _, _ = multipletests(
        res_df.loc[valid, "pvalue"], method="fdr_bh"
    )
    padj[valid.values] = padj_valid
    res_df["padj"] = padj

    if not quiet:
        n_sig = (res_df["padj"] < alpha).sum()
        print(f"[edgeR-like] Significant genes (padj < {alpha}): {n_sig}")

    return res_df
