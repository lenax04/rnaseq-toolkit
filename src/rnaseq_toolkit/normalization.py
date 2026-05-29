"""
normalization.py — Unified normalization interface for RNA-seq count data.

Supported methods:
  - 'deseq2'  : Median-of-ratios (via PyDESeq2 size factors)
  - 'tmm'     : Trimmed Mean of M-values (approximated in Python)
  - 'tpm'     : Transcripts Per Million (requires gene lengths)
  - 'rpkm'    : Reads Per Kilobase per Million (requires gene lengths)
  - 'cpm'     : Counts Per Million
  - 'vst'     : Variance Stabilizing Transformation (via PyDESeq2)
  - 'rlog'    : Regularized log (approximated)
"""

import numpy as np
import pandas as pd
from typing import Optional, Literal


NormMethod = Literal["deseq2", "tmm", "tpm", "rpkm", "cpm", "vst", "rlog"]


def normalize_counts(
    counts: pd.DataFrame,
    method: NormMethod = "deseq2",
    metadata: Optional[pd.DataFrame] = None,
    design: Optional[str] = None,
    gene_lengths: Optional[pd.Series] = None,
) -> pd.DataFrame:
    """
    Normalize a raw count matrix using the specified method.

    Parameters
    ----------
    counts : pd.DataFrame
        Raw count matrix with genes as rows and samples as columns.
        All values must be non-negative integers.
    method : str
        Normalization method. One of: 'deseq2', 'tmm', 'tpm', 'rpkm',
        'cpm', 'vst', 'rlog'. Default: 'deseq2'.
    metadata : pd.DataFrame, optional
        Sample metadata (required for 'deseq2' and 'vst').
    design : str, optional
        Design formula string, e.g. '~condition' (required for 'deseq2'/'vst').
    gene_lengths : pd.Series, optional
        Gene lengths in base pairs (required for 'tpm' and 'rpkm').

    Returns
    -------
    pd.DataFrame
        Normalized count matrix with same shape as input.

    Examples
    --------
    >>> import pandas as pd
    >>> from rnaseq_toolkit import normalize_counts
    >>> counts = pd.read_csv("counts.csv", index_col=0)
    >>> meta = pd.read_csv("metadata.csv", index_col=0)
    >>> norm = normalize_counts(counts, method="deseq2",
    ...                         metadata=meta, design="~condition")
    """
    method = method.lower()

    if method == "cpm":
        return _normalize_cpm(counts)
    elif method == "tpm":
        if gene_lengths is None:
            raise ValueError("gene_lengths required for TPM normalization.")
        return _normalize_tpm(counts, gene_lengths)
    elif method == "rpkm":
        if gene_lengths is None:
            raise ValueError("gene_lengths required for RPKM normalization.")
        return _normalize_rpkm(counts, gene_lengths)
    elif method == "tmm":
        return _normalize_tmm(counts)
    elif method in ("deseq2", "vst", "rlog"):
        if metadata is None or design is None:
            raise ValueError(
                f"metadata and design are required for '{method}' normalization."
            )
        return _normalize_deseq2_family(counts, metadata, design, method)
    else:
        raise ValueError(
            f"Unknown normalization method: '{method}'. "
            f"Choose from: deseq2, tmm, tpm, rpkm, cpm, vst, rlog."
        )


# ── Private helpers ────────────────────────────────────────────────────────────

def _normalize_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    """Counts Per Million."""
    lib_sizes = counts.sum(axis=0)
    return counts.div(lib_sizes, axis=1) * 1e6


def _normalize_tpm(counts: pd.DataFrame, gene_lengths: pd.Series) -> pd.DataFrame:
    """Transcripts Per Million."""
    # Align gene lengths with count matrix rows
    lengths_kb = gene_lengths.reindex(counts.index).fillna(1000) / 1000.0
    rpk = counts.div(lengths_kb, axis=0)
    scaling = rpk.sum(axis=0) / 1e6
    return rpk.div(scaling, axis=1)


def _normalize_rpkm(counts: pd.DataFrame, gene_lengths: pd.Series) -> pd.DataFrame:
    """Reads Per Kilobase per Million."""
    lib_sizes = counts.sum(axis=0) / 1e6
    lengths_kb = gene_lengths.reindex(counts.index).fillna(1000) / 1000.0
    cpm = counts.div(lib_sizes, axis=1)
    return cpm.div(lengths_kb, axis=0)


def _normalize_tmm(counts: pd.DataFrame) -> pd.DataFrame:
    """
    Approximate TMM normalization (Robinson & Oshlack, 2010).
    Computes trimmed mean of M-values relative to a reference sample,
    then scales library sizes accordingly.
    """
    # Use the sample with 75th-percentile count as reference
    ref_idx = int(np.argmin(np.abs(counts.sum(axis=0).values -
                                   np.percentile(counts.sum(axis=0).values, 75))))
    ref = counts.iloc[:, ref_idx].astype(float)
    ref = ref.replace(0, np.nan)

    norm_factors = []
    for col in counts.columns:
        sample = counts[col].astype(float).replace(0, np.nan)
        # M-values (log2 fold-change)
        m = np.log2(sample / ref).dropna()
        # A-values (average log expression)
        a = (np.log2(sample) + np.log2(ref)).dropna() / 2
        common = m.index.intersection(a.index)
        m, a = m[common], a[common]
        # Trim 30% from M and 5% from A
        m_trim = _trim_series(m, 0.30)
        a_trim = _trim_series(a, 0.05)
        shared = m_trim.index.intersection(a_trim.index)
        if len(shared) == 0:
            norm_factors.append(1.0)
        else:
            norm_factors.append(2 ** m_trim[shared].mean())

    norm_factors = np.array(norm_factors)
    norm_factors /= np.exp(np.mean(np.log(norm_factors)))  # normalize to geometric mean

    lib_sizes = counts.sum(axis=0)
    effective_lib = lib_sizes * norm_factors
    return counts.div(effective_lib, axis=1) * 1e6


def _trim_series(s: pd.Series, proportion: float) -> pd.Series:
    n = len(s)
    lo = int(np.floor(n * proportion))
    hi = int(np.ceil(n * (1 - proportion)))
    sorted_idx = s.sort_values().index
    return s[sorted_idx[lo:hi]]


def _normalize_deseq2_family(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    design: str,
    method: str,
) -> pd.DataFrame:
    """DESeq2 median-of-ratios, VST, or rlog via PyDESeq2."""
    from pydeseq2.dds import DeseqDataSet

    design_col = design.replace("~", "").strip()
    # PyDESeq2 expects samples as rows, genes as columns
    counts_t = counts.T.copy()
    counts_t = counts_t.loc[metadata.index.intersection(counts_t.index)]
    meta_aligned = metadata.loc[counts_t.index]

    dds = DeseqDataSet(
        counts=counts_t,
        metadata=meta_aligned,
        design_factors=design_col,
        refit_cooks=True,
        n_cpus=2,
    )
    dds.fit_size_factors()

    if method == "deseq2":
        # Return size-factor-normalized counts (samples × genes → genes × samples)
        sf = dds.obsm["size_factors"] if "size_factors" in dds.obsm else dds.size_factors
        norm = counts_t.div(pd.Series(sf, index=counts_t.index), axis=0)
        return norm.T  # back to genes × samples

    elif method == "vst":
        dds.vst(use_design=False)
        vst_data = dds.layers["vst_counts"]
        return pd.DataFrame(vst_data.T, index=counts.index, columns=counts.columns)

    elif method == "rlog":
        # Approximate rlog: log2(norm_counts + 0.5)
        sf = dds.size_factors
        norm = counts_t.div(pd.Series(sf, index=counts_t.index), axis=0)
        rlog = np.log2(norm + 0.5)
        return rlog.T
