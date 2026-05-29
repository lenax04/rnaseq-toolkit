"""
Tests for normalization module.
"""

import pytest
import numpy as np
import pandas as pd
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from rnaseq_toolkit.normalization import normalize_counts


# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture
def small_counts():
    """Small synthetic count matrix: 100 genes × 6 samples."""
    np.random.seed(42)
    n_genes, n_samples = 100, 6
    counts = np.random.negative_binomial(10, 0.5, size=(n_genes, n_samples))
    genes = [f"GENE{i:04d}" for i in range(n_genes)]
    samples = [f"S{i:02d}" for i in range(n_samples)]
    return pd.DataFrame(counts, index=genes, columns=samples)


@pytest.fixture
def small_metadata():
    """Metadata for 6 samples: 3 control, 3 treated."""
    return pd.DataFrame(
        {"condition": ["control"] * 3 + ["treated"] * 3},
        index=[f"S{i:02d}" for i in range(6)],
    )


# ── CPM tests ─────────────────────────────────────────────────────────────────

def test_cpm_shape(small_counts):
    norm = normalize_counts(small_counts, method="cpm")
    assert norm.shape == small_counts.shape


def test_cpm_column_sums(small_counts):
    norm = normalize_counts(small_counts, method="cpm")
    # Each column should sum to ~1e6
    for col in norm.columns:
        assert abs(norm[col].sum() - 1e6) < 1.0, f"Column {col} sum != 1e6"


def test_cpm_non_negative(small_counts):
    norm = normalize_counts(small_counts, method="cpm")
    assert (norm >= 0).all().all()


# ── TMM tests ─────────────────────────────────────────────────────────────────

def test_tmm_shape(small_counts):
    norm = normalize_counts(small_counts, method="tmm")
    assert norm.shape == small_counts.shape


def test_tmm_non_negative(small_counts):
    norm = normalize_counts(small_counts, method="tmm")
    assert (norm >= 0).all().all()


# ── TPM tests ─────────────────────────────────────────────────────────────────

def test_tpm_requires_lengths(small_counts):
    with pytest.raises(ValueError, match="gene_lengths"):
        normalize_counts(small_counts, method="tpm")


def test_tpm_column_sums(small_counts):
    lengths = pd.Series(
        np.random.randint(500, 5000, size=len(small_counts)),
        index=small_counts.index,
    )
    norm = normalize_counts(small_counts, method="tpm", gene_lengths=lengths)
    for col in norm.columns:
        assert abs(norm[col].sum() - 1e6) < 1.0, f"Column {col} sum != 1e6"


# ── Unknown method ─────────────────────────────────────────────────────────────

def test_unknown_method_raises(small_counts):
    with pytest.raises(ValueError, match="Unknown normalization method"):
        normalize_counts(small_counts, method="invalid_method")


# ── Index / column preservation ───────────────────────────────────────────────

def test_cpm_preserves_index(small_counts):
    norm = normalize_counts(small_counts, method="cpm")
    assert list(norm.index) == list(small_counts.index)
    assert list(norm.columns) == list(small_counts.columns)
