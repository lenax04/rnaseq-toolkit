"""
conftest.py — Shared pytest configuration and fixtures.
"""
import pytest
import numpy as np
import pandas as pd


@pytest.fixture(scope="session")
def synthetic_counts():
    """Session-scoped synthetic count matrix (500 genes × 8 samples)."""
    rng = np.random.default_rng(0)
    n_genes, n_samples = 500, 8
    counts = rng.negative_binomial(15, 0.5, size=(n_genes, n_samples)).astype(int)
    genes = [f"GENE{i:05d}" for i in range(n_genes)]
    samples = [f"SAMPLE_{i:02d}" for i in range(n_samples)]
    return pd.DataFrame(counts, index=genes, columns=samples)


@pytest.fixture(scope="session")
def synthetic_metadata():
    """Session-scoped metadata for 8 samples."""
    return pd.DataFrame(
        {
            "condition": ["control"] * 4 + ["treated"] * 4,
            "batch": ["A", "A", "B", "B"] * 2,
        },
        index=[f"SAMPLE_{i:02d}" for i in range(8)],
    )
