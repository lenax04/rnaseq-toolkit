"""
Tests for visualization module.
"""

import pytest
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
import os
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from rnaseq_toolkit.visualization import plot_volcano, plot_ma, plot_pca


# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture
def mock_dea_results():
    """Synthetic DEA results DataFrame."""
    np.random.seed(42)
    n = 500
    lfc = np.random.normal(0, 2, n)
    pval = np.random.uniform(0, 1, n)
    padj = np.minimum(pval * n / np.arange(1, n + 1)[np.argsort(pval).argsort()], 1.0)
    base_mean = np.abs(np.random.normal(100, 50, n))
    genes = [f"GENE{i:04d}" for i in range(n)]
    return pd.DataFrame({
        "baseMean": base_mean,
        "log2FoldChange": lfc,
        "lfcSE": np.abs(np.random.normal(0.3, 0.1, n)),
        "stat": lfc / 0.3,
        "pvalue": pval,
        "padj": padj,
    }, index=genes)


@pytest.fixture
def mock_norm_counts():
    """Synthetic normalized count matrix."""
    np.random.seed(42)
    n_genes, n_samples = 200, 6
    counts = np.abs(np.random.normal(100, 30, (n_genes, n_samples)))
    genes = [f"GENE{i:04d}" for i in range(n_genes)]
    samples = [f"S{i:02d}" for i in range(n_samples)]
    return pd.DataFrame(counts, index=genes, columns=samples)


@pytest.fixture
def mock_metadata():
    return pd.DataFrame(
        {"condition": ["control"] * 3 + ["treated"] * 3},
        index=[f"S{i:02d}" for i in range(6)],
    )


# ── Volcano tests ─────────────────────────────────────────────────────────────

def test_volcano_returns_figure(mock_dea_results):
    fig = plot_volcano(mock_dea_results)
    assert isinstance(fig, plt.Figure)
    plt.close("all")


def test_volcano_saves_file(mock_dea_results):
    with tempfile.TemporaryDirectory() as tmpdir:
        out = os.path.join(tmpdir, "volcano.png")
        plot_volcano(mock_dea_results, output_path=out)
        assert os.path.exists(out)
        assert os.path.getsize(out) > 1000  # non-empty file
    plt.close("all")


def test_volcano_with_nan(mock_dea_results):
    """Should handle NaN values gracefully."""
    res = mock_dea_results.copy()
    res.loc[res.index[:50], "padj"] = np.nan
    fig = plot_volcano(res)
    assert isinstance(fig, plt.Figure)
    plt.close("all")


# ── MA plot tests ─────────────────────────────────────────────────────────────

def test_ma_returns_figure(mock_dea_results):
    fig = plot_ma(mock_dea_results)
    assert isinstance(fig, plt.Figure)
    plt.close("all")


def test_ma_saves_file(mock_dea_results):
    with tempfile.TemporaryDirectory() as tmpdir:
        out = os.path.join(tmpdir, "ma.png")
        plot_ma(mock_dea_results, output_path=out)
        assert os.path.exists(out)
    plt.close("all")


# ── PCA tests ─────────────────────────────────────────────────────────────────

def test_pca_returns_figure(mock_norm_counts, mock_metadata):
    fig = plot_pca(mock_norm_counts, mock_metadata, color_by="condition")
    assert isinstance(fig, plt.Figure)
    plt.close("all")


def test_pca_saves_file(mock_norm_counts, mock_metadata):
    with tempfile.TemporaryDirectory() as tmpdir:
        out = os.path.join(tmpdir, "pca.png")
        plot_pca(mock_norm_counts, mock_metadata, color_by="condition",
                 output_path=out)
        assert os.path.exists(out)
    plt.close("all")
