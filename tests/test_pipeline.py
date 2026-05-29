"""
Integration tests for RNAseqPipeline.
Uses small synthetic data to test the full pipeline without heavy computation.
"""

import pytest
import numpy as np
import pandas as pd
import os
import tempfile
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))


def _make_synthetic_data(n_genes=200, n_samples=6, seed=42):
    """Create synthetic count data with a known DE signal."""
    rng = np.random.default_rng(seed)
    # Base counts
    counts = rng.negative_binomial(20, 0.5, size=(n_genes, n_samples)).astype(int)
    # Inject DE signal into first 20 genes
    counts[:10, 3:] = (counts[:10, 3:] * 4).astype(int)   # upregulated
    counts[10:20, 3:] = (counts[10:20, 3:] // 4 + 1).astype(int)  # downregulated

    genes = [f"GENE{i:04d}" for i in range(n_genes)]
    samples = [f"S{i:02d}" for i in range(n_samples)]
    counts_df = pd.DataFrame(counts, index=genes, columns=samples)
    meta_df = pd.DataFrame(
        {"condition": ["control"] * 3 + ["treated"] * 3},
        index=samples,
    )
    return counts_df, meta_df


@pytest.fixture
def data_dir(tmp_path):
    """Write synthetic data to temp CSV files and return paths."""
    counts, meta = _make_synthetic_data()
    counts_path = str(tmp_path / "counts.csv")
    meta_path = str(tmp_path / "metadata.csv")
    counts.to_csv(counts_path)
    meta.to_csv(meta_path)
    return counts_path, meta_path, str(tmp_path / "results")


# ── Pipeline init tests ───────────────────────────────────────────────────────

def test_pipeline_loads_data(data_dir):
    from rnaseq_toolkit.pipeline import RNAseqPipeline
    counts_path, meta_path, out_dir = data_dir
    pipe = RNAseqPipeline(counts_path, meta_path, "~condition", out_dir)
    assert pipe.counts.shape[1] == 6
    assert len(pipe.metadata) == 6


def test_pipeline_filters_low_counts(data_dir):
    from rnaseq_toolkit.pipeline import RNAseqPipeline
    counts_path, meta_path, out_dir = data_dir
    pipe = RNAseqPipeline(counts_path, meta_path, "~condition", out_dir,
                          min_count=5, min_samples=2)
    assert len(pipe.counts) <= 200


def test_pipeline_mismatched_samples_raises(tmp_path):
    from rnaseq_toolkit.pipeline import RNAseqPipeline
    counts = pd.DataFrame(
        np.ones((10, 3), dtype=int),
        index=[f"G{i}" for i in range(10)],
        columns=["A", "B", "C"],
    )
    meta = pd.DataFrame(
        {"condition": ["x", "y"]},
        index=["X", "Y"],  # different names
    )
    counts.to_csv(str(tmp_path / "c.csv"))
    meta.to_csv(str(tmp_path / "m.csv"))
    with pytest.raises(ValueError, match="No common samples"):
        RNAseqPipeline(str(tmp_path / "c.csv"), str(tmp_path / "m.csv"),
                       "~condition", str(tmp_path / "out"))


# ── Normalization tests ───────────────────────────────────────────────────────

def test_pipeline_normalization_cpm(data_dir):
    from rnaseq_toolkit.pipeline import RNAseqPipeline
    counts_path, meta_path, out_dir = data_dir
    pipe = RNAseqPipeline(counts_path, meta_path, "~condition", out_dir)
    norm = pipe.run_normalization(method="cpm")
    assert norm.shape == pipe.counts.shape
    # Check output file created
    assert os.path.exists(os.path.join(out_dir, "normalized_counts.csv"))


# ── DEA tests (lightweight — uses small data) ─────────────────────────────────

def test_pipeline_dea_deseq2_runs(data_dir):
    """Test that DESeq2 pipeline runs end-to-end on small data."""
    from rnaseq_toolkit.pipeline import RNAseqPipeline
    counts_path, meta_path, out_dir = data_dir
    pipe = RNAseqPipeline(counts_path, meta_path, "~condition", out_dir)
    pipe.run_normalization(method="cpm")
    results = pipe.run_dea(
        method="deseq2",
        contrast=["condition", "treated", "control"],
    )
    assert isinstance(results, pd.DataFrame)
    assert "log2FoldChange" in results.columns
    assert "padj" in results.columns
    assert len(results) > 0


def test_pipeline_results_file_created(data_dir):
    from rnaseq_toolkit.pipeline import RNAseqPipeline
    counts_path, meta_path, out_dir = data_dir
    pipe = RNAseqPipeline(counts_path, meta_path, "~condition", out_dir)
    pipe.run_normalization(method="cpm")
    pipe.run_dea(method="deseq2",
                 contrast=["condition", "treated", "control"])
    assert os.path.exists(os.path.join(out_dir, "deseq2_results.csv"))
