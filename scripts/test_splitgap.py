"""
test_splitgap.py
================
Unit test suite for the splitgap diagnostic library and CLI.
"""

import sys
import os
import subprocess
import numpy as np
import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from splitgap.core import (
    compute_nearest_train_distances,
    audit_cv_splits,
    partition_modulo,
    partition_contiguous_blocks,
)
from splitgap.recommend import recommend_block_partition


def test_nearest_train_distances():
    test_pos = [10, 20, 30]
    train_pos = [8, 22, 50]
    # |10 - 8| = 2, |20 - 22| = 2, |30 - 22| = 8
    expected = np.array([2, 2, 8])
    res = compute_nearest_train_distances(test_pos, train_pos)
    np.testing.assert_array_equal(res, expected)


def test_input_validation():
    with pytest.raises(ValueError):
        compute_nearest_train_distances([1, 2], [])

    with pytest.raises(ValueError):
        audit_cv_splits([1, 2, 3], [0, 1])  # Length mismatch

    with pytest.raises(ValueError):
        audit_cv_splits([1, 2, 3], [0, 0, 0])  # Only 1 fold


def test_modulo_leakage_detection():
    # Dense consecutive positions 1..100
    positions = np.arange(1, 101)
    mod_folds = partition_modulo(positions, k_folds=5)
    audit = audit_cv_splits(positions, mod_folds)

    # In consecutive positions, modulo-5 leaves test position p with train positions p-1 and p+1, so gap = 1 aa!
    assert audit["median_gap"] == 1.0
    assert audit["pct_le_1aa"] == 100.0
    assert audit["risk_tier"] == "CRITICAL_SPATIAL_LEAKAGE"


def test_contiguous_block_separation():
    positions = np.arange(1, 101)
    blk_folds = partition_contiguous_blocks(positions, k_folds=5)
    audit = audit_cv_splits(positions, blk_folds)

    # In contiguous blocks [1..20], [21..40], ...
    # Test residues at block edges have gap 1, but interior residues have gaps up to 10 aa!
    assert audit["median_gap"] >= 5.0
    assert audit["pct_le_2aa"] < 25.0


def test_recommend_block_partition():
    # Sparsely distributed positions
    positions = [10, 25, 40, 55, 70, 85, 100, 115, 130, 145]
    rec = recommend_block_partition(positions, target_median_gap=15.0, target_p10_gap=5.0)
    assert rec["status"] in ["OPTIMAL_SEPARATION_ACHIEVED", "PARTIAL_SEPARATION_ACHIEVED"]
    assert rec["recommended_k_folds"] >= 3


def test_cli_execution():
    python_bin = sys.executable
    # Test CLI audit with comma-separated list of consecutive integers
    cmd = [python_bin, "-m", "splitgap.cli", "audit", "-p", "10,11,12,13,14,15,16,17,18,19,20", "--json"]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=ROOT)
    assert proc.returncode == 0, f"CLI error: {proc.stderr}"
    assert "modulo_split" in proc.stdout
    assert "contiguous_block_split" in proc.stdout

    # Test CLI recommend
    cmd_rec = [python_bin, "-m", "splitgap.cli", "recommend", "-p", "10,20,30,40,50,60,70,80,90,100", "--json"]
    proc_rec = subprocess.run(cmd_rec, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=ROOT)
    assert proc_rec.returncode == 0, f"CLI error: {proc_rec.stderr}"
    assert "recommended_k_folds" in proc_rec.stdout


if __name__ == "__main__":
    test_nearest_train_distances()
    test_input_validation()
    test_modulo_leakage_detection()
    test_contiguous_block_separation()
    test_recommend_block_partition()
    test_cli_execution()
    print("ALL SPLITGAP UNIT TESTS PASSED DETERMINISTICALLY (100%).")
