"""
splitgap.core
=============
Core vectorized algorithms for auditing Sequence-Neighbour Gap spatial leakage
in biological cross-validation splits.
"""

from typing import Dict, List, Union, Tuple, Any
import numpy as np


def compute_nearest_train_distances(
    test_positions: Union[List[int], np.ndarray],
    train_positions: Union[List[int], np.ndarray]
) -> np.ndarray:
    """Compute the absolute sequence distance from each test residue to the nearest training residue.

    Parameters
    ----------
    test_positions : array-like of shape (N_test,)
        Residue indices in the evaluation / test set.
    train_positions : array-like of shape (N_train,)
        Residue indices in the training set.

    Returns
    -------
    np.ndarray of shape (N_test,)
        Vector of minimum sequence distances: min_{q in train} |p - q|.
    """
    test_arr = np.asarray(test_positions, dtype=np.int64)
    train_arr = np.asarray(train_positions, dtype=np.int64)

    if len(test_arr) == 0:
        return np.array([], dtype=np.int64)
    if len(train_arr) == 0:
        raise ValueError("train_positions cannot be empty.")

    # Vectorized broadcasting: (N_test, 1) - (1, N_train)
    diffs = np.abs(test_arr[:, np.newaxis] - train_arr[np.newaxis, :])
    return np.min(diffs, axis=1)


def audit_cv_splits(
    positions: Union[List[int], np.ndarray],
    fold_assignments: Union[List[int], np.ndarray]
) -> Dict[str, Any]:
    """Audit an existing cross-validation assignment for Sequence-Neighbour Gap leakage.

    Parameters
    ----------
    positions : array-like of shape (N,)
        Residue sequence indices.
    fold_assignments : array-like of shape (N,)
        Integer fold identifiers (e.g. 0..K-1) corresponding to each position.

    Returns
    -------
    dict
        Comprehensive summary metrics including percentiles, leakage fractions,
        and diagnostic risk tier.
    """
    pos_arr = np.asarray(positions, dtype=np.int64)
    folds = np.asarray(fold_assignments, dtype=np.int64)

    if len(pos_arr) != len(folds):
        raise ValueError(f"Length mismatch: {len(pos_arr)} positions vs {len(folds)} fold assignments.")

    unique_folds = np.unique(folds)
    if len(unique_folds) < 2:
        raise ValueError("Must provide at least 2 distinct folds to audit cross-validation.")

    all_gaps: List[int] = []
    per_fold_metrics: Dict[int, Dict[str, float]] = {}

    for k in unique_folds:
        test_mask = (folds == k)
        test_pos = pos_arr[test_mask]
        train_pos = pos_arr[~test_mask]

        if len(test_pos) == 0:
            continue

        gaps = compute_nearest_train_distances(test_pos, train_pos)
        all_gaps.extend(gaps.tolist())

        per_fold_metrics[int(k)] = {
            "n_test": int(len(test_pos)),
            "median_gap": float(np.median(gaps)),
            "pct_le_2aa": float(np.mean(gaps <= 2) * 100),
            "pct_le_5aa": float(np.mean(gaps <= 5) * 100),
        }

    gaps_arr = np.array(all_gaps, dtype=np.int64)

    pct_le_1aa = float(np.mean(gaps_arr <= 1) * 100)
    pct_le_2aa = float(np.mean(gaps_arr <= 2) * 100)
    pct_le_3aa = float(np.mean(gaps_arr <= 3) * 100)
    pct_le_5aa = float(np.mean(gaps_arr <= 5) * 100)
    pct_le_10aa = float(np.mean(gaps_arr <= 10) * 100)

    median_gap = float(np.median(gaps_arr))

    # Risk tier determination
    if pct_le_2aa >= 50.0 or median_gap <= 2.0:
        risk_tier = "CRITICAL_SPATIAL_LEAKAGE"
        recommendation = "Reject split. >50% of test residues sit within 2 aa of training data."
    elif pct_le_5aa >= 50.0 or median_gap <= 5.0:
        risk_tier = "MODERATE_SPATIAL_LEAKAGE"
        recommendation = "Warning: Moderate proximity. >50% of test residues sit within 5 aa of training data."
    else:
        risk_tier = "SPATIALLY_SEPARATED"
        recommendation = "Acceptable spatial separation. Generalization evaluated across distinct sequence blocks."

    return {
        "n_total_positions": int(len(pos_arr)),
        "n_folds": int(len(unique_folds)),
        "median_gap": median_gap,
        "mean_gap": float(np.mean(gaps_arr)),
        "std_gap": float(np.std(gaps_arr)),
        "p10": float(np.percentile(gaps_arr, 10)),
        "p25": float(np.percentile(gaps_arr, 25)),
        "p50": median_gap,
        "p75": float(np.percentile(gaps_arr, 75)),
        "p90": float(np.percentile(gaps_arr, 90)),
        "pct_le_1aa": pct_le_1aa,
        "pct_le_2aa": pct_le_2aa,
        "pct_le_3aa": pct_le_3aa,
        "pct_le_5aa": pct_le_5aa,
        "pct_le_10aa": pct_le_10aa,
        "risk_tier": risk_tier,
        "recommendation": recommendation,
        "per_fold": per_fold_metrics,
        "gaps": gaps_arr.tolist()
    }


def partition_modulo(positions: Union[List[int], np.ndarray], k_folds: int = 5) -> np.ndarray:
    """Generate Modulo-K fold assignments: fold = pos % k_folds."""
    pos_arr = np.asarray(positions, dtype=np.int64)
    return pos_arr % k_folds


def partition_contiguous_blocks(positions: Union[List[int], np.ndarray], k_folds: int = 5) -> np.ndarray:
    """Generate Contiguous-Block fold assignments after sorting positions."""
    pos_arr = np.asarray(positions, dtype=np.int64)
    sorted_indices = np.argsort(pos_arr)
    n = len(pos_arr)
    block_size = int(np.ceil(n / k_folds))

    folds = np.zeros(n, dtype=np.int64)
    for k in range(k_folds):
        start = k * block_size
        end = min((k + 1) * block_size, n)
        folds[sorted_indices[start:end]] = k

    return folds

