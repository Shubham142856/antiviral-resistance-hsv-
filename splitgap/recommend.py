"""
splitgap.recommend
==================
Algorithmic recommendation engine for designing biologically valid,
spatially separated cross-validation splits.
"""

from typing import Dict, List, Union, Tuple, Any
import numpy as np
from .core import audit_cv_splits, compute_nearest_train_distances


def recommend_block_partition(
    positions: Union[List[int], np.ndarray],
    target_median_gap: float = 10.0,
    target_p10_gap: float = 5.0,
    max_folds: int = 10,
    min_folds: int = 3
) -> Dict[str, Any]:
    """Find the optimal number of contiguous blocks and partition that satisfies
    the requested minimum spatial separation.

    Parameters
    ----------
    positions : array-like of shape (N,)
        Residue sequence indices.
    target_median_gap : float, default=10.0
        Desired minimum median sequence separation (in amino acids) to training set.
    target_p10_gap : float, default=5.0
        Desired minimum 10th percentile separation (in amino acids).
    max_folds : int, default=10
        Maximum number of folds to evaluate.
    min_folds : int, default=3
        Minimum number of folds to evaluate.

    Returns
    -------
    dict
        Recommended configuration, achieved gap metrics, and fold assignments.
    """
    pos_arr = np.asarray(sorted(set(positions)), dtype=np.int64)
    n = len(pos_arr)

    evaluations: List[Dict[str, Any]] = []
    best_candidate = None

    for k in range(min_folds, max_folds + 1):
        block_size = int(np.ceil(n / k))
        folds = np.zeros(n, dtype=np.int64)
        for fold_idx in range(k):
            start = fold_idx * block_size
            end = min((fold_idx + 1) * block_size, n)
            folds[start:end] = fold_idx

        audit = audit_cv_splits(pos_arr, folds)
        meets_criteria = (audit["median_gap"] >= target_median_gap) and (audit["p10"] >= target_p10_gap)

        cand = {
            "k_folds": k,
            "block_size_residues": block_size,
            "median_gap": audit["median_gap"],
            "p10": audit["p10"],
            "pct_le_2aa": audit["pct_le_2aa"],
            "pct_le_5aa": audit["pct_le_5aa"],
            "meets_criteria": meets_criteria,
            "audit": audit,
            "folds": folds.tolist()
        }
        evaluations.append(cand)

        if meets_criteria and best_candidate is None:
            best_candidate = cand

    # If no candidate strictly met criteria, choose candidate with highest median gap
    if best_candidate is None:
        best_candidate = max(evaluations, key=lambda x: x["median_gap"])
        status = "PARTIAL_SEPARATION_ACHIEVED"
        note = (
            f"Unable to achieve target p10 >= {target_p10_gap} aa due to tight pocket clustering. "
            f"Selected k={best_candidate['k_folds']} for maximum achievable separation."
        )
    else:
        status = "OPTIMAL_SEPARATION_ACHIEVED"
        note = f"Selected k={best_candidate['k_folds']} satisfying median >= {target_median_gap} and p10 >= {target_p10_gap}."

    return {
        "status": status,
        "note": note,
        "recommended_k_folds": best_candidate["k_folds"],
        "block_size_residues": best_candidate["block_size_residues"],
        "achieved_median_gap": best_candidate["median_gap"],
        "achieved_p10_gap": best_candidate["p10"],
        "achieved_pct_le_2aa": best_candidate["pct_le_2aa"],
        "achieved_pct_le_5aa": best_candidate["pct_le_5aa"],
        "all_evaluations": [
            {
                "k_folds": e["k_folds"],
                "median_gap": e["median_gap"],
                "p10": e["p10"],
                "pct_le_2aa": e["pct_le_2aa"],
                "pct_le_5aa": e["pct_le_5aa"],
                "meets_criteria": e["meets_criteria"]
            }
            for e in evaluations
        ]
    }
