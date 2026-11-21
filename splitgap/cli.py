"""
splitgap.cli
============
Command-line interface for the splitgap diagnostic utility.
"""

import sys
import os
import argparse
import json
import pandas as pd
import numpy as np
from .core import (
    audit_cv_splits,
    partition_modulo,
    partition_contiguous_blocks,
)
from .recommend import recommend_block_partition


def _parse_positions_arg(pos_arg: str, col_name: str = "pos") -> np.ndarray:
    if os.path.exists(pos_arg):
        sep = "\t" if pos_arg.endswith((".tsv", ".tab")) else ","
        df = pd.read_csv(pos_arg, sep=sep)
        col = col_name if col_name in df.columns else ("position" if "position" in df.columns else None)
        if col is None:
            raise ValueError(f"Column '{col_name}' not found in {pos_arg}. Available: {list(df.columns)}")
        return df[col].dropna().astype(int).values
    else:
        # Try parsing comma-separated integers
        parts = [p.strip() for p in pos_arg.split(",") if p.strip()]
        return np.array([int(p) for p in parts], dtype=np.int64)


def main():
    parser = argparse.ArgumentParser(
        prog="splitgap",
        description="Audit Sequence-Neighbour Gap spatial leakage in biological cross-validation splits."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand: audit
    audit_parser = subparsers.add_parser("audit", help="Audit an existing or synthetic split for spatial leakage.")
    audit_parser.add_argument("--positions", "-p", required=True, help="File path (CSV/TSV) or comma-separated position list.")
    audit_parser.add_argument("--pos-col", default="pos", help="Column name containing residue indices (default: 'pos').")
    audit_parser.add_argument("--fold-col", default=None, help="Optional column name containing fold assignments.")
    audit_parser.add_argument("--strategy", choices=["modulo", "contiguous", "both"], default="both",
                              help="Partitioning strategy to test if fold-col is not provided.")
    audit_parser.add_argument("--k-folds", "-k", type=int, default=5, help="Number of folds (default: 5).")
    audit_parser.add_argument("--json", action="store_true", help="Output raw JSON results.")

    # Subcommand: recommend
    rec_parser = subparsers.add_parser("recommend", help="Recommend optimal contiguous block partitioning.")
    rec_parser.add_argument("--positions", "-p", required=True, help="File path (CSV/TSV) or comma-separated position list.")
    rec_parser.add_argument("--pos-col", default="pos", help="Column name containing residue indices (default: 'pos').")
    rec_parser.add_argument("--target-median", type=float, default=10.0, help="Target median sequence gap (default: 10.0 aa).")
    rec_parser.add_argument("--target-p10", type=float, default=5.0, help="Target 10th percentile sequence gap (default: 5.0 aa).")
    rec_parser.add_argument("--json", action="store_true", help="Output raw JSON results.")

    args = parser.parse_args()

    positions = _parse_positions_arg(args.positions, col_name=args.pos_col)

    if args.command == "audit":
        results = {}
        if args.fold_col is not None and os.path.exists(args.positions):
            sep = "\t" if args.positions.endswith((".tsv", ".tab")) else ","
            df = pd.read_csv(args.positions, sep=sep)
            folds = df[args.fold_col].values
            results["custom_split"] = audit_cv_splits(positions, folds)
        else:
            if args.strategy in ["modulo", "both"]:
                mod_folds = partition_modulo(positions, k_folds=args.k_folds)
                results["modulo_split"] = audit_cv_splits(positions, mod_folds)
            if args.strategy in ["contiguous", "both"]:
                blk_folds = partition_contiguous_blocks(positions, k_folds=args.k_folds)
                results["contiguous_block_split"] = audit_cv_splits(positions, blk_folds)

        if args.json:
            print(json.dumps(results, indent=2))
        else:
            print("=" * 80)
            print("SPLITGAP AUDIT REPORT")
            print("=" * 80)
            print(f"Total Unique Positions Evaluated: {len(np.unique(positions))}")
            for split_name, res in results.items():
                print(f"\n--- Stratum: {split_name.upper()} ---")
                print(f"  Risk Tier          : {res['risk_tier']}")
                print(f"  Median Gap         : {res['median_gap']:.1f} aa")
                print(f"  25th Percentile    : {res['p25']:.1f} aa")
                print(f"  10th Percentile    : {res['p10']:.1f} aa")
                print(f"  Residues <= 2 aa   : {res['pct_le_2aa']:.1f}%")
                print(f"  Residues <= 5 aa   : {res['pct_le_5aa']:.1f}%")
                print(f"  Recommendation     : {res['recommendation']}")
            print("=" * 80)

    elif args.command == "recommend":
        rec = recommend_block_partition(
            positions,
            target_median_gap=args.target_median,
            target_p10_gap=args.target_p10
        )
        if args.json:
            print(json.dumps(rec, indent=2))
        else:
            print("=" * 80)
            print("SPLITGAP RECOMMENDATION REPORT")
            print("=" * 80)
            print(f"Status              : {rec['status']}")
            print(f"Recommended K-Folds : {rec['recommended_k_folds']}")
            print(f"Block Size          : {rec['block_size_residues']} residues per fold")
            print(f"Achieved Median Gap : {rec['achieved_median_gap']:.1f} aa")
            print(f"Achieved P10 Gap    : {rec['achieved_p10_gap']:.1f} aa")
            print(f"Achieved <= 2 aa    : {rec['achieved_pct_le_2aa']:.1f}%")
            print(f"Note                : {rec['note']}")
            print("=" * 80)


if __name__ == "__main__":
    main()

