"""
qc_split_separation_diagnostic.py
=================================
Label-free split-quality diagnostic for the HSV-1 leakage ladder.

An AUROC collapse under a grouped split tells you that leakage existed, but not which
channel was open. This script measures the channel directly and without labels:

  For every held-out TEST position in a fold, the sequence distance to the nearest
  TRAIN position. If that distance is 1-2 residues, the sequence-neighbourhood channel
  is still open and any model with local sequence sensitivity can interpolate across
  the split, however the groups were defined.

Schemes compared (all pre-computed in data/processed/primary_3d_spatial_folds.tsv):
  cv_fold_id          position-grouped 5-fold
  seq_block_id        contiguous sequence block 5-fold
  spatial_cluster_id  3D k-means on C-alpha coordinates, 5 clusters

Also reports the cluster-agreement between the 3D and contiguous-block schemes, to test
whether a spatial split is an independent control or a noisy restatement of a sequence
split.

CRITICAL CONSTRAINT: strictly zero resistance labels, fold-changes or phenotypes are read.
Only positions and fold/cluster assignments are used.
"""

import json
import os

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_mutual_info_score, normalized_mutual_info_score

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

SCHEMES = [("position-grouped", "cv_fold_id"),
           ("contiguous-block", "seq_block_id"),
           ("3D k-means", "spatial_cluster_id")]


def nearest_train_distances(sub, col):
    """Min |pos_test - pos_train| for every test position, pooled over folds."""
    out = []
    for f in sorted(sub[col].unique()):
        te = sub.loc[sub[col] == f, "pos"].values
        tr = sub.loc[sub[col] != f, "pos"].values
        if len(te) == 0 or len(tr) == 0:
            continue
        out += [int(np.min(np.abs(tr - p))) for p in te]
    return np.asarray(out)


def main():
    path = os.path.join(ROOT, "data", "processed", "primary_3d_spatial_folds.tsv")
    df = pd.read_csv(path, sep="\t")
    print("=" * 78)
    print("LABEL-FREE SPLIT SEPARATION DIAGNOSTIC")
    print("=" * 78)
    print("source: %s" % os.path.relpath(path, ROOT))
    print("rows: %d   rows with 3D coordinates: %d" % (len(df), int(df.has_3d_coord.sum())))
    print("")
    print("Sequence distance from each TEST position to its nearest TRAIN position.")
    print("Low median => sequence-neighbourhood channel still open => split does not")
    print("control local-sequence interpolation, regardless of how groups were defined.")

    out = {"scope": "LABEL-FREE: no resistance labels, fold-changes or phenotypes read",
           "source": os.path.relpath(path, ROOT),
           "n_rows": int(len(df)),
           "n_rows_with_3d": int(df.has_3d_coord.sum()),
           "metric": ("min |pos_test - pos_train| per held-out position, pooled over "
                      "the 5 folds of each scheme; unique positions only"),
           "targets": {}}

    for (gene, drug), s in df.groupby(["gene", "drug"]):
        s = s[s.has_3d_coord == 1].drop_duplicates("pos")
        if len(s) < 20:
            continue
        key = "%s-%s" % (gene, drug)
        print("")
        print("== %s   unique positions = %d" % (key, len(s)))
        print("   %-18s %8s %7s %8s %9s %9s" %
              ("scheme", "median", "p25", "mean", "frac<=2", "frac<=5"))
        rec = {"n_unique_positions": int(len(s)), "schemes": {}}
        for name, col in SCHEMES:
            a = nearest_train_distances(s, col)
            print("   %-18s %8.1f %7.1f %8.2f %9.3f %9.3f" %
                  (name, np.median(a), np.percentile(a, 25), a.mean(),
                   (a <= 2).mean(), (a <= 5).mean()))
            rec["schemes"][name] = {
                "column": col,
                "median_nearest_train_distance": float(np.median(a)),
                "p25": float(np.percentile(a, 25)),
                "mean": round(float(a.mean()), 3),
                "frac_within_2_residues": round(float((a <= 2).mean()), 4),
                "frac_within_5_residues": round(float((a <= 5).mean()), 4)}

        up = s
        rec["spatial_vs_contiguous_agreement"] = {
            "normalized_mutual_info": round(float(normalized_mutual_info_score(
                up.spatial_cluster_id, up.seq_block_id)), 4),
            "adjusted_mutual_info": round(float(adjusted_mutual_info_score(
                up.spatial_cluster_id, up.seq_block_id)), 4)}
        print("   spatial vs contiguous-block cluster agreement: NMI=%.3f  AMI=%.3f"
              % (rec["spatial_vs_contiguous_agreement"]["normalized_mutual_info"],
                 rec["spatial_vs_contiguous_agreement"]["adjusted_mutual_info"]))
        out["targets"][key] = rec

    log_path = os.path.join(ROOT, "logs", "split_separation_diagnostic.json")
    with open(log_path, "w") as f:
        json.dump(out, f, indent=1)
    print("")
    print("Wrote %s" % os.path.relpath(log_path, ROOT))
    print("")
    print("NO LABELS, FOLD-CHANGES OR PHENOTYPES WERE READ.")


if __name__ == "__main__":
    main()
