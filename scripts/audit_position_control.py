"""Reproduce the split audit on matched records and report within-fold AUROCs."""

import hashlib
import json
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import sklearn
from sklearn.metrics import roc_auc_score

from run_split_degradation_benchmark import PROTEIN_LENGTHS, SPLIT_SCHEMES, TARGETS, run_oof_cv

matplotlib.use("Agg")
from matplotlib import pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/processed/primary_3d_spatial_folds.tsv"


def evaluate(frame, column):
    if frame.groupby("pos")[column].nunique().max() != 1:
        raise ValueError(f"Position assigned to multiple folds in {column}")
    pooled, scores = run_oof_cv(frame, ["norm_pos"], column)
    training_prevalence = np.empty(len(frame))
    gaps = []
    folds = []
    for fold, group in frame.groupby(column):
        test_pos = group.pos.unique()
        train_pos = frame.loc[frame[column] != fold, "pos"].unique()
        gaps.extend(np.abs(test_pos[:, None] - train_pos).min(axis=1))
        training_prevalence[group.index] = frame.loc[frame[column] != fold, "target_binary"].mean()
        labels = group.target_binary
        auc = roc_auc_score(labels, scores[group.index]) if labels.nunique() == 2 else None
        folds.append(dict(fold=int(fold), records=len(group), positive=int(labels.sum()),
                          positions=len(test_pos), auc=auc))
    values = [row["auc"] for row in folds if row["auc"] is not None]
    positive = scores[frame.target_binary == 1]
    negative = scores[frame.target_binary == 0]
    same_fold = (frame.loc[frame.target_binary == 1, column].to_numpy()[:, None]
                 == frame.loc[frame.target_binary == 0, column].to_numpy())
    concordance = (positive[:, None] > negative) + 0.5 * (positive[:, None] == negative)
    pair_summary = {
        group: dict(count=int(mask.sum()), concordance=float(concordance[mask].mean()))
        for group, mask in [("within", same_fold), ("between", ~same_fold)] if mask.any()}
    if not np.isclose(concordance.mean(), pooled):
        raise ValueError("Pair concordance disagrees with pooled AUROC")
    return dict(pair_summary=pair_summary, records=len(frame), positions=len(gaps), median_gap=float(np.median(gaps)),
                prevalence_only_pooled=float(roc_auc_score(frame.target_binary, training_prevalence)),
                pooled=pooled, within_fold_mean=float(np.mean(values)), folds=folds), scores


def paired_position_interval(frame, grouped, blocked):
    """Conditional on fitted scores, resample whole positions in both schemes together."""
    sites, inverse = np.unique(frame.pos, return_inverse=True)
    rng = np.random.default_rng(42)
    differences = []
    for _ in range(2000):
        counts = np.bincount(rng.integers(len(sites), size=len(sites)), minlength=len(sites))
        weights = counts[inverse]
        if frame.target_binary[weights > 0].nunique() < 2:
            continue
        differences.append(roc_auc_score(frame.target_binary, grouped, sample_weight=weights)
                           - roc_auc_score(frame.target_binary, blocked, sample_weight=weights))
    return dict(delta=float(roc_auc_score(frame.target_binary, grouped)
                            - roc_auc_score(frame.target_binary, blocked)),
                ci95=np.quantile(differences, [0.025, 0.975]).tolist(),
                valid_draws=len(differences), requested_draws=2000, seed=42)


def figure(results):
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.6), sharey=True)
    for ax, (arm, cohorts) in zip(axes, results.items()):
        for x, (_, result) in enumerate(cohorts["coordinate_matched"].items()):
            vals = [r["auc"] for r in result["folds"] if r["auc"] is not None]
            ax.scatter(x + np.linspace(-0.09, 0.09, len(vals)), vals, s=23,
                       color="#63788b", label="Individual fold" if x == 0 else None)
            ax.scatter(x, result["pooled"], marker="D", s=49, color="#ac3b18",
                       label="Pooled" if x == 0 else None)
        ax.axhline(0.5, color="gray", linestyle=":", linewidth=1)
        ax.set(title=arm, xticks=range(3), xticklabels=["Position", "3D", "Block"],
               xlim=(-0.45, 2.45), ylim=(-0.03, 1.03))
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Position-only AUROC")
    fig.legend(*axes[0].get_legend_handles_labels(), loc="lower center", ncol=2, frameon=False)
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    fig.savefig(ROOT / "manuscript/figures/position_control_folds.png", dpi=300)
    plt.close(fig)


def main():
    df = pd.read_csv(SOURCE, sep="\t", usecols=[
        "gene", "drug", "pos", "target_binary", "has_3d_coord", *[s[1] for s in SPLIT_SCHEMES]])
    results = {}
    intervals = {}
    for arm, gene, drug in TARGETS:
        cohort = df.loc[(df.gene == gene) & (df.drug == drug)].copy()
        cohort["norm_pos"] = cohort.pos / PROTEIN_LENGTHS[gene]
        results[arm] = {}
        for key, subset in [("all", cohort), ("coordinate_matched", cohort[cohort.has_3d_coord == 1])]:
            frame = subset.reset_index(drop=True)
            summaries, scores = {}, {}
            for name, column in SPLIT_SCHEMES:
                if key == "all" and column == "spatial_cluster_id":
                    continue
                summaries[name], scores[column] = evaluate(frame, column)
            results[arm][key] = summaries
            if key == "coordinate_matched":
                intervals[arm] = paired_position_interval(frame, scores["cv_fold_id"], scores["seq_block_id"])
    output = dict(source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                  sklearn=sklearn.__version__, targets=results, matched_position_intervals=intervals)
    (ROOT / "logs/position_control_folds.json").write_text(
        json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    figure(results)


if __name__ == "__main__":
    main()
