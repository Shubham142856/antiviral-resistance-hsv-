"""
run_paired_family_contrast.py
=============================
Paired feature-family contrast under position-clustered bootstrap.

Why this exists. Two of the three proposed log-mining findings cannot be settled from
the frozen degradation log:

  * Comparing whether the conservation family's degradation is "significant" while the
    language-model family's is not is the difference-of-significance fallacy
    (Gelman & Stern, 2006) - the same fallacy this manuscript already names when
    discussing the radicalness dissociation. Two families can differ in significance
    purely because their effect sizes differ. The admissible statistic is the PAIRED
    difference of the two degradations, evaluated inside one shared resample.

  * Ranking families by smallest degradation delta is confounded: the unfitted BLOSUM62
    control has delta 0.0000 by construction. Absolute out-of-block AUROC is the
    unconfounded quantity, and its comparison across families must also be paired.

This script therefore computes, within each position-clustered bootstrap resample and
with the SAME resampled positions applied to every family:

  1. block AUROC difference  (PLM  -  conservation)             -> absolute performance
  2. paired delta-of-deltas  (delta_conservation - delta_PLM)    -> degradation contrast

Degradation delta is defined exactly as in the existing benchmark:
    delta = AUROC(position-grouped) - AUROC(contiguous-block)

Models are fitted once on the real data to obtain out-of-fold predictions; the bootstrap
then resamples those stored predictions, so nothing is refitted inside the loop. Runtime
is a few minutes on CPU and no GPU is required.

Output: logs/paired_family_contrast_results.json
"""

import json
import os
import warnings

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

DATA = os.path.join(ROOT, "data", "processed", "primary_3d_spatial_folds.tsv")
OUT = os.path.join(ROOT, "logs", "paired_family_contrast_results.json")

# Reference lengths used only to form norm_pos, matching the existing benchmark.
GENE_LEN = {"UL23": 376, "UL30": 1235}

FAMILIES = {
    "Positional Identity": ["norm_pos"],
    "Evolutionary Conservation": ["score_blosum62", "score_neg_entropy",
                                  "score_weighted_cons", "score_pssm"],
    "Protein Language Models": ["score_esm2_650m", "score_protbert", "score_esm2_35m"],
    "Multi-Modal Combined": ["norm_pos", "score_blosum62", "score_neg_entropy",
                              "score_weighted_cons", "score_pssm",
                              "score_esm2_650m", "score_protbert", "score_esm2_35m"],
}
SPLITS = {"position_grouped": "cv_fold_id", "contiguous_block": "seq_block_id"}
N_BOOT = 2000
SEED = 42


def oof_predictions(X, y, groups, model_kind):
    """Out-of-fold predicted probabilities, one pass over the pre-assigned folds."""
    pred = np.full(len(y), np.nan)
    for f in np.unique(groups):
        te = groups == f
        tr = ~te
        if tr.sum() == 0 or te.sum() == 0 or len(np.unique(y[tr])) < 2:
            continue
        Xtr, Xte = X[tr], X[te]
        if model_kind == "logistic_regression":
            sc = StandardScaler().fit(Xtr)
            Xtr, Xte = sc.transform(Xtr), sc.transform(Xte)
            clf = LogisticRegression(max_iter=1000, random_state=SEED)
        else:
            clf = RandomForestClassifier(n_estimators=100, max_depth=5,
                                         random_state=SEED, n_jobs=-1)
        clf.fit(Xtr, y[tr])
        pred[te] = clf.predict_proba(Xte)[:, 1]
    return pred


def main():
    df = pd.read_csv(DATA, sep="\t")
    df["norm_pos"] = [p / float(GENE_LEN[g]) for p, g in zip(df.pos, df.gene)]
    print("=" * 84)
    print("PAIRED FEATURE-FAMILY CONTRAST  (B=%d position-clustered, seed %d)"
          % (N_BOOT, SEED))
    print("=" * 84)

    out = {
        "purpose": ("paired contrasts that the frozen degradation log cannot supply: "
                    "absolute out-of-block PLM-minus-conservation difference, and the "
                    "paired delta-of-deltas between the two families' degradations"),
        "rationale": ("comparing one family's p-value against another's is the "
                      "difference-of-significance fallacy; and ranking by smallest "
                      "degradation is confounded because the unfitted control has "
                      "delta 0.0000 by construction"),
        "degradation_definition": "AUROC(position_grouped) - AUROC(contiguous_block)",
        "n_boot": N_BOOT, "seed": SEED, "targets": {},
    }

    for (gene, drug), sub0 in df.groupby(["gene", "drug"]):
        key = "%s-%s" % (gene, drug)
        cols = sorted({c for v in FAMILIES.values() for c in v})
        sub = sub0.dropna(subset=cols + ["target_binary"]).reset_index(drop=True)
        if len(sub) < 40:
            continue
        y = sub.target_binary.values.astype(int)
        print("")
        print("== %s   N=%d  R=%d" % (key, len(sub), int(y.sum())))

        rec = {"n_rows": int(len(sub)), "n_resistant": int(y.sum()), "models": {}}

        for model_kind, mlab in [("random_forest", "RF"),
                                 ("logistic_regression", "LR")]:
            # Stored OOF predictions per family per split; nothing refits in bootstrap.
            preds = {}
            for fam, feats in FAMILIES.items():
                X = sub[feats].values.astype(float)
                for sname, scol in SPLITS.items():
                    preds[(fam, sname)] = oof_predictions(
                        X, y, sub[scol].values, model_kind)

            point = {}
            for fam in FAMILIES:
                pg = roc_auc_score(y, preds[(fam, "position_grouped")])
                bk = roc_auc_score(y, preds[(fam, "contiguous_block")])
                point[fam] = {"position_grouped_auc": round(float(pg), 4),
                              "contiguous_block_auc": round(float(bk), 4),
                              "degradation_delta": round(float(pg - bk), 4)}

            # One shared sequence of resampled position clusters for every family.
            rng = np.random.RandomState(SEED)
            positions = sub.pos.unique()
            idx_by_pos = {p: np.where(sub.pos.values == p)[0] for p in positions}

            b_block_diff, b_delta_of_deltas = [], []
            for _ in range(N_BOOT):
                samp = rng.choice(positions, size=len(positions), replace=True)
                rows = np.concatenate([idx_by_pos[p] for p in samp])
                yb = y[rows]
                if len(np.unique(yb)) < 2:
                    continue
                try:
                    cb = roc_auc_score(yb, preds[("Evolutionary Conservation",
                                                  "contiguous_block")][rows])
                    pb = roc_auc_score(yb, preds[("Protein Language Models",
                                                  "contiguous_block")][rows])
                    cp = roc_auc_score(yb, preds[("Evolutionary Conservation",
                                                  "position_grouped")][rows])
                    pp = roc_auc_score(yb, preds[("Protein Language Models",
                                                  "position_grouped")][rows])
                except ValueError:
                    continue
                b_block_diff.append(pb - cb)
                b_delta_of_deltas.append((cp - cb) - (pp - pb))

            def summarize(v, label):
                a = np.asarray(v)
                lo, hi = np.percentile(a, [2.5, 97.5])
                p_le0 = float(np.mean(a <= 0))
                floor = 1.0 / N_BOOT
                excl = bool(lo > 0 or hi < 0)
                print("   %-6s %-34s %+.4f [%+.4f, %+.4f]  excludes 0: %-5s  p%s"
                      % (mlab, label, a.mean(), lo, hi, excl,
                         (" = %.4f" % p_le0) if p_le0 >= floor
                         else (" < %.4f" % floor)))
                return {"mean": round(float(a.mean()), 4),
                        "ci_95": [round(float(lo), 4), round(float(hi), 4)],
                        "excludes_zero": excl,
                        "p_one_sided": (round(p_le0, 4) if p_le0 >= floor
                                        else "< %.4f" % floor),
                        "n_valid_resamples": int(len(a))}

            res = {"point_estimates": point,
                   "paired_block_auroc_plm_minus_conservation":
                       summarize(b_block_diff, "block AUROC: PLM - conservation"),
                   "paired_delta_of_deltas_conservation_minus_plm":
                       summarize(b_delta_of_deltas,
                                 "degradation: conserv - PLM")}
            rec["models"][model_kind] = res
        out["targets"][key] = rec

    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)
    print("")
    print("Wrote %s" % os.path.relpath(OUT, ROOT))
    print("")
    print("Read the two paired rows per model. If 'block AUROC: PLM - conservation'")
    print("excludes zero, PLM genuinely outperforms conservation out of block on that")
    print("arm. If 'degradation: conserv - PLM' excludes zero, the two families really")
    print("do differ in how much they degrade; if it spans zero, they do not, and the")
    print("3-of-6 significance pattern was an artefact of comparing p-values.")


if __name__ == "__main__":
    main()
