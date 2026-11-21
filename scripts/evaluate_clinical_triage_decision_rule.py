"""Retrospective cohort accounting for the historical triage analysis."""
import os
import json
import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

def compute_binary_metrics(y_true, y_pred):
    y_t = np.asarray(y_true, dtype=int)
    y_p = np.asarray(y_pred, dtype=int)

    tp = int(np.sum((y_t == 1) & (y_p == 1)))
    fp = int(np.sum((y_t == 0) & (y_p == 1)))
    tn = int(np.sum((y_t == 0) & (y_p == 0)))
    fn = int(np.sum((y_t == 1) & (y_p == 0)))

    sens = float(tp / (tp + fn)) if (tp + fn) > 0 else 0.0
    spec = float(tn / (tn + fp)) if (tn + fp) > 0 else 0.0
    ppv = float(tp / (tp + fp)) if (tp + fp) > 0 else 0.0
    npv = float(tn / (tn + fn)) if (tn + fn) > 0 else 0.0
    bacc = float((sens + spec) / 2.0)
    acc = float((tp + tn) / len(y_t)) if len(y_t) > 0 else 0.0

    return {
        "N": len(y_t),
        "TP": tp,
        "FP": fp,
        "TN": tn,
        "FN": fn,
        "sensitivity": sens,
        "specificity": spec,
        "ppv": ppv,
        "npv": npv,
        "balanced_accuracy": bacc,
        "accuracy": acc
    }


def main():
    print("=" * 80)
    print("EVALUATING HIERARCHICAL CLINICAL DECISION TRIAGE RULE")
    print("=" * 80)

    # 1. Load active site definitions
    with open(os.path.join(ROOT, "logs", "active_site_dark_space_surveillance.json")) as f:
        dark_data = json.load(f)

    ul30_as_pos = set([r["pos"] for r in dark_data["ul30_catalytic_shell"]])
    ul23_as_pos = set([r["pos"] for r in dark_data["ul23_catalytic_shell"]])

    # 2. Load Missense Features
    feat_path = os.path.join(ROOT, "data", "processed", "where_what_features.tsv")
    df_feat = pd.read_csv(feat_path, sep="\t")

    df_feat = df_feat.loc[df_feat["virus"] == "HSV-1"].copy()

    # Cohort A: UL23-Aciclovir (missense)
    df_u23_ms = df_feat[(df_feat["gene"] == "UL23") & (df_feat["drug"].str.lower() == "aciclovir")].copy()
    y_true_u23_ms = df_u23_ms["target_binary"].values
    y_pred_u23_ms = []
    for _, r in df_u23_ms.iterrows():
        # Tier 2: Pocket proximity (<= 5.0 A)
        if r["dist_to_site_angstrom"] <= 5.0 or (r["pos"] in ul23_as_pos):
            y_pred_u23_ms.append(1)
        # Tier 3: Radical missense (severity_blosum >= 2)
        elif r["severity_blosum"] >= 2:
            y_pred_u23_ms.append(1)
        else:
            y_pred_u23_ms.append(0)
    metrics_u23_ms = compute_binary_metrics(y_true_u23_ms, y_pred_u23_ms)

    # Cohort B: UL23-Aciclovir (full clinical cohort including truncations)
    df_unified = pd.read_csv(os.path.join(ROOT, "data", "processed", "unified_modeling_dataset.tsv"), sep="\t")
    ul23_eval = df_unified.loc[
        df_unified["virus"].eq("HSV-1") & df_unified["gene"].eq("UL23")
        & df_unified["drug"].str.casefold().isin(["aciclovir", "acv"])
        & df_unified["target_binary"].notna()
    ].copy()
    keys = ["virus", "gene", "aa_change"]
    missense = ul23_eval.loc[ul23_eval.mutation_category.eq("single_substitution")]
    recorded = df_u23_ms.assign(recorded_prediction=y_pred_u23_ms)
    paired = missense.merge(recorded[keys + ["target_binary", "recorded_prediction"]],
                            on=keys, how="outer", validate="one_to_one", indicator=True)
    if not paired["_merge"].eq("both").all() or not paired.target_binary_x.eq(paired.target_binary_y).all():
        raise ValueError("Missense cohorts differ in membership or phenotype")
    other = ul23_eval.loc[ul23_eval.mutation_category.isin(["nonsense_stop", "frameshift", "indel"])]
    metrics_u23_full = compute_binary_metrics(
        np.r_[paired.target_binary_x, other.target_binary],
        np.r_[paired.recorded_prediction, np.ones(len(other), dtype=int)])
    excluded = ul23_eval.loc[~ul23_eval.index.isin(missense.index.union(other.index))]

    # Cohort C: UL30-Aciclovir
    df_u30_acv = df_feat[(df_feat["gene"] == "UL30") & (df_feat["drug"].str.lower() == "aciclovir")].copy()
    y_true_u30_acv = df_u30_acv["target_binary"].values
    y_pred_u30_acv = []
    for _, r in df_u30_acv.iterrows():
        if r["dist_to_site_angstrom"] <= 5.0 or (r["pos"] in ul30_as_pos):
            y_pred_u30_acv.append(1)
        elif r["severity_blosum"] >= 2:
            y_pred_u30_acv.append(1)
        else:
            y_pred_u30_acv.append(0)
    metrics_u30_acv = compute_binary_metrics(y_true_u30_acv, y_pred_u30_acv)

    # Cohort D: UL30-Foscarnet
    df_u30_fos = df_feat[(df_feat["gene"] == "UL30") & (df_feat["drug"].str.lower() == "foscarnet")].copy()
    y_true_u30_fos = df_u30_fos["target_binary"].values
    y_pred_u30_fos = []
    for _, r in df_u30_fos.iterrows():
        if r["dist_to_site_angstrom"] <= 5.0 or (r["pos"] in ul30_as_pos):
            y_pred_u30_fos.append(1)
        elif r["severity_blosum"] >= 2:
            y_pred_u30_fos.append(1)
        else:
            y_pred_u30_fos.append(0)
    metrics_u30_fos = compute_binary_metrics(y_true_u30_fos, y_pred_u30_fos)

    results = {
        "triage_architecture": {
            "tier_1": "Loss-of-Function Inactivation (UL23 frameshift / nonsense / indel) -> positive under historical rule",
            "tier_2": "3D Catalytic Pocket Proximity (dist <= 5.0 A in cryo-EM 8EXX / 8V1T / 1OF1) -> Resistant",
            "tier_3": "Radical Substitution Filter (severity_blosum >= 2 <=> BLOSUM62 <= -2) -> Resistant / High Risk"
        },
        "ul23_hsv1_labeled_records": len(ul23_eval),
        "excluded_categories": excluded.mutation_category.value_counts().to_dict(),
        "evaluations": {
            "UL23_Aciclovir_Evaluable_HSV1": metrics_u23_full,
            "UL23_Aciclovir_Missense_Cohort": metrics_u23_ms,
            "UL30_Aciclovir_Cohort": metrics_u30_acv,
            "UL30_Foscarnet_Cohort": metrics_u30_fos
        }
    }

    out_file = os.path.join(ROOT, "logs", "clinical_triage_evaluation.json")
    with open(out_file, "w") as f:
        json.dump(results, f, indent=2)
    print(f"Saved clinical triage evaluation to: {out_file}")

    print("\n" + "=" * 115)
    print(f"| {'Target Cohort':<35} | {'N':<4} | {'Sens':<7} | {'Spec':<7} | {'PPV':<7} | {'NPV':<7} | {'BAcc':<7} | {'TP/FP/TN/FN':<14} |")
    print("|" + "-" * 37 + "|" + "-" * 6 + "|" + "-" * 9 + "|" + "-" * 9 + "|" + "-" * 9 + "|" + "-" * 9 + "|" + "-" * 9 + "|" + "-" * 16 + "|")
    for name, m in results["evaluations"].items():
        cm = f"{m['TP']}/{m['FP']}/{m['TN']}/{m['FN']}"
        print(f"| {name:<35} | {m['N']:<4} | {m['sensitivity']*100:>5.1f}% | {m['specificity']*100:>5.1f}% | {m['ppv']*100:>5.1f}% | {m['npv']*100:>5.1f}% | {m['balanced_accuracy']*100:>5.1f}% | {cm:<14} |")
    print("=" * 115)

if __name__ == "__main__":
    main()
