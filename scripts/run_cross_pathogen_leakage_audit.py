"""
run_cross_pathogen_leakage_audit.py
===================================
Diagnostic audit quantifying Sequence-Neighbour Gap spatial leakage across
compact viral enzymes in multiple human pathogens:
  1. HSV-1 UL23 (Thymidine Kinase, Herpesviridae)
  2. HSV-1 UL30 (DNA Polymerase, Herpesviridae)
  3. HCMV UL54 (DNA Polymerase, Herpesviridae)
  4. HIV-1 Protease (PR, Retroviridae)
  5. HIV-1 Reverse Transcriptase (RT, Retroviridae)
  6. SARS-CoV-2 Main Protease (Mpro, Coronaviridae)

Evaluates:
  - Modulo-5 (position-grouped / modulo-K) cross-validation
  - Contiguous-Block (5 blocks, block size ceil(N/5)) cross-validation

Outputs:
  - logs/cross_pathogen_leakage_audit.json
"""

import os
import re
import json
import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

def load_targets():
    # 1. HSV-1 targets from processed dataset
    hsv_path = os.path.join(ROOT, "data", "processed", "where_what_features.tsv")
    df_hsv = pd.read_csv(hsv_path, sep="\t")
    ul23_pos = sorted(df_hsv[df_hsv["gene"] == "UL23"]["pos"].dropna().astype(int).unique())
    ul30_pos = sorted(df_hsv[df_hsv["gene"] == "UL30"]["pos"].dropna().astype(int).unique())

    # 2. HCMV UL54 target from HerpesDRG
    herpesdrg_path = os.path.join(ROOT, "data", "raw", "herpesdrg-db", "herpesdrg-db.tsv")
    df_raw = pd.read_csv(herpesdrg_path, sep="\t")
    df_hcmv = df_raw[(df_raw["virus"] == "HCMV") & (df_raw["gene"] == "UL54")].copy()
    def parse_pos(s):
        m = re.match(r"^([A-Z])(\d+)([A-Z])$", str(s).strip())
        return int(m.group(2)) if m else None
    df_hcmv["pos"] = [parse_pos(x) for x in df_hcmv["aa_change"]]
    hcmv_pos = sorted(df_hcmv["pos"].dropna().astype(int).unique())

    # 3. HIV-1 Protease (PR) major DRMs (Stanford HIVdb / IAS-USA 2022)
    hiv_pr_pos = [30, 32, 46, 47, 48, 50, 54, 76, 82, 84, 88, 90]

    # 4. HIV-1 Reverse Transcriptase (RT) major NRTI/NNRTI DRMs (Stanford HIVdb / IAS-USA 2022)
    hiv_rt_pos = [41, 65, 67, 70, 74, 100, 103, 106, 115, 151, 181, 184, 188, 190, 210, 215, 219, 227, 230, 236]

    # 5. SARS-CoV-2 Main Protease (Mpro) nirmatrelvir/Paxlovid DRMs (Moghadasi et al. 2023, Iketani et al. 2023)
    sars2_mpro_pos = [50, 140, 143, 144, 165, 166, 167, 172, 189, 192]

    return {
        "HSV-1 UL23 (Thymidine Kinase)": {
            "family": "Herpesviridae",
            "enzyme_class": "Kinase",
            "len_aa": 376,
            "positions": ul23_pos,
            "source": "HerpesDRG curated clinical cohort"
        },
        "HSV-1 UL30 (DNA Polymerase)": {
            "family": "Herpesviridae",
            "enzyme_class": "Polymerase",
            "len_aa": 1235,
            "positions": ul30_pos,
            "source": "HerpesDRG curated clinical cohort"
        },
        "HCMV UL54 (DNA Polymerase)": {
            "family": "Herpesviridae",
            "enzyme_class": "Polymerase",
            "len_aa": 1242,
            "positions": hcmv_pos,
            "source": "HerpesDRG curated clinical cohort"
        },
        "HIV-1 Protease (PR)": {
            "family": "Retroviridae",
            "enzyme_class": "Protease",
            "len_aa": 99,
            "positions": hiv_pr_pos,
            "source": "Stanford HIVdb / IAS-USA 2022 major DRMs"
        },
        "HIV-1 Reverse Transcriptase (RT)": {
            "family": "Retroviridae",
            "enzyme_class": "Polymerase",
            "len_aa": 440,
            "positions": hiv_rt_pos,
            "source": "Stanford HIVdb / IAS-USA 2022 major DRMs"
        },
        "SARS-CoV-2 Main Protease (Mpro)": {
            "family": "Coronaviridae",
            "enzyme_class": "Protease",
            "len_aa": 306,
            "positions": sars2_mpro_pos,
            "source": "Moghadasi et al. 2023 / Iketani et al. 2023"
        }
    }

def compute_split_gaps(positions, k_folds=5):
    pos_arr = np.array(sorted(positions))
    N = len(pos_arr)

    # Modulo-K
    mod_gaps = []
    for k in range(k_folds):
        test_mask = (pos_arr % k_folds == k)
        test_pos = pos_arr[test_mask]
        train_pos = pos_arr[~test_mask]
        if len(test_pos) == 0 or len(train_pos) == 0:
            continue
        for tp in test_pos:
            gaps = np.abs(train_pos - tp)
            mod_gaps.append(int(np.min(gaps)))

    # Contiguous-Block (5 contiguous partitions)
    block_size = int(np.ceil(N / k_folds))
    block_gaps = []
    for k in range(k_folds):
        start = k * block_size
        end = min((k + 1) * block_size, N)
        test_pos = pos_arr[start:end]
        train_pos = np.concatenate([pos_arr[:start], pos_arr[end:]])
        if len(test_pos) == 0 or len(train_pos) == 0:
            continue
        for tp in test_pos:
            gaps = np.abs(train_pos - tp)
            block_gaps.append(int(np.min(gaps)))

    mod_arr = np.array(mod_gaps)
    blk_arr = np.array(block_gaps)

    # Statistical test: Wilcoxon signed-rank test on paired gaps per test residue
    # Since test positions are identical in both evaluations, we compare paired gaps:
    w_stat, w_pval = stats.wilcoxon(blk_arr, mod_arr, alternative="greater")

    return {
        "N_positions": N,
        "modulo": {
            "median": float(np.median(mod_arr)),
            "mean": float(np.mean(mod_arr)),
            "std": float(np.std(mod_arr)),
            "p25": float(np.percentile(mod_arr, 25)),
            "p10": float(np.percentile(mod_arr, 10)),
            "pct_le_1aa": float(np.mean(mod_arr <= 1) * 100),
            "pct_le_2aa": float(np.mean(mod_arr <= 2) * 100),
            "pct_le_3aa": float(np.mean(mod_arr <= 3) * 100),
            "pct_le_5aa": float(np.mean(mod_arr <= 5) * 100),
            "gaps": [int(x) for x in mod_arr]
        },
        "contiguous_block": {
            "median": float(np.median(blk_arr)),
            "mean": float(np.mean(blk_arr)),
            "std": float(np.std(blk_arr)),
            "p25": float(np.percentile(blk_arr, 25)),
            "p10": float(np.percentile(blk_arr, 10)),
            "pct_le_1aa": float(np.mean(blk_arr <= 1) * 100),
            "pct_le_2aa": float(np.mean(blk_arr <= 2) * 100),
            "pct_le_3aa": float(np.mean(blk_arr <= 3) * 100),
            "pct_le_5aa": float(np.mean(blk_arr <= 5) * 100),
            "gaps": [int(x) for x in blk_arr]
        },
        "wilcoxon_stat": float(w_stat),
        "wilcoxon_pval": float(w_pval),
        "leakage_mitigation_ratio": float(np.median(blk_arr) / max(np.median(mod_arr), 1e-6))
    }

def main():
    print("=" * 80)
    print("RUNNING CROSS-PATHOGEN SEQUENCE-NEIGHBOUR GAP AUDIT")
    print("=" * 80)

    targets = load_targets()
    results = {}

    for name, info in targets.items():
        print(f"Auditing: {name} (N = {len(info['positions'])} positions, {info['source']})...")
        audit_res = compute_split_gaps(info["positions"], k_folds=5)
        results[name] = {
            "family": info["family"],
            "enzyme_class": info["enzyme_class"],
            "len_aa": info["len_aa"],
            "source": info["source"],
            "metrics": audit_res
        }

    # Print summary table
    print("\n" + "=" * 120)
    print(f"| {'Target':<32} | {'Family':<14} | {'N':<4} | {'Mod Med':<8} | {'Mod <=2aa':<10} | {'Mod <=5aa':<10} | {'Blk Med':<8} | {'Blk <=2aa':<10} | {'Blk <=5aa':<10} | {'P-value':<10} |")
    print("|" + "-" * 34 + "|" + "-" * 16 + "|" + "-" * 6 + "|" + "-" * 10 + "|" + "-" * 12 + "|" + "-" * 12 + "|" + "-" * 10 + "|" + "-" * 12 + "|" + "-" * 12 + "|" + "-" * 12 + "|")
    for name, r in results.items():
        m = r["metrics"]
        print(f"| {name:<32} | {r['family']:<14} | {m['N_positions']:<4} | {m['modulo']['median']:>5.1f} aa | {m['modulo']['pct_le_2aa']:>8.1f}% | {m['modulo']['pct_le_5aa']:>8.1f}% | {m['contiguous_block']['median']:>5.1f} aa | {m['contiguous_block']['pct_le_2aa']:>8.1f}% | {m['contiguous_block']['pct_le_5aa']:>8.1f}% | {m['wilcoxon_pval']:>10.2e} |")
    print("=" * 120)

    out_file = os.path.join(ROOT, "logs", "cross_pathogen_leakage_audit.json")
    with open(out_file, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved cross-pathogen audit to: {out_file}")

if __name__ == "__main__":
    main()
