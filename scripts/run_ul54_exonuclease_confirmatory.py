"""
run_ul54_exonuclease_confirmatory.py
====================================
Confirmatory execution of the pre-registered Exonuclease Active-Site Distance
and Target-Mechanism Double Dissociation hypothesis on HCMV UL54.

Anchored strictly to:
  docs/UL54_EXONUCLEASE_CONFIRMATORY_PREREGISTRATION.md (Version 1.0)

Evaluates:
  1. Primary Structural Features:
     - Dist_pol: Minimum distance to the 19-residue polymerase active-site shell
     - Dist_exo: Minimum distance to the 7-residue exonuclease catalytic pocket
       {301, 302, 303, 407, 408, 412, 413} (mapped from 8EXX MG 1303 coordination shell)
  2. Domain Boundary Sensitivity:
     - Proximity to linear domain windows: (300, 545), (300, 520), (280, 560)
  3. Pre-Registered Confirmatory Contrast:
     - Double Dissociation Crossover:
       Delta_Crossover = [S_0(GCV, exo) - S_0(GCV, pol)] - [S_0(FOS, exo) - S_0(FOS, pol)]
  4. Resampling:
     - Shared study-clustered bootstrap (B=2000, seed 42)

Output:
  - logs/step2_hcmv_exonuclease_confirmatory_results.json
"""

import os
import sys
import json
import re
import warnings
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser

warnings.filterwarnings('ignore')

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

EXO_CATALYTIC_RESIDUES = [301, 302, 303, 407, 408, 412, 413]
SENSITIVITY_WINDOWS = {
    "standard_exo_300_545": (300, 545),
    "tight_exo_300_520": (300, 520),
    "expanded_exo_280_560": (280, 560)
}
DRUGS = ["Ganciclovir", "Foscarnet", "Cidofovir"]
N_BOOT = 2000
RNG_SEED = 42

def main():
    print("=" * 80)
    print("CONFIRMATORY TEST: UL54 EXONUCLEASE ACTIVE-SITE DISTANCE & DOUBLE DISSOCIATION")
    print("Protocol: docs/UL54_EXONUCLEASE_CONFIRMATORY_PREREGISTRATION.md")
    print("=" * 80)

    # 1. Structural Distance Extraction on Swiss-Model P08546
    pdb_path = os.path.join(ROOT, "data", "structures", "P08546_swissmodel_8oj7.pdb")
    print(f"Loading structure: {os.path.relpath(pdb_path, ROOT)}")
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("UL54", pdb_path)
    model = list(struct.get_models())[0]

    # Load 19-residue polymerase active site
    with open(os.path.join(ROOT, "logs", "hcmv_ul54_active_site_mapping.json")) as f:
        pol_shell = sorted(set(json.load(f)["hcmv_active_site_positions"]))
    print(f"  Polymerase active-site shell (19 residues): {pol_shell}")
    print(f"  Exonuclease catalytic active site (7 residues): {EXO_CATALYTIC_RESIDUES}")

    ca_coords = {}
    pol_atoms = []
    exo_atoms = []

    for chain in model:
        for res in chain:
            r_id = res.id[1]
            if "CA" in res:
                ca_coords[r_id] = res["CA"].get_coord()
            if r_id in pol_shell:
                for a in res:
                    pol_atoms.append(a.get_coord())
            if r_id in EXO_CATALYTIC_RESIDUES:
                for a in res:
                    exo_atoms.append(a.get_coord())

    pol_atoms = np.array(pol_atoms)
    exo_atoms = np.array(exo_atoms)

    # Primary Distances
    d_pol = {p: float(np.min(np.linalg.norm(pol_atoms - c, axis=1))) for p, c in ca_coords.items()}
    d_exo = {p: float(np.min(np.linalg.norm(exo_atoms - c, axis=1))) for p, c in ca_coords.items()}

    # Sensitivity window distances (distance to nearest CA in window)
    window_dists = {}
    for win_name, (w_start, w_end) in SENSITIVITY_WINDOWS.items():
        w_cas = np.array([c for p, c in ca_coords.items() if w_start <= p <= w_end])
        window_dists[win_name] = {p: float(np.min(np.linalg.norm(w_cas - c, axis=1))) for p, c in ca_coords.items()}

    print(f"Extracted distances for {len(ca_coords)} residues.")
    print(f"  Median dist to pol shell: {np.median(list(d_pol.values())):.1f} A")
    print(f"  Median dist to exo site:  {np.median(list(d_exo.values())):.1f} A")

    # 2. Phenotyping & Cohort Extraction from HerpesDRG
    raw_path = os.path.join(ROOT, "data", "raw", "herpesdrg-db", "herpesdrg-db.tsv")
    raw_df = pd.read_csv(raw_path, sep="\t")
    hcmv_df = raw_df[(raw_df["virus"] == "HCMV") & (raw_df["gene"] == "UL54")].copy()

    def parse_mut(s):
        m = re.match(r"^([A-Z])(\d+)([A-Z])$", str(s).strip())
        return (m.group(1), int(m.group(2)), m.group(3)) if m else (None, None, None)

    p = [parse_mut(x) for x in hcmv_df["aa_change"]]
    hcmv_df["wt"] = [x[0] for x in p]
    hcmv_df["pos"] = [x[1] for x in p]
    hcmv_df["mut"] = [x[2] for x in p]
    hcmv_df = hcmv_df[hcmv_df["wt"].notna()].copy()
    hcmv_df["pos"] = hcmv_df["pos"].astype(int)

    def get_study_id(row):
        doi = str(row["ref_doi"]).strip()
        if doi and doi not in ["nan", "None", ""]:
            return f"DOI:{doi}"
        t = str(row["ref_title"]).strip()
        if t and t not in ["nan", "None", ""]:
            return f"TITLE:{t[:45]}"
        return "Unassigned"

    hcmv_df["study_id"] = [get_study_id(r) for _, r in hcmv_df.iterrows()]

    def is_res(val, thresh=2.0):
        s = str(val).strip()
        if s == "Resistant":
            return 1
        if s in ["Polymorphism", "Sensitive"]:
            return 0
        try:
            return 1 if float(s) >= thresh else 0
        except:
            return None

    # Prepare datasets per drug
    max_d_pol = max(d_pol.values())
    max_d_exo = max(d_exo.values())

    drug_dfs = {}
    for d in DRUGS:
        sub = hcmv_df.copy()
        sub["y"] = [is_res(x) for x in sub[d]]
        sub = sub[sub["y"].notna()].copy()
        sub["y"] = sub["y"].astype(int)
        sub = sub.drop_duplicates(subset=["aa_change", "study_id"]).copy().reset_index(drop=True)

        # Attach scores (higher score = closer = predicted resistant)
        sub["score_pol"] = [- d_pol.get(p, max_d_pol) for p in sub["pos"]]
        sub["score_exo"] = [- d_exo.get(p, max_d_exo) for p in sub["pos"]]

        for win_name, win_d in window_dists.items():
            max_w = max(win_d.values())
            sub[f"score_{win_name}"] = [- win_d.get(p, max_w) for p in sub["pos"]]

        drug_dfs[d] = sub
        print(f"Loaded {d}: {len(sub)} rows, {int(sub['y'].sum())} R / {int((sub['y'] == 0).sum())} S")

    # 3. Helper for Within-Study Concordance
    def compute_study_pair_stats(df_sub, score_col):
        st_groups = df_sub.groupby("study_id")
        stats = {}
        for st_id, grp in st_groups:
            r_sc = grp.loc[grp["y"] == 1, score_col].values
            s_sc = grp.loc[grp["y"] == 0, score_col].values
            n_r = len(r_sc)
            n_s = len(s_sc)
            if n_r > 0 and n_s > 0:
                n_p = n_r * n_s
                diffs = r_sc[:, None] - s_sc[None, :]
                concordant = float(np.sum(diffs > 0) + 0.5 * np.sum(diffs == 0))
                stats[st_id] = {
                    "n_res": n_r,
                    "n_sens": n_s,
                    "n_pairs": n_p,
                    "concordant": concordant,
                    "auroc": concordant / n_p
                }
        return stats

    def pair_weighted_from_stats(stats_dict):
        tot_p = sum(s["n_pairs"] for s in stats_dict.values())
        tot_c = sum(s["concordant"] for s in stats_dict.values())
        return (tot_c / tot_p) if tot_p > 0 else np.nan, tot_p

    def study_mean_from_stats(stats_dict):
        if not stats_dict:
            return np.nan
        return float(np.mean([s["auroc"] for s in stats_dict.values()]))

    # Calculate observed statistics
    observed = {}
    for d in DRUGS:
        df_d = drug_dfs[d]
        stats_pol = compute_study_pair_stats(df_d, "score_pol")
        stats_exo = compute_study_pair_stats(df_d, "score_exo")

        s0_pol, np_pol = pair_weighted_from_stats(stats_pol)
        s0_exo, np_exo = pair_weighted_from_stats(stats_exo)

        mean_pol = study_mean_from_stats(stats_pol)
        mean_exo = study_mean_from_stats(stats_exo)

        sens_results = {}
        for win_name in SENSITIVITY_WINDOWS:
            st_win = compute_study_pair_stats(df_d, f"score_{win_name}")
            s0_w, _ = pair_weighted_from_stats(st_win)
            sens_results[win_name] = round(float(s0_w), 4)

        observed[d] = {
            "s0_pol": float(s0_pol),
            "s0_exo": float(s0_exo),
            "study_mean_pol": float(mean_pol),
            "study_mean_exo": float(mean_exo),
            "delta_exo_minus_pol": float(s0_exo - s0_pol),
            "n_pairs": int(np_pol),
            "m_mixed_studies": len(stats_pol),
            "sensitivity_windows_s0": sens_results
        }

    # 4. Joint Study-Clustered Bootstrap Engine (B=2000, seed 42)
    print("\nRunning joint study-clustered bootstrap (B=2,000 resamples, seed 42)...")
    rng = np.random.RandomState(RNG_SEED)

    # Pre-extract per-study stats for fast bootstrapping
    pre_stats = {}
    for d in DRUGS:
        df_d = drug_dfs[d]
        pre_stats[d] = {
            "studies": list(df_d["study_id"].unique()),
            "pol": compute_study_pair_stats(df_d, "score_pol"),
            "exo": compute_study_pair_stats(df_d, "score_exo")
        }

    boot_records = {d: {"s0_pol": [], "s0_exo": [], "delta": []} for d in DRUGS}
    crossover_boot_gcv = []
    crossover_boot_cdv = []

    for b in range(N_BOOT):
        # We bootstrap studies per drug
        b_s0 = {}
        for d in DRUGS:
            st_list = pre_stats[d]["studies"]
            sample_st = rng.choice(st_list, size=len(st_list), replace=True)
            
            tot_p_pol, tot_c_pol = 0, 0.0
            tot_p_exo, tot_c_exo = 0, 0.0

            for st_id in sample_st:
                sp = pre_stats[d]["pol"].get(st_id)
                se = pre_stats[d]["exo"].get(st_id)
                if sp:
                    tot_p_pol += sp["n_pairs"]
                    tot_c_pol += sp["concordant"]
                if se:
                    tot_p_exo += se["n_pairs"]
                    tot_c_exo += se["concordant"]

            s0_p = (tot_c_pol / tot_p_pol) if tot_p_pol > 0 else np.nan
            s0_e = (tot_c_exo / tot_p_exo) if tot_p_exo > 0 else np.nan
            d_val = s0_e - s0_p if (not np.isnan(s0_e) and not np.isnan(s0_p)) else np.nan

            boot_records[d]["s0_pol"].append(s0_p)
            boot_records[d]["s0_exo"].append(s0_e)
            boot_records[d]["delta"].append(d_val)
            b_s0[d] = {"pol": s0_p, "exo": s0_e, "delta": d_val}

        # Crossover deltas
        # Delta_Crossover = [Delta(GCV)] - [Delta(FOS)]
        if not np.isnan(b_s0["Ganciclovir"]["delta"]) and not np.isnan(b_s0["Foscarnet"]["delta"]):
            crossover_boot_gcv.append(b_s0["Ganciclovir"]["delta"] - b_s0["Foscarnet"]["delta"])
        if not np.isnan(b_s0["Cidofovir"]["delta"]) and not np.isnan(b_s0["Foscarnet"]["delta"]):
            crossover_boot_cdv.append(b_s0["Cidofovir"]["delta"] - b_s0["Foscarnet"]["delta"])

    # Compute CIs
    final_drug_results = {}
    for d in DRUGS:
        pol_ci = [float(np.nanpercentile(boot_records[d]["s0_pol"], 2.5)),
                  float(np.nanpercentile(boot_records[d]["s0_pol"], 97.5))]
        exo_ci = [float(np.nanpercentile(boot_records[d]["s0_exo"], 2.5)),
                  float(np.nanpercentile(boot_records[d]["s0_exo"], 97.5))]
        del_ci = [float(np.nanpercentile(boot_records[d]["delta"], 2.5)),
                  float(np.nanpercentile(boot_records[d]["delta"], 97.5))]
        p_del = float(np.mean(np.array(boot_records[d]["delta"]) <= 0)) if observed[d]["delta_exo_minus_pol"] > 0 else float(np.mean(np.array(boot_records[d]["delta"]) >= 0))

        final_drug_results[d] = {
            "s0_polymerase": {
                "point": round(observed[d]["s0_pol"], 4),
                "ci_95": [round(pol_ci[0], 4), round(pol_ci[1], 4)],
                "study_mean": round(observed[d]["study_mean_pol"], 4)
            },
            "s0_exonuclease": {
                "point": round(observed[d]["s0_exo"], 4),
                "ci_95": [round(exo_ci[0], 4), round(exo_ci[1], 4)],
                "study_mean": round(observed[d]["study_mean_exo"], 4)
            },
            "delta_exo_minus_pol": {
                "observed": round(observed[d]["delta_exo_minus_pol"], 4),
                "ci_95": [round(del_ci[0], 4), round(del_ci[1], 4)],
                "p_value": round(p_del, 4)
            },
            "n_pairs": observed[d]["n_pairs"],
            "m_mixed_studies": observed[d]["m_mixed_studies"],
            "sensitivity_windows_s0": observed[d]["sensitivity_windows_s0"]
        }

    # Crossover results
    obs_xover_gcv = observed["Ganciclovir"]["delta_exo_minus_pol"] - observed["Foscarnet"]["delta_exo_minus_pol"]
    ci_xover_gcv = [float(np.nanpercentile(crossover_boot_gcv, 2.5)),
                    float(np.nanpercentile(crossover_boot_gcv, 97.5))]
    p_xover_gcv = float(np.mean(np.array(crossover_boot_gcv) <= 0))

    obs_xover_cdv = observed["Cidofovir"]["delta_exo_minus_pol"] - observed["Foscarnet"]["delta_exo_minus_pol"]
    ci_xover_cdv = [float(np.nanpercentile(crossover_boot_cdv, 2.5)),
                    float(np.nanpercentile(crossover_boot_cdv, 97.5))]
    p_xover_cdv = float(np.mean(np.array(crossover_boot_cdv) <= 0))

    # Evaluate Pre-Registered Rules
    # Rule 1: GCV Exo Recovery
    gcv_pt = final_drug_results["Ganciclovir"]["s0_exonuclease"]["point"]
    gcv_lo = final_drug_results["Ganciclovir"]["s0_exonuclease"]["ci_95"][0]
    if gcv_pt >= 0.550 and gcv_lo > 0.500:
        rule1_verdict = "CONFIRMED (Point >= 0.550 and CI excludes 0.50)"
    elif gcv_pt >= 0.550:
        rule1_verdict = "INCONCLUSIVE (Point >= 0.550 but CI spans 0.50)"
    else:
        rule1_verdict = "NOT SUPPORTED (Point < 0.550)"

    # Rule 2: CDV Exo Recovery
    cdv_pt = final_drug_results["Cidofovir"]["s0_exonuclease"]["point"]
    cdv_lo = final_drug_results["Cidofovir"]["s0_exonuclease"]["ci_95"][0]
    if cdv_pt >= 0.550 and cdv_lo > 0.500:
        rule2_verdict = "CONFIRMED (Point >= 0.550 and CI excludes 0.50)"
    elif cdv_pt >= 0.550:
        rule2_verdict = "INCONCLUSIVE (Point >= 0.550 but CI spans 0.50)"
    else:
        rule2_verdict = "NOT SUPPORTED (Point < 0.500)"

    # Rule 3: FOS Degradation & Crossover
    fos_exo_pt = final_drug_results["Foscarnet"]["s0_exonuclease"]["point"]
    fos_pol_pt = final_drug_results["Foscarnet"]["s0_polymerase"]["point"]
    if fos_exo_pt < fos_pol_pt and obs_xover_gcv > 0 and ci_xover_gcv[0] > 0:
        rule3_verdict = "CONFIRMED (Double Dissociation Crossover Statistically Significant)"
    elif fos_exo_pt < fos_pol_pt and obs_xover_gcv > 0:
        rule3_verdict = "INCONCLUSIVE (Crossover Direction Supported, CI Spans Zero)"
    else:
        rule3_verdict = "VIOLATED (No Exonuclease Degradation on FOS)"

    print("\n" + "=" * 80)
    print("EMPIRICAL CONFIRMATORY RESULTS:")
    print("=" * 80)
    print(f"1. Ganciclovir (Primary Replication Drug):")
    print(f"   - Polymerase Active-Site S0 : {final_drug_results['Ganciclovir']['s0_polymerase']['point']} {final_drug_results['Ganciclovir']['s0_polymerase']['ci_95']}")
    print(f"   - Exonuclease Active-Site S0: {final_drug_results['Ganciclovir']['s0_exonuclease']['point']} {final_drug_results['Ganciclovir']['s0_exonuclease']['ci_95']}")
    print(f"   - Delta (Exo - Pol)         : {final_drug_results['Ganciclovir']['delta_exo_minus_pol']['observed']:+} {final_drug_results['Ganciclovir']['delta_exo_minus_pol']['ci_95']} (p = {final_drug_results['Ganciclovir']['delta_exo_minus_pol']['p_value']:.4f})")
    print(f"   - Rule 1 Verdict            : {rule1_verdict}")

    print(f"\n2. Cidofovir (Secondary Nucleotide Analogue):")
    print(f"   - Polymerase Active-Site S0 : {final_drug_results['Cidofovir']['s0_polymerase']['point']} {final_drug_results['Cidofovir']['s0_polymerase']['ci_95']}")
    print(f"   - Exonuclease Active-Site S0: {final_drug_results['Cidofovir']['s0_exonuclease']['point']} {final_drug_results['Cidofovir']['s0_exonuclease']['ci_95']}")
    print(f"   - Delta (Exo - Pol)         : {final_drug_results['Cidofovir']['delta_exo_minus_pol']['observed']:+} {final_drug_results['Cidofovir']['delta_exo_minus_pol']['ci_95']} (p = {final_drug_results['Cidofovir']['delta_exo_minus_pol']['p_value']:.4f})")
    print(f"   - Rule 2 Verdict            : {rule2_verdict}")

    print(f"\n3. Foscarnet (Pyrophosphate Analogue Counter-Control):")
    print(f"   - Polymerase Active-Site S0 : {final_drug_results['Foscarnet']['s0_polymerase']['point']} {final_drug_results['Foscarnet']['s0_polymerase']['ci_95']}")
    print(f"   - Exonuclease Active-Site S0: {final_drug_results['Foscarnet']['s0_exonuclease']['point']} {final_drug_results['Foscarnet']['s0_exonuclease']['ci_95']}")
    print(f"   - Delta (Exo - Pol)         : {final_drug_results['Foscarnet']['delta_exo_minus_pol']['observed']:+} {final_drug_results['Foscarnet']['delta_exo_minus_pol']['ci_95']}")

    print(f"\n4. Formal Double Dissociation Crossover Tests:")
    print(f"   - GCV vs FOS Crossover Delta: {obs_xover_gcv:+.4f} [{ci_xover_gcv[0]:.4f}, {ci_xover_gcv[1]:.4f}] (p = {p_xover_gcv:.4f})")
    print(f"   - CDV vs FOS Crossover Delta: {obs_xover_cdv:+.4f} [{ci_xover_cdv[0]:.4f}, {ci_xover_cdv[1]:.4f}] (p = {p_xover_cdv:.4f})")
    print(f"   - Rule 3 Verdict            : {rule3_verdict}")

    print(f"\n5. Sensitivity Windows S0 (Exo Boundary Robustness):")
    for d in DRUGS:
        print(f"   - {d}: {final_drug_results[d]['sensitivity_windows_s0']}")

    out = {
        "protocol": {
            "document": "docs/UL54_EXONUCLEASE_CONFIRMATORY_PREREGISTRATION.md",
            "version": "1.0",
            "status": "CONFIRMATORY EXECUTION",
            "n_boot": N_BOOT,
            "seed": RNG_SEED
        },
        "definitions": {
            "polymerase_active_site_shell_19": pol_shell,
            "exonuclease_catalytic_site_7": EXO_CATALYTIC_RESIDUES,
            "sensitivity_windows": SENSITIVITY_WINDOWS
        },
        "targets": final_drug_results,
        "crossover_contrasts": {
            "gcv_vs_fos": {
                "observed_delta": round(obs_xover_gcv, 4),
                "ci_95": [round(ci_xover_gcv[0], 4), round(ci_xover_gcv[1], 4)],
                "p_value": round(p_xover_gcv, 4)
            },
            "cdv_vs_fos": {
                "observed_delta": round(obs_xover_cdv, 4),
                "ci_95": [round(ci_xover_cdv[0], 4), round(ci_xover_cdv[1], 4)],
                "p_value": round(p_xover_cdv, 4)
            }
        },
        "verdicts": {
            "rule1_gcv_exonuclease_recovery": rule1_verdict,
            "rule2_cdv_exonuclease_recovery": rule2_verdict,
            "rule3_double_dissociation_crossover": rule3_verdict
        }
    }

    out_path = os.path.join(ROOT, "logs", "step2_hcmv_exonuclease_confirmatory_results.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nWrote results to: {os.path.relpath(out_path, ROOT)}")

if __name__ == "__main__":
    main()

