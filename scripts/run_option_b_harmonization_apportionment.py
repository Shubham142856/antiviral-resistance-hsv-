"""
run_option_b_harmonization_apportionment.py
===========================================
Formal Option B Harmonization and Artifact Apportionment Experiment.

Objectives:
  1. Harmonize Option B (PDB 8EXX coordinates) to the identical 19-residue
     catalytic shell used in Option A (B_shell19), filtering out distal MG 1303.
  2. Factorially isolate and apportion the two confounded causes of the old Option B divergence:
     - Factor 1: Reference-set specification (Distal MG 1303 + triad vs. Clean 19-residue shell)
     - Factor 2: Missingness imputation policy (Complete cases N=133 vs. -max_dist imputation N=180)
  3. Log full reproducible results to:
     logs/option_b_harmonization_and_apportionment_results.json
"""

import os
import sys
import json
import re
import warnings
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, MMCIFParser
from Bio.Align import PairwiseAligner, substitution_matrices

warnings.filterwarnings('ignore')

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from scripts.map_active_site_to_hcmv import get_gene_seq

DRUGS = ["Ganciclovir", "Foscarnet", "Cidofovir"]
HSV1_SHELL_19 = [717, 718, 719, 720, 721, 722, 723, 785, 789, 811,
                 812, 815, 816, 818, 819, 887, 888, 925, 927]
HSV1_TRIAD = [717, 888, 889]

def ca_of(model, chain_id=None):
    out = {}
    for chain in model:
        if chain_id is not None and chain.id != chain_id:
            continue
        for res in chain:
            if "CA" in res:
                out[res.id[1]] = res["CA"].get_coord()
    return out

def min_dist(ref_atoms, ca_coords):
    ref = np.asarray(ref_atoms)
    return {p: float(np.min(np.linalg.norm(ref - c, axis=1))) for p, c in ca_coords.items()}

def collect(model, chain_id, residues=(), ligands=()):
    atoms = []
    for res in model.get_residues():
        if res.get_resname() in ligands:
            atoms += [a.get_coord() for a in res]
        if res.get_parent().id == chain_id and res.id[0] == " " and res.id[1] in residues:
            atoms += [a.get_coord() for a in res]
    return atoms

def main():
    print("=" * 80)
    print("OPTION B HARMONIZATION & ARTIFACT APPORTIONMENT BENCHMARK")
    print("=" * 80)

    # 1. Alignment HSV-1 UL30 <-> HCMV UL54
    seq_hsv = get_gene_seq("NC_001806.2", "UL30")
    seq_hcmv = get_gene_seq("NC_006273.2", "UL54")
    al = PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score, al.mode = -10.0, -0.5, "global"
    aln = al.align(seq_hsv, seq_hcmv)[0]
    hi = ci = 0
    hcmv_to_hsv, hsv_to_hcmv = {}, {}
    for r, q in zip(aln[0], aln[1]):
        if r != "-": hi += 1
        if q != "-": ci += 1
        if r != "-" and q != "-":
            hcmv_to_hsv[ci] = hi
            hsv_to_hcmv[hi] = ci

    # 2. Coordinates
    sm = list(PDBParser(QUIET=True).get_structure("UL54", os.path.join(ROOT, "data", "structures", "P08546_swissmodel_8oj7.pdb")).get_models())[0]
    ca_sm = ca_of(sm)

    cif = MMCIFParser(QUIET=True)
    exx = list(cif.get_structure("8EXX", os.path.join(ROOT, "data", "structures", "8exx.cif")).get_models())[0]
    ca_exx = ca_of(exx, "A")

    with open(os.path.join(ROOT, "logs", "hcmv_ul54_active_site_mapping.json")) as f:
        hcmv_shell_19 = sorted(set(json.load(f)["hcmv_active_site_positions"]))

    # Option A (Swiss-Model + shell19)
    d_A = min_dist(collect(sm, list(sm)[0].id, residues=hcmv_shell_19), ca_sm)

    # Option B Clean (8EXX + shell19 mapped)
    d_sh_exx = min_dist(collect(exx, "A", residues=HSV1_SHELL_19), ca_exx)
    d_B_shell19 = {p: d_sh_exx[h] for p, h in hcmv_to_hsv.items() if h in d_sh_exx}

    # Option B Old (8EXX + PPF + all MGs including distal MG 1303 + triad {717,888,889})
    d_old_exx = min_dist(collect(exx, "A", residues=HSV1_TRIAD, ligands=("PPF", "MG")), ca_exx)
    d_B_old = {p: d_old_exx[h] for p, h in hcmv_to_hsv.items() if h in d_old_exx}

    # 3. Load HerpesDRG HCMV UL54
    raw_df = pd.read_csv(os.path.join(ROOT, "data", "raw", "herpesdrg-db", "herpesdrg-db.tsv"), sep="\t")
    hcmv_df = raw_df[(raw_df["virus"] == "HCMV") & (raw_df["gene"] == "UL54")].copy()

    def parse_mut(s):
        m = re.match(r"^([A-Z])(\d+)([A-Z])$", str(s).strip())
        return (m.group(1), int(m.group(2)), m.group(3)) if m else (None, None, None)

    p = [parse_mut(x) for x in hcmv_df["aa_change"]]
    hcmv_df["pos"] = [x[1] for x in p]
    hcmv_df = hcmv_df[hcmv_df["pos"].notna()].copy()
    hcmv_df["pos"] = hcmv_df["pos"].astype(int)

    def get_study_id(row):
        doi = str(row["ref_doi"]).strip()
        if doi and doi not in ["nan", "None", ""]: return f"DOI:{doi}"
        t = str(row["ref_title"]).strip()
        if t and t not in ["nan", "None", ""]: return f"TITLE:{t[:45]}"
        return "Unassigned"

    hcmv_df["study_id"] = [get_study_id(r) for _, r in hcmv_df.iterrows()]

    def is_res(val, thresh=2.0):
        s = str(val).strip()
        if s == "Resistant": return 1
        if s in ["Polymorphism", "Sensitive"]: return 0
        try: return 1 if float(s) >= thresh else 0
        except: return None

    def calc_study_stats(df_sub, score_col):
        st_groups = df_sub.groupby("study_id")
        tot_c, tot_p = 0.0, 0
        per_study = []
        for st_id, grp in st_groups:
            r_sc = grp.loc[grp["y"] == 1, score_col].values
            s_sc = grp.loc[grp["y"] == 0, score_col].values
            if len(r_sc) > 0 and len(s_sc) > 0:
                npairs = len(r_sc) * len(s_sc)
                diffs = r_sc[:, None] - s_sc[None, :]
                c = float(np.sum(diffs > 0) + 0.5 * np.sum(diffs == 0))
                tot_c += c
                tot_p += npairs
                per_study.append(c / npairs)
        s0 = (tot_c / tot_p) if tot_p > 0 else np.nan
        s_mean = float(np.mean(per_study)) if per_study else np.nan
        return s0, s_mean, tot_p, len(per_study)

    # 4. Factorial Matrix per Drug
    max_d_A = max(d_A.values())
    max_d_B_sh = max(d_B_shell19.values())
    max_d_B_old = max(d_B_old.values())

    results = {
        "status": "COMPLETED",
        "description": "Apportionment of reference set specification (distal MG) vs imputation policy (complete cases vs -max_dist)",
        "drugs": {}
    }

    print("\nFACTORIAL APPORTIONMENT TABLE (S_0 Pair-Weighted AUROC):")
    print("%-12s %-32s %-20s %8s %8s %8s" % ("Drug", "Condition", "Reference Set", "S_0", "Mean", "N_pairs"))
    print("-" * 92)

    for d in DRUGS:
        sub = hcmv_df.copy()
        sub["y"] = [is_res(x) for x in sub[d]]
        sub = sub[sub["y"].notna()].copy()
        sub["y"] = sub["y"].astype(int)
        sub = sub.drop_duplicates(subset=["aa_change", "study_id"]).copy().reset_index(drop=True)

        # Baseline: Option A (Swiss-Model + shell19) Full Cohort
        sub["sc_A_full"] = [- d_A.get(p, max_d_A) for p in sub["pos"]]
        s0_A_full, mean_A_full, p_A_full, m_A_full = calc_study_stats(sub, "sc_A_full")

        # Baseline: Option A on Common Positions Only (N=133)
        sub_common = sub[sub["pos"].isin(d_B_shell19.keys())].copy()
        sub_common["sc_A_common"] = [- d_A[p] for p in sub_common["pos"]]
        s0_A_common, mean_A_common, p_A_common, m_A_common = calc_study_stats(sub_common, "sc_A_common")

        # Cell 1: Clean Shell19 + Complete Cases (Pure Provenance Test)
        sub_common["sc_B_clean_common"] = [- d_B_shell19[p] for p in sub_common["pos"]]
        s0_cell1, mean_cell1, p_cell1, m_cell1 = calc_study_stats(sub_common, "sc_B_clean_common")

        # Cell 2: Clean Shell19 + Imputed (Missingness Effect on Clean Shell)
        sub["sc_B_clean_imputed"] = [- d_B_shell19.get(p, max_d_B_sh) for p in sub["pos"]]
        s0_cell2, mean_cell2, p_cell2, m_cell2 = calc_study_stats(sub, "sc_B_clean_imputed")

        # Cell 3: Distal MG + Triad + Complete Cases (Reference Effect Alone, No Imputation)
        sub_common["sc_B_old_common"] = [- d_B_old[p] for p in sub_common["pos"]]
        s0_cell3, mean_cell3, p_cell3, m_cell3 = calc_study_stats(sub_common, "sc_B_old_common")

        # Cell 4: Distal MG + Triad + Imputed (Original Buggy Option B)
        sub["sc_B_old_imputed"] = [- d_B_old.get(p, max_d_B_old) for p in sub["pos"]]
        s0_cell4, mean_cell4, p_cell4, m_cell4 = calc_study_stats(sub, "sc_B_old_imputed")

        # Effects Apportionment
        total_divergence = s0_cell4 - s0_A_full
        provenance_effect = s0_cell1 - s0_A_common
        reference_effect_clean_cases = s0_cell3 - s0_cell1
        imputation_effect_clean_shell = s0_cell2 - s0_cell1
        imputation_effect_buggy_shell = s0_cell4 - s0_cell3

        print("%-12s %-32s %-20s %8.4f %8.4f %8d" % (d, "Option A (Swiss-Model)", "Shell19 (Full)", s0_A_full, mean_A_full, p_A_full))
        print("%-12s %-32s %-20s %8.4f %8.4f %8d" % (d, "Option A (Swiss-Model)", "Shell19 (Common N=133)", s0_A_common, mean_A_common, p_A_common))
        print("%-12s %-32s %-20s %8.4f %8.4f %8d" % (d, "Cell 1: B_harmonized (Clean)", "Shell19 (Common N=133)", s0_cell1, mean_cell1, p_cell1))
        print("%-12s %-32s %-20s %8.4f %8.4f %8d" % (d, "Cell 2: B_clean + Imputed", "Shell19 (Imputed N=180)", s0_cell2, mean_cell2, p_cell2))
        print("%-12s %-32s %-20s %8.4f %8.4f %8d" % (d, "Cell 3: B_distal_MG (Common)", "Distal MG+Triad (Common)", s0_cell3, mean_cell3, p_cell3))
        print("%-12s %-32s %-20s %8.4f %8.4f %8d" % (d, "Cell 4: Old Option B (Original)", "Distal MG+Triad (Imputed)", s0_cell4, mean_cell4, p_cell4))
        print("  -> Provenance Delta (Cell 1 - A_common) : %+.4f" % provenance_effect)
        print("  -> Reference Set Shift (Cell 3 - Cell 1): %+.4f" % reference_effect_clean_cases)
        print("  -> Imputation Shift (Cell 4 - Cell 3)   : %+.4f" % imputation_effect_buggy_shell)
        print("  -> Total Divergence (Old B - Opt A)     : %+.4f\n" % total_divergence)

        results["drugs"][d] = {
            "option_a_full_cohort": {
                "s0_pair_weighted": round(s0_A_full, 4),
                "study_mean": round(mean_A_full, 4),
                "n_pairs": p_A_full,
                "m_studies": m_A_full
            },
            "option_a_common_cases": {
                "s0_pair_weighted": round(s0_A_common, 4),
                "study_mean": round(mean_A_common, 4),
                "n_pairs": p_A_common,
                "m_studies": m_A_common
            },
            "cell1_b_clean_shell19_common": {
                "s0_pair_weighted": round(s0_cell1, 4),
                "study_mean": round(mean_cell1, 4),
                "n_pairs": p_cell1,
                "m_studies": m_cell1
            },
            "cell2_b_clean_shell19_imputed": {
                "s0_pair_weighted": round(s0_cell2, 4),
                "study_mean": round(mean_cell2, 4),
                "n_pairs": p_cell2,
                "m_studies": m_cell2
            },
            "cell3_b_distal_mg_triad_common": {
                "s0_pair_weighted": round(s0_cell3, 4),
                "study_mean": round(mean_cell3, 4),
                "n_pairs": p_cell3,
                "m_studies": m_cell3
            },
            "cell4_old_option_b_imputed": {
                "s0_pair_weighted": round(s0_cell4, 4),
                "study_mean": round(mean_cell4, 4),
                "n_pairs": p_cell4,
                "m_studies": m_cell4
            },
            "effects_apportionment": {
                "pure_provenance_delta": round(provenance_effect, 4),
                "reference_set_shift_common": round(reference_effect_clean_cases, 4),
                "imputation_shift_buggy_shell": round(imputation_effect_buggy_shell, 4),
                "total_divergence": round(total_divergence, 4)
            }
        }

    out_path = os.path.join(ROOT, "logs", "option_b_harmonization_and_apportionment_results.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print("Wrote results to: %s" % os.path.relpath(out_path, ROOT))

if __name__ == "__main__":
    main()

