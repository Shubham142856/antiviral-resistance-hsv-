"""
explore_ul54_drug_localization.py
=================================
EXPLORATORY (POST HOC) analysis of the already-executed, frozen Step 2 HCMV UL54 arm.

Status: this is NOT a pre-registered test. The Step 2 primary arm (UL54-Ganciclovir
active-site distance) failed its pre-registered criterion at S_0 = 0.4135, and that
verdict stands unchanged. The analysis below asks why, and is reported as exploratory.

Question: the single active-site distance feature gives S_0 = 0.8822 on Foscarnet,
0.4135 on Ganciclovir and 0.2507 on Cidofovir. The distance vector is identical for all
three drugs, so the swing must come from where each drug's resistant variants sit.

Premise tested here: Foscarnet resistance localizes to the polymerase catalytic core,
while Ganciclovir and Cidofovir resistance localizes to the 3'->5' exonuclease
proofreading region, which is far from the polymerase active site.

Phenotyping follows the frozen Step 2 pre-registration exactly: fold-change >= 2.0 or a
reported 'Resistant' label => resistant; fold-change < 2.0 or 'Polymorphism'/'Sensitive'
=> sensitive. Cohort counts are asserted against the frozen log.
"""

import json
import os
import re
import warnings

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser

warnings.filterwarnings("ignore")
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

EXO_REGION = (300, 545)      # Exo I / II / III, HCMV UL54
POL_CORE = (690, 960)        # polymerase catalytic core spanning the mapped active site
DRUGS = ["Foscarnet", "Ganciclovir", "Cidofovir"]

# frozen distance-only S_0 (Option A), logs/step2_hcmv_composite_results.json
FROZEN_S0 = {"Foscarnet": 0.8822, "Ganciclovir": 0.4135, "Cidofovir": 0.2507}
# frozen cohort sizes, logs/cross_virus_replication_results.json
FROZEN_N = {"Ganciclovir": (230, 76, 154), "Foscarnet": (160, 50, 110),
            "Cidofovir": (149, 57, 92)}


def phenotype(v):
    s = str(v).strip()
    if s in ("nan", "", "None"):
        return np.nan
    low = s.lower()
    if low.startswith("resist"):
        return 1
    if low.startswith(("polymorph", "sensit")):
        return 0
    try:
        return 1 if float(s) >= 2.0 else 0
    except ValueError:
        return np.nan


def active_site_distances():
    """Option A distances: Swiss-Model UL54 C-alpha -> 19-residue active-site shell.
    These coordinates are validated against experimental 8EXX at rho = 0.989,
    median |delta| = 0.39 A (logs/option_ab_distance_agreement.json)."""
    with open(os.path.join(ROOT, "logs", "hcmv_ul54_active_site_mapping.json")) as f:
        shell = sorted(set(json.load(f)["hcmv_active_site_positions"]))
    model = list(PDBParser(QUIET=True).get_structure(
        "UL54", os.path.join(ROOT, "data", "structures",
                             "P08546_swissmodel_8oj7.pdb")).get_models())[0]
    ca, ref = {}, []
    for chain in model:
        for res in chain:
            if "CA" in res:
                ca[res.id[1]] = res["CA"].get_coord()
            if res.id[1] in shell:
                ref += [a.get_coord() for a in res]
    ref = np.asarray(ref)
    return {p: float(np.min(np.linalg.norm(ref - c, axis=1))) for p, c in ca.items()}


def main():
    print("=" * 78)
    print("EXPLORATORY: UL54 drug-dependent structural localization")
    print("Pre-registered Step 2 primary verdict (GCV S_0 = 0.4135, NOT SUPPORTED)")
    print("is unchanged by this analysis.")
    print("=" * 78)

    dist = active_site_distances()
    exo_d = [v for p, v in dist.items() if EXO_REGION[0] <= p <= EXO_REGION[1]]
    core_d = [v for p, v in dist.items() if POL_CORE[0] <= p <= POL_CORE[1]]
    print("")
    print("Geometric separation of the two catalytic regions:")
    print("  exo %d-%d      : n=%d  median %.1f A to polymerase shell"
          % (EXO_REGION[0], EXO_REGION[1], len(exo_d), np.median(exo_d)))
    print("  pol core %d-%d : n=%d  median %.1f A to polymerase shell"
          % (POL_CORE[0], POL_CORE[1], len(core_d), np.median(core_d)))
    print("  separation of medians: %.1f A" % (np.median(exo_d) - np.median(core_d)))

    raw = pd.read_csv(os.path.join(ROOT, "data", "raw", "herpesdrg-db", "herpesdrg-db.tsv"),
                      sep="\t")
    u = raw[(raw.virus == "HCMV") & (raw.gene == "UL54")].copy()
    m = [re.match(r"^([A-Z])(\d+)([A-Z])$", str(s).strip()) for s in u.aa_change]
    u["pos"] = [int(x.group(2)) if x else np.nan for x in m]
    u = u[u.pos.notna()].copy()
    u["pos"] = u.pos.astype(int)

    print("")
    print("%-13s %5s %5s %5s %12s %12s %10s %10s %9s"
          % ("drug", "N", "R", "S", "%R in exo", "%S in exo",
             "medR pos", "medS pos", "dist S_0"))
    out = {
        "status": ("EXPLORATORY / POST HOC on the frozen Step 2 HCMV arm; "
                   "pre-registered GCV verdict S_0=0.4135 NOT SUPPORTED is unchanged"),
        "phenotyping_rule": "fold-change >= 2.0 or 'Resistant' => R (frozen Step 2 rule)",
        "exo_region": list(EXO_REGION),
        "pol_core_region": list(POL_CORE),
        "coordinate_source": ("data/structures/P08546_swissmodel_8oj7.pdb, validated vs "
                              "8EXX at spearman 0.989 / median 0.39 A"),
        "geometric_separation": {
            "exo_median_distance_to_pol_shell_angstrom": round(float(np.median(exo_d)), 1),
            "pol_core_median_distance_to_pol_shell_angstrom": round(float(np.median(core_d)), 1),
            "separation_of_medians_angstrom": round(float(np.median(exo_d) - np.median(core_d)), 1)},
        "drugs": {}}

    for d in DRUGS:
        s = u.copy()
        s["y"] = [phenotype(v) for v in s[d]]
        s = s[s.y.notna()].drop_duplicates(subset=["aa_change"])
        R, S = s[s.y == 1], s[s.y == 0]
        r_exo = int(R.pos.between(*EXO_REGION).sum())
        s_exo = int(S.pos.between(*EXO_REGION).sum())
        pr, ps = 100.0 * r_exo / len(R), 100.0 * s_exo / len(S)
        print("%-13s %5d %5d %5d %11.1f%% %11.1f%% %10.1f %10.1f %9.4f"
              % (d, len(s), len(R), len(S), pr, ps, R.pos.median(), S.pos.median(),
                 FROZEN_S0[d]))
        exp_n = FROZEN_N[d]
        match = (len(s), len(R), len(S)) == exp_n
        if not match:
            print("   WARNING cohort mismatch vs frozen log: got %s expected %s"
                  % ((len(s), len(R), len(S)), exp_n))
        out["drugs"][d] = {
            "n_total": len(s), "n_resistant": len(R), "n_sensitive": len(S),
            "matches_frozen_cohort": bool(match),
            "frozen_cohort_expected": list(exp_n),
            "pct_resistant_in_exo": round(pr, 1),
            "pct_sensitive_in_exo": round(ps, 1),
            "n_resistant_in_exo": r_exo, "n_sensitive_in_exo": s_exo,
            "median_pos_resistant": float(R.pos.median()),
            "median_pos_sensitive": float(S.pos.median()),
            "frozen_distance_only_s0_option_a": FROZEN_S0[d]}

    order = sorted(DRUGS, key=lambda d: out["drugs"][d]["pct_resistant_in_exo"])
    s0s = [FROZEN_S0[d] for d in order]
    monotone = all(s0s[i] > s0s[i + 1] for i in range(len(s0s) - 1))
    out["monotonicity_check"] = {
        "drugs_sorted_by_pct_resistant_in_exo": order,
        "corresponding_distance_s0": s0s,
        "s0_strictly_decreasing": bool(monotone),
        "n_points": len(order),
        "caveat": ("3 drugs only; ordering is exact but no inferential claim is made, "
                   "and the GCV S_0 95% CI [0.1901, 0.6941] is wide")}
    print("")
    print("Monotonicity: %s -> S_0 %s  strictly decreasing = %s  (n=3, descriptive only)"
          % (order, [round(x, 4) for x in s0s], monotone))

    log_path = os.path.join(ROOT, "logs", "ul54_drug_localization_exploratory.json")
    with open(log_path, "w") as f:
        json.dump(out, f, indent=1)
    print("")
    print("Wrote %s" % os.path.relpath(log_path, ROOT))


if __name__ == "__main__":
    main()
