"""
run_ul54_leave_site_out_robustness.py
=====================================
Robustness check for the UL54 exonuclease / polymerase double dissociation.

Concern: the 7-residue exonuclease reference set {301,302,303,407,408,412,413} contains
positions that are themselves frequent resistance sites (D301N, N408D/K/S, D413E/N).
A variant AT a reference residue has distance 0 by construction, so S_0(exo) could be
carried by a handful of hot positions rather than by a spatial gradient.

The same concern applies symmetrically to the 19-residue polymerase shell, so both sides
are tested identically. Four cohorts per drug:

  full              all variants (reproduces the confirmatory run)
  drop_exo_sites    variants at the 7 exonuclease reference residues removed
  drop_pol_sites    variants at the 19 polymerase reference residues removed
  drop_both         both reference sets removed

If S_0(exo) on GCV/CDV and S_0(pol) on FOS survive removal of their own defining
residues, the signal is a genuine spatial gradient. If a value collapses toward 0.5, that
arm was driven by the defining positions themselves and must be reported as such.

Estimator, phenotyping, study assignment, seed and B are identical to
scripts/run_ul54_exonuclease_confirmatory.py so the numbers are directly comparable.

Output: logs/ul54_leave_site_out_robustness.json
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

EXO_SITES = [301, 302, 303, 407, 408, 412, 413]
DRUGS = ["Ganciclovir", "Foscarnet", "Cidofovir"]
N_BOOT = 2000
SEED = 42


def phenotype(val, thresh=2.0):
    s = str(val).strip()
    if s == "Resistant":
        return 1
    if s in ("Polymorphism", "Sensitive"):
        return 0
    try:
        return 1 if float(s) >= thresh else 0
    except ValueError:
        return None


def study_id(row):
    doi = str(row["ref_doi"]).strip()
    if doi and doi not in ("nan", "None", ""):
        return "DOI:%s" % doi
    t = str(row["ref_title"]).strip()
    if t and t not in ("nan", "None", ""):
        return "TITLE:%s" % t[:45]
    return "Unassigned"


def study_pair_stats(df, score_col):
    stats = {}
    for st, grp in df.groupby("study_id"):
        r = grp.loc[grp.y == 1, score_col].values
        s = grp.loc[grp.y == 0, score_col].values
        if len(r) and len(s):
            diff = r[:, None] - s[None, :]
            stats[st] = {"n_pairs": len(r) * len(s),
                         "concordant": float(np.sum(diff > 0) + 0.5 * np.sum(diff == 0))}
    return stats


def pair_weighted(stats):
    tp = sum(v["n_pairs"] for v in stats.values())
    tc = sum(v["concordant"] for v in stats.values())
    return (tc / tp) if tp else np.nan, tp


def main():
    parser = PDBParser(QUIET=True)
    model = list(parser.get_structure(
        "UL54", os.path.join(ROOT, "data", "structures",
                             "P08546_swissmodel_8oj7.pdb")).get_models())[0]
    with open(os.path.join(ROOT, "logs", "hcmv_ul54_active_site_mapping.json")) as f:
        pol_shell = sorted(set(json.load(f)["hcmv_active_site_positions"]))

    ca, pol_atoms, exo_atoms = {}, [], []
    for chain in model:
        for res in chain:
            rid = res.id[1]
            if "CA" in res:
                ca[rid] = res["CA"].get_coord()
            if rid in pol_shell:
                pol_atoms += [a.get_coord() for a in res]
            if rid in EXO_SITES:
                exo_atoms += [a.get_coord() for a in res]
    pol_atoms, exo_atoms = np.asarray(pol_atoms), np.asarray(exo_atoms)
    d_pol = {p: float(np.min(np.linalg.norm(pol_atoms - c, axis=1))) for p, c in ca.items()}
    d_exo = {p: float(np.min(np.linalg.norm(exo_atoms - c, axis=1))) for p, c in ca.items()}
    max_pol, max_exo = max(d_pol.values()), max(d_exo.values())

    print("=" * 80)
    print("LEAVE-SITE-OUT ROBUSTNESS: UL54 exonuclease / polymerase double dissociation")
    print("=" * 80)
    print("polymerase shell (19): %s" % pol_shell)
    print("exonuclease site  (7): %s" % EXO_SITES)

    raw = pd.read_csv(os.path.join(ROOT, "data", "raw", "herpesdrg-db", "herpesdrg-db.tsv"),
                      sep="\t")
    h = raw[(raw.virus == "HCMV") & (raw.gene == "UL54")].copy()
    m = [re.match(r"^([A-Z])(\d+)([A-Z])$", str(s).strip()) for s in h.aa_change]
    h["pos"] = [int(x.group(2)) if x else np.nan for x in m]
    h = h[h.pos.notna()].copy()
    h["pos"] = h.pos.astype(int)
    h["study_id"] = [study_id(r) for _, r in h.iterrows()]

    COHORTS = {
        "full": lambda s: s,
        "drop_exo_sites": lambda s: s[~s.pos.isin(EXO_SITES)],
        "drop_pol_sites": lambda s: s[~s.pos.isin(pol_shell)],
        "drop_both": lambda s: s[~s.pos.isin(EXO_SITES) & ~s.pos.isin(pol_shell)],
    }

    out = {
        "purpose": ("robustness of S_0(exo) and S_0(pol) to removal of the variants that "
                    "sit on the reference residues defining each feature"),
        "estimator": "within-study pair-weighted S_0, study-clustered bootstrap",
        "n_boot": N_BOOT, "seed": SEED,
        "polymerase_shell_19": pol_shell, "exonuclease_site_7": EXO_SITES,
        "drugs": {}}

    rng = np.random.RandomState(SEED)
    for d in DRUGS:
        sub0 = h.copy()
        sub0["y"] = [phenotype(x) for x in sub0[d]]
        sub0 = sub0[sub0.y.notna()].copy()
        sub0["y"] = sub0.y.astype(int)
        sub0 = sub0.drop_duplicates(subset=["aa_change", "study_id"]).reset_index(drop=True)
        sub0["score_pol"] = [-d_pol.get(p, max_pol) for p in sub0.pos]
        sub0["score_exo"] = [-d_exo.get(p, max_exo) for p in sub0.pos]

        n_at_exo = int(sub0.pos.isin(EXO_SITES).sum())
        n_at_pol = int(sub0.pos.isin(pol_shell).sum())
        r_at_exo = int(sub0[(sub0.y == 1) & sub0.pos.isin(EXO_SITES)].shape[0])
        r_at_pol = int(sub0[(sub0.y == 1) & sub0.pos.isin(pol_shell)].shape[0])

        print("")
        print("== %s   N=%d  R=%d  S=%d" % (d, len(sub0), int(sub0.y.sum()),
                                            int((sub0.y == 0).sum())))
        print("   rows at exo reference residues: %d (of which R=%d)" % (n_at_exo, r_at_exo))
        print("   rows at pol reference residues: %d (of which R=%d)" % (n_at_pol, r_at_pol))
        print("   %-16s %6s %6s %10s %10s %24s" %
              ("cohort", "N", "pairs", "S0_pol", "S0_exo", "delta exo-pol [CI]"))

        rec = {"n_total": len(sub0), "n_resistant": int(sub0.y.sum()),
               "n_rows_at_exo_sites": n_at_exo, "n_resistant_at_exo_sites": r_at_exo,
               "n_rows_at_pol_sites": n_at_pol, "n_resistant_at_pol_sites": r_at_pol,
               "cohorts": {}}

        for cname, filt in COHORTS.items():
            sub = filt(sub0).copy()
            st_pol = study_pair_stats(sub, "score_pol")
            st_exo = study_pair_stats(sub, "score_exo")
            s0p, npairs = pair_weighted(st_pol)
            s0e, _ = pair_weighted(st_exo)

            studies = sorted(set(st_pol) | set(st_exo))
            deltas, s0ps, s0es = [], [], []
            for _ in range(N_BOOT):
                samp = rng.choice(studies, size=len(studies), replace=True)
                tp = tc = tpe = tce = 0.0
                for st in samp:
                    a, b = st_pol.get(st), st_exo.get(st)
                    if a:
                        tp += a["n_pairs"]; tc += a["concordant"]
                    if b:
                        tpe += b["n_pairs"]; tce += b["concordant"]
                if tp > 0 and tpe > 0:
                    bp, be = tc / tp, tce / tpe
                    s0ps.append(bp); s0es.append(be); deltas.append(be - bp)
            ci = [float(np.percentile(deltas, 2.5)), float(np.percentile(deltas, 97.5))] \
                if deltas else [float("nan")] * 2
            pval = float(np.mean(np.asarray(deltas) <= 0)) if deltas else float("nan")

            print("   %-16s %6d %6d %10.4f %10.4f   %+.4f [%+.4f,%+.4f]" %
                  (cname, len(sub), npairs, s0p, s0e, s0e - s0p, ci[0], ci[1]))

            rec["cohorts"][cname] = {
                "n_rows": len(sub), "n_pairs": int(npairs),
                "m_mixed_studies": len(st_pol),
                "s0_polymerase": round(float(s0p), 4),
                "s0_exonuclease": round(float(s0e), 4),
                "s0_polymerase_ci_95": [round(float(np.percentile(s0ps, 2.5)), 4),
                                        round(float(np.percentile(s0ps, 97.5)), 4)] if s0ps else None,
                "s0_exonuclease_ci_95": [round(float(np.percentile(s0es, 2.5)), 4),
                                         round(float(np.percentile(s0es, 97.5)), 4)] if s0es else None,
                "delta_exo_minus_pol": round(float(s0e - s0p), 4),
                "delta_ci_95": [round(ci[0], 4), round(ci[1], 4)],
                "delta_p_one_sided": round(pval, 4)}
        out["drugs"][d] = rec

    # verdict: does each arm's own signal survive removal of its defining residues
    v = {}
    for d in DRUGS:
        c = out["drugs"][d]["cohorts"]
        key = "s0_exonuclease" if d in ("Ganciclovir", "Cidofovir") else "s0_polymerase"
        full, dropped = c["full"][key], c["drop_both"][key]
        v[d] = {"primary_feature": key,
                "s0_full": full, "s0_drop_both": dropped,
                "retained_delta": round(dropped - full, 4),
                "survives_above_0.55": bool(dropped >= 0.55),
                "ci_excludes_0.50": bool(
                    (c["drop_both"][key + "_ci_95"] or [0, 0])[0] > 0.50)}
    out["survival_verdicts"] = v
    print("")
    print("SURVIVAL (own defining residues removed):")
    for d, r in v.items():
        print("  %-12s %-16s %.4f -> %.4f   >=0.55: %-5s  CI excl 0.50: %s"
              % (d, r["primary_feature"], r["s0_full"], r["s0_drop_both"],
                 r["survives_above_0.55"], r["ci_excludes_0.50"]))

    p = os.path.join(ROOT, "logs", "ul54_leave_site_out_robustness.json")
    with open(p, "w") as f:
        json.dump(out, f, indent=1)
    print("")
    print("Wrote %s" % os.path.relpath(p, ROOT))


if __name__ == "__main__":
    main()
