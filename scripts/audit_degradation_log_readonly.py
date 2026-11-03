"""
audit_degradation_log_readonly.py
=================================
READ-ONLY audit of logs/split_degradation_benchmark_results.json.

No model is fitted, no resampling is performed, no log is written. This script only
prints what the frozen log does and does not contain, so that three proposed findings
can be checked before any of them is written into the manuscript.

It answers four questions:

  Q1  Does the log carry a confidence interval for the 3D k-means column?
      If not, any claim comparing the 3D AUROC against the position-grouped AUROC is
      being made without an uncertainty estimate.

  Q2  Under the 3D split, is the evolutionary-conservation family also flat?
      If it is, "PLM is split-agnostic" is not a PLM-specific property.

  Q3  Where does the unfitted BLOSUM62 negative control rank if families are ordered by
      smallest degradation delta?
      Its delta is 0.0000 by construction, so if it ranks first then "smallest delta"
      is confounded with being uninformative and cannot serve as a robustness ranking.

  Q4  Ranked by absolute contiguous-block AUROC instead - out-of-block performance,
      which is not confounded with informativeness - what is the ordering, and what is
      the PLM-minus-conservation margin on each arm?
"""

import json
import os

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
LOG = os.path.join(ROOT, "logs", "split_degradation_benchmark_results.json")

FAMILIES = ["Positional Identity", "Evolutionary Conservation",
            "Protein Language Models", "Multi-Modal Combined"]


def main():
    with open(LOG, "r", encoding="utf-8") as f:
        d = json.load(f)
    targets = d["targets"]
    print("=" * 84)
    print("READ-ONLY AUDIT: logs/split_degradation_benchmark_results.json")
    print("=" * 84)

    # ---- Q1: is there any interval on the 3D column? -------------------------
    print("")
    print("Q1  Confidence intervals present per column")
    print("-" * 84)
    sample = targets["UL30-Aciclovir"]["models"]["random_forest"]["Protein Language Models"]
    keys = sorted(sample.keys())
    print("    keys on one family record:")
    for k in keys:
        print("      %s" % k)
    has_3d_ci = any("spatial" in k and "ci" in k for k in keys)
    print("")
    print("    position_grouped_ci_95 present : %s" % ("position_grouped_ci_95" in keys))
    print("    contiguous_block_ci_95 present : %s" % ("contiguous_block_ci_95" in keys))
    print("    spatial_3d  CI present         : %s" % has_3d_ci)
    if not has_3d_ci:
        print("    => The 3D k-means column carries NO interval. A claim that PLM AUROC is")
        print("       'essentially unchanged' from position-grouped to 3D cannot be")
        print("       supported from this log, in either direction.")

    # ---- Q2: is conservation also flat under the 3D split? ------------------
    print("")
    print("Q2  Change from position-grouped to 3D k-means, PLM vs conservation")
    print("-" * 84)
    print("    %-26s %-6s %10s %10s %10s" % ("family", "model", "pos-grp", "3D", "change"))
    for t in targets:
        for m_id, m_lab in [("random_forest", "RF"), ("logistic_regression", "LR")]:
            mods = targets[t]["models"].get(m_id, {})
            for fam in ["Evolutionary Conservation", "Protein Language Models"]:
                if fam not in mods:
                    continue
                r = mods[fam]
                pg, sp = r["position_grouped_auc"], r.get("spatial_3d_auc")
                if sp is None:
                    continue
                print("    %-26s %-6s %10.4f %10.4f %+10.4f"
                      % (fam[:26] + " | " + t[:11], m_lab, pg, sp, sp - pg))
    print("")
    print("    => If conservation changes by a comparable amount, flatness under the 3D")
    print("       split is not specific to language models.")

    # ---- Q3: where does the zero-delta negative control rank? ---------------
    print("")
    print("Q3  Ranking by SMALLEST degradation delta, negative control included")
    print("-" * 84)
    for t in targets:
        nc = targets[t].get("negative_control_zero_shot_blosum", {})
        nc_delta = nc.get("degradation_delta")
        for m_id, m_lab in [("random_forest", "RF"), ("logistic_regression", "LR")]:
            mods = targets[t]["models"].get(m_id, {})
            rows = [(fam, mods[fam]["degradation_delta_pg_minus_block"])
                    for fam in FAMILIES if fam in mods]
            if nc_delta is not None:
                rows.append(("NEGATIVE CONTROL BLOSUM62 (unfitted)", nc_delta))
            rows.sort(key=lambda x: abs(x[1]))
            print("    %s / %s" % (t, m_lab))
            for i, (fam, dl) in enumerate(rows, 1):
                flag = "  <== unfitted control ranks here" if "NEGATIVE" in fam else ""
                print("       %d. %-42s delta = %+.4f%s" % (i, fam[:42], dl, flag))
    print("")
    print("    => A score that is never fitted has delta 0.0000 by construction. If it")
    print("       ranks first, 'smallest delta' measures insensitivity, not robustness,")
    print("       and must not be reported as a robustness ranking.")

    # ---- Q4: absolute out-of-block performance ------------------------------
    print("")
    print("Q4  Ranking by ABSOLUTE contiguous-block AUROC (not confounded)")
    print("-" * 84)
    margins = []
    for t in targets:
        for m_id, m_lab in [("random_forest", "RF"), ("logistic_regression", "LR")]:
            mods = targets[t]["models"].get(m_id, {})
            rows = [(fam, mods[fam]["contiguous_block_auc"],
                     mods[fam].get("contiguous_block_ci_95"))
                    for fam in FAMILIES if fam in mods]
            rows.sort(key=lambda x: -x[1])
            print("    %s / %s" % (t, m_lab))
            for i, (fam, auc, ci) in enumerate(rows, 1):
                cis = (" [%.4f, %.4f]" % (ci[0], ci[1])) if ci else " [no CI in log]"
                print("       %d. %-30s block AUROC = %.4f%s" % (i, fam[:30], auc, cis))
            cons = mods.get("Evolutionary Conservation", {}).get("contiguous_block_auc")
            plm = mods.get("Protein Language Models", {}).get("contiguous_block_auc")
            if cons is not None and plm is not None:
                margins.append((t, m_lab, cons, plm, plm - cons))

    print("")
    print("    PLM minus conservation, out-of-block:")
    print("    %-22s %-6s %10s %10s %10s" % ("target", "model", "cons", "PLM", "margin"))
    for t, m, c, p, mg in margins:
        print("    %-22s %-6s %10.4f %10.4f %+10.4f" % (t, m, c, p, mg))
    wins = sum(1 for *_, mg in margins if mg > 0)
    print("")
    print("    PLM ahead in %d of %d combinations; margin range %+.4f to %+.4f"
          % (wins, len(margins), min(m[-1] for m in margins), max(m[-1] for m in margins)))
    print("")
    print("    NOTE: these block-AUROC intervals are NOT paired, so overlapping intervals")
    print("    do not settle the comparison. The paired test is in the companion script")
    print("    scripts/run_paired_family_contrast.py, which requires a re-run.")
    print("")
    print("Nothing was fitted, resampled or written. END.")


if __name__ == "__main__":
    main()
