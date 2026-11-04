#!/usr/bin/env python3
"""
dump_ul23_scope_contingency.py
==============================
Computes and dumps the exact clinical scope contingency table for UL23
from data/processed/unified_modeling_dataset.tsv to logs/ul23_clinical_scope_contingency.json.

This settles the exact provenance of:
  - 186 / 419 (44.4%) non-missense resistant variants across the unified UL23 database
  - 233 / 419 (55.6%) missense resistant variants
  - 98 / 255 (38.4%) non-missense resistant variants under aciclovir-specific filtering
"""

import os
import json
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_PATH = os.path.join(ROOT, "data", "processed", "unified_modeling_dataset.tsv")
OUT_PATH = os.path.join(ROOT, "logs", "ul23_clinical_scope_contingency.json")

def main():
    df = pd.read_csv(DATA_PATH, sep="\t")
    ul23 = df[df["gene"] == "UL23"].copy()
    
    # 1. Unified UL23 (all drugs cataloged)
    total_rows = len(ul23)
    evaluable = ul23[ul23["target_binary"].notna()].copy()
    resistant_all = evaluable[evaluable["target_binary"] == 1.0].copy()
    
    all_drugs_breakdown = resistant_all["mutation_category"].value_counts().to_dict()
    n_res_all = len(resistant_all)
    n_trunc_all = int(resistant_all["mutation_category"].isin(["nonsense_stop", "frameshift", "indel"]).sum())
    n_missense_all = int((resistant_all["mutation_category"] == "single_substitution").sum())
    
    # 2. UL23 Aciclovir-specific subset
    ul23_acv = ul23[ul23["drug"].str.upper().str.contains("ACV|ACICLOVIR|ACYCLOVIR")].copy()
    evaluable_acv = ul23_acv[ul23_acv["target_binary"].notna()].copy()
    resistant_acv = evaluable_acv[evaluable_acv["target_binary"] == 1.0].copy()
    
    acv_breakdown = resistant_acv["mutation_category"].value_counts().to_dict()
    n_res_acv = len(resistant_acv)
    n_trunc_acv = int(resistant_acv["mutation_category"].isin(["nonsense_stop", "frameshift", "indel"]).sum())
    n_missense_acv = int((resistant_acv["mutation_category"] == "single_substitution").sum())
    
    payload = {
        "purpose": "Provenance and exact contingency counts for UL23 clinical scope quantification",
        "source_dataset": "data/processed/unified_modeling_dataset.tsv",
        "total_ul23_entries": total_rows,
        "evaluable_ul23_entries": len(evaluable),
        "all_drugs_ul23_catalog": {
            "description": "All cataloged resistant UL23 entries across all 7 antiviral drugs (aciclovir, foscarnet, penciclovir, brivudine, cidofovir, pritelivir, amenamevir)",
            "total_resistant_denominator": n_res_all,
            "non_missense_numerator": n_trunc_all,
            "missense_count": n_missense_all,
            "pct_non_missense": round(100.0 * n_trunc_all / n_res_all, 1),
            "pct_missense": round(100.0 * n_missense_all / n_res_all, 1),
            "exact_quotient_non_missense": f"{n_trunc_all} / {n_res_all}",
            "exact_quotient_missense": f"{n_missense_all} / {n_res_all}",
            "mutation_category_counts": all_drugs_breakdown
        },
        "aciclovir_specific_ul23_catalog": {
            "description": "UL23 entries filtered strictly to aciclovir assays",
            "total_resistant_denominator": n_res_acv,
            "non_missense_numerator": n_trunc_acv,
            "missense_count": n_missense_acv,
            "pct_non_missense": round(100.0 * n_trunc_acv / n_res_acv, 1),
            "pct_missense": round(100.0 * n_missense_acv / n_res_acv, 1),
            "exact_quotient_non_missense": f"{n_trunc_acv} / {n_res_acv}",
            "exact_quotient_missense": f"{n_missense_acv} / {n_res_acv}",
            "mutation_category_counts": acv_breakdown
        }
    }
    
    with open(OUT_PATH, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
        
    print(f"Wrote {OUT_PATH}")
    print(f"All-drugs UL23: non-missense = {n_trunc_all} / {n_res_all} ({payload['all_drugs_ul23_catalog']['pct_non_missense']}%)")
    print(f"Aciclovir-specific UL23: non-missense = {n_trunc_acv} / {n_res_acv} ({payload['aciclovir_specific_ul23_catalog']['pct_non_missense']}%)")

if __name__ == "__main__":
    main()

