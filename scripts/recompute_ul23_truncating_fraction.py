"""
recompute_ul23_truncating_fraction.py
=====================================
Recomputes the truncating fraction of cataloged UL23 aciclovir resistance.

Why. The manuscript states "44.4% (186/419) of cataloged resistance-conferring variants
in thymidine kinase (UL23)" in four places, including the abstract. Neither 186 nor 419
appears in any file under logs/. The only occurrence of 44.4 in the logs is
stage_b2_study_counts.json -> pct_resistant 44.4444 for a single study (4 resistant of
9 total), which is an unrelated quantity. logs/ul23_label_composition.json records
nonsense_stop 168 and frameshift 68 over a denominator of 941 total UL23 rows, and it
does not stratify by resistance at all, so it cannot produce the stated figure either.

This script prints the fraction under every defensible specification, so that one can be
pre-specified, cited, and the unsupported figure replaced. It writes a log; it changes no
manuscript text.

The specification axes are:
  source      HerpesDRG only, or HerpesDRG merged with the CMI 2025 Table S3 supplement
  threshold   fold-change cutoff for calling resistance (2.0x and 3.0x are both used
              elsewhere in this project, so both are reported)
  dedup       all rows, or deduplicated on distinct amino-acid change
  indels      whether insertion/deletion records are counted as truncating

Output: logs/ul23_truncating_fraction_recomputed.json
"""

import json
import os
import re

import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "logs", "ul23_truncating_fraction_recomputed.json")


def categorize(s):
    """Classify an amino-acid change string. 'other' is kept separate and never
    silently folded into either numerator or denominator category."""
    t = str(s).strip()
    if re.match(r"^[A-Z]\d+[A-Z]$", t):
        return "missense"
    if "*" in t or re.match(r"^[A-Z]\d+X$", t, re.I):
        return "nonsense"
    if re.search(r"fs|del|ins|dup", t, re.I):
        return "indel_fs"
    return "other"


def phenotype(val, thresh):
    s = str(val).strip()
    low = s.lower()
    if low.startswith("resist"):
        return 1
    if low.startswith(("polymorph", "sensit")):
        return 0
    try:
        return 1 if float(s) >= thresh else 0
    except ValueError:
        return np.nan


def load_herpesdrg():
    p = os.path.join(ROOT, "data", "raw", "herpesdrg-db", "herpesdrg-db.tsv")
    d = pd.read_csv(p, sep="\t")
    d = d[d.gene == "UL23"][["aa_change", "Aciclovir"]].copy()
    d["source"] = "HerpesDRG"
    return d


def load_cmi():
    """CMI 2025 Table S3 supplement, if its columns can be identified."""
    p = os.path.join(ROOT, "data", "processed", "cmi_2025_table_s3.csv")
    if not os.path.exists(p):
        return None, "file absent"
    d = pd.read_csv(p)
    gene_c = next((c for c in d.columns if "gene" in c.lower()), None)
    mut_c = next((c for c in d.columns
                  if "mutation" in c.lower() or "aa" in c.lower()
                  or "change" in c.lower() or "substitution" in c.lower()), None)
    acv_c = next((c for c in d.columns
                  if "aciclovir" in c.lower() or "acyclovir" in c.lower()
                  or c.strip().upper() == "ACV"), None)
    if not all([gene_c, mut_c, acv_c]):
        return None, ("columns not identifiable; available = %s" % list(d.columns))
    s = d[d[gene_c].astype(str).str.upper().str.contains("UL23", na=False)]
    s = s[[mut_c, acv_c]].rename(columns={mut_c: "aa_change", acv_c: "Aciclovir"}).copy()
    s["source"] = "CMI_S3"
    return s, "ok"


def summarize(df, thresh, dedup, count_indels):
    d = df.copy()
    d["y"] = [phenotype(v, thresh) for v in d["Aciclovir"]]
    d = d[d.y.notna()].copy()
    if dedup:
        d = d.drop_duplicates(subset=["aa_change"])
    d["cat"] = [categorize(s) for s in d.aa_change]
    R = d[d.y == 1]
    trunc_cats = ["nonsense", "indel_fs"] if count_indels else ["nonsense"]
    n_tr = int(R.cat.isin(trunc_cats).sum())
    n_tot = int(len(R))
    n_ms = int((R.cat == "missense").sum())
    n_other = int((R.cat == "other").sum())
    return {
        "n_resistant_total": n_tot,
        "n_truncating": n_tr,
        "n_missense": n_ms,
        "n_other_unclassified": n_other,
        "pct_truncating": round(100.0 * n_tr / n_tot, 1) if n_tot else None,
        "pct_missense": round(100.0 * n_ms / n_tot, 1) if n_tot else None,
    }


def main():
    hdrg = load_herpesdrg()
    cmi, cmi_status = load_cmi()
    print("=" * 88)
    print("RECOMPUTE: truncating fraction of cataloged UL23 aciclovir resistance")
    print("=" * 88)
    print("HerpesDRG UL23 rows: %d" % len(hdrg))
    print("CMI Table S3 status: %s%s"
          % (cmi_status, ("" if cmi is None else "  (%d UL23 rows)" % len(cmi))))
    print("")
    print("Manuscript currently asserts 44.4%% (186/419). Neither 186 nor 419 appears")
    print("in any file under logs/. Specifications below; none is privileged until one")
    print("is pre-specified.")
    print("")

    sources = {"HerpesDRG only": hdrg}
    if cmi is not None:
        sources["HerpesDRG + CMI S3"] = pd.concat([hdrg, cmi], ignore_index=True)

    results = {}
    hdr = "%-20s %6s %6s %7s %10s %10s %8s %8s"
    print(hdr % ("source", "thresh", "dedup", "indels",
                 "n_resist", "n_trunc", "%trunc", "%miss"))
    print("-" * 88)
    for sname, sdf in sources.items():
        for thresh in (2.0, 3.0):
            for dedup in (False, True):
                for ind in (True, False):
                    r = summarize(sdf, thresh, dedup, ind)
                    key = "%s | FC>=%.1f | %s | indels_%s" % (
                        sname, thresh, "dedup" if dedup else "all rows",
                        "in" if ind else "out")
                    results[key] = r
                    print(hdr % (sname[:20], "%.1f" % thresh,
                                 "yes" if dedup else "no",
                                 "in" if ind else "out",
                                 r["n_resistant_total"], r["n_truncating"],
                                 "%.1f%%" % r["pct_truncating"],
                                 "%.1f%%" % r["pct_missense"]))

    near = {k: v for k, v in results.items()
            if v["pct_truncating"] is not None and abs(v["pct_truncating"] - 44.4) <= 1.5}
    print("")
    if near:
        print("Specifications reproducing 44.4%% within 1.5 points:")
        for k, v in near.items():
            print("   %s -> %.1f%% (%d/%d)"
                  % (k, v["pct_truncating"], v["n_truncating"], v["n_resistant_total"]))
    else:
        print("NO specification reproduces 44.4%. The figure must be replaced with a")
        print("pre-specified value from this table, or the claim reworded to one that")
        print("the data support. It appears in the abstract, so this blocks submission.")

    payload = {
        "purpose": ("recomputation of the UL23 truncating fraction, because the stated "
                    "44.4% (186/419) traces to no file under logs/"),
        "manuscript_asserted": {"pct": 44.4, "numerator": 186, "denominator": 419,
                                "traces_to_log": False,
                                "appears_in": ["abstract", "introduction",
                                               "results", "discussion"]},
        "cmi_supplement_status": cmi_status,
        "specifications": results,
        "specifications_matching_44_4_within_1_5_points": list(near.keys()),
        "note": ("'other_unclassified' records are reported separately and are never "
                 "folded into the truncating or missense counts; a specification whose "
                 "n_other_unclassified is large should not be used."),
    }
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=1)
    print("")
    print("Wrote %s" % os.path.relpath(OUT, ROOT))


if __name__ == "__main__":
    main()
