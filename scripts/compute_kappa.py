"""Compute Cohen's kappa between consensus_phenotype and annotator_phenotype columns
in logs/kappa_sample_50rows.tsv. Writes logs/kappa_result.json.

Usage:
    python scripts/compute_kappa.py
"""

import json
import pathlib
import csv
import sys
from collections import defaultdict

ROOT = pathlib.Path(__file__).resolve().parent.parent
TSV  = ROOT / "logs" / "kappa_sample_50rows.tsv"
OUT  = ROOT / "logs" / "kappa_result.json"


def cohen_kappa(labels_a, labels_b):
    """Compute Cohen's kappa for two annotators over matched label lists."""
    assert len(labels_a) == len(labels_b), "Label lists must be the same length."
    classes = sorted(set(labels_a) | set(labels_b))
    n = len(labels_a)

    # Observed agreement
    p_o = sum(1 for a, b in zip(labels_a, labels_b) if a == b) / n

    # Expected agreement
    freq_a = defaultdict(int)
    freq_b = defaultdict(int)
    for a, b in zip(labels_a, labels_b):
        freq_a[a] += 1
        freq_b[b] += 1

    p_e = sum((freq_a[c] / n) * (freq_b[c] / n) for c in classes)

    if p_e == 1.0:
        kappa = 1.0  # Perfect agreement, degenerate case
    else:
        kappa = (p_o - p_e) / (1.0 - p_e)

    # Confusion matrix
    confusion = defaultdict(lambda: defaultdict(int))
    for a, b in zip(labels_a, labels_b):
        confusion[a][b] += 1

    return kappa, p_o, p_e, dict(confusion), classes


def main():
    if not TSV.exists():
        print(f"ERROR: {TSV} not found.", file=sys.stderr)
        sys.exit(1)

    rows = []
    with open(TSV, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)

    print(f"Loaded {len(rows)} rows from {TSV.name}")

    # Validate that annotator_phenotype is filled
    missing = [i + 2 for i, r in enumerate(rows)
               if r.get("annotator_phenotype", "").strip() == ""]
    if missing:
        print(f"\nWARNING: {len(missing)} rows have empty annotator_phenotype:")
        print(f"  Row numbers (1-indexed, including header): {missing}")
        print("  Complete annotation before computing kappa.")

    # Only compute on rows where both columns are filled
    paired = [(r["consensus_phenotype"].strip(), r["annotator_phenotype"].strip())
              for r in rows
              if r.get("consensus_phenotype", "").strip()
              and r.get("annotator_phenotype", "").strip()]

    if len(paired) < 2:
        print(f"\nERROR: Insufficient paired labels ({len(paired)}). "
              "Second annotator must complete the TSV first.")
        sys.exit(1)

    labels_consensus = [p[0] for p in paired]
    labels_annotator  = [p[1] for p in paired]

    kappa, p_o, p_e, confusion, classes = cohen_kappa(labels_consensus, labels_annotator)

    threshold_met = kappa >= 0.70
    verdict = "PASSED" if threshold_met else "FAILED — resolve disagreements before submission"

    print(f"\n{'='*55}")
    print(f"  Dual-Annotator Kappa Report")
    print(f"{'='*55}")
    print(f"  Rows evaluated:       {len(paired)} / {len(rows)}")
    print(f"  Observed agreement:   {p_o:.4f}")
    print(f"  Expected agreement:   {p_e:.4f}")
    print(f"  Cohen's kappa (κ):    {kappa:.4f}")
    print(f"  Threshold (κ ≥ 0.70): {verdict}")
    print(f"\n  Confusion matrix (consensus_phenotype rows × annotator_phenotype cols):")
    header = [""] + classes
    print("  " + "\t".join(header))
    for c_row in classes:
        line = [c_row] + [str(confusion.get(c_row, {}).get(c_col, 0)) for c_col in classes]
        print("  " + "\t".join(line))
    print(f"{'='*55}\n")

    result = {
        "tsv_path": str(TSV),
        "rows_total": len(rows),
        "rows_evaluated": len(paired),
        "rows_missing_annotator": len(missing),
        "observed_agreement": round(p_o, 6),
        "expected_agreement": round(p_e, 6),
        "cohen_kappa": round(kappa, 6),
        "threshold": 0.70,
        "threshold_met": threshold_met,
        "verdict": verdict,
        "classes": classes,
        "confusion_matrix": {k: dict(v) for k, v in confusion.items()},
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
