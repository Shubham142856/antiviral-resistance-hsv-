# Dual-Annotator Kappa Curation Instructions

**File to annotate:** `logs/kappa_sample_50rows.tsv`  
**Purpose:** Independent second-annotator review for Cohen's κ computation  
**Required before:** Any journal submission or preprint posting  
**Governing protocol:** `AUDIT_PROTOCOL.md` Strict Version 4

---

## What You Are Being Asked To Do

You will independently label each of the 50 rows in `logs/kappa_sample_50rows.tsv` as either **Resistant** or **Sensitive** based solely on the information provided (variant, drug, mutation category, and cited sources).

**You must NOT look at the `consensus_phenotype` column while making your decision.**  
Fill only the three annotator columns: `annotator_phenotype`, `annotator_confidence`, and `annotator_notes`.

---

## Step-by-Step Instructions

### Step 1 — Open the TSV file

Open `logs/kappa_sample_50rows.tsv` in Excel, LibreOffice Calc, or any TSV-capable editor.

The columns are:

| Column | Description |
|:-------|:------------|
| `tuple_key` | Unique identifier (Virus\|Gene\|Drug\|Mutation) |
| `virus` | HSV-1 or HSV-2 |
| `gene` | UL23 (Thymidine Kinase) or UL30 (DNA Polymerase) |
| `drug` | Aciclovir, Foscarnet, Penciclovir, Brivudine, Pritelivir, Amenamevir, or Cidofovir |
| `aa_change` | Amino acid change (e.g. M121K = Met→Lys at position 121; N306* = stop codon) |
| `mutation_category` | `single_substitution`, `nonsense_stop`, `frameshift`, `indel` |
| `consensus_phenotype` | **DO NOT READ THIS COLUMN until after you finish all 50 rows** |
| `sources` | PubMed IDs or DOIs. You may consult these. |
| `num_sources` | Number of independent sources |
| `annotator_phenotype` | **YOU FILL THIS: Resistant or Sensitive** |
| `annotator_confidence` | **YOU FILL THIS: High / Medium / Low** |
| `annotator_notes` | **YOU FILL THIS: Optional brief note if you are uncertain** |

---

### Step 2 — Labelling Rules

Use these rules in order:

1. **Nonsense/Frameshift (nonsense_stop, frameshift, indel) in UL23:**  
   → Label **Resistant** if drug is Aciclovir, Penciclovir, or Brivudine (requires TK activation).  
   → Label **Sensitive** if drug is Foscarnet or Cidofovir (does NOT require TK activation).

2. **Single substitution — consult source(s) listed:**  
   → Look up the PMID(s) or DOI(s) in PubMed or the publisher site.  
   → Find the phenotypic susceptibility data for this specific variant and drug.  
   → Apply the consensus threshold: **Aciclovir IC₅₀ ≥ 2.0 µg/mL = Resistant**; **Foscarnet IC₅₀ ≥ 100 µg/mL = Resistant**; or use the qualitative label if reported.

3. **If you cannot determine from the source:** Label as the phenotype stated in the source abstract/table, or write `UNCLEAR` in `annotator_notes` and label your best guess in `annotator_phenotype`.

4. **Confidence scale:**
   - `High` — clearly stated in a primary source table with the specific IC₅₀ or qualitative call
   - `Medium` — inferred from context or stated without an exact value
   - `Low` — source ambiguous, conflicting, or not accessible

---

### Step 3 — Save and Return

After completing all 50 rows, save the file as TSV (Tab-separated) with the same name and return it to the primary investigator. Do not change any other columns.

---

### Step 4 — Cohen's Kappa Computation (Primary Investigator)

Once both annotators have completed their columns independently, run:

```bash
cd d:\antiviral-resistance-hsv--1
python scripts/compute_kappa.py
```

This script will:
1. Read `logs/kappa_sample_50rows.tsv` (columns `consensus_phenotype` vs `annotator_phenotype`)
2. Compute Cohen's κ
3. Report agreement table, κ value, and whether κ ≥ 0.70 threshold is met
4. Write result to `logs/kappa_result.json`

**Threshold:** κ ≥ 0.70 is required. If κ < 0.70, resolve disagreements using the protocol in `docs/pi_signoff_form.md` Section 1.3.

---

## Resistance Threshold Reference Card

| Drug | Assay | Resistant threshold |
|:-----|:------|:--------------------|
| Aciclovir | PRA / Dye Uptake | IC₅₀ ≥ 2.0 µg/mL OR fold-change ≥ 3× |
| Foscarnet | PRA / Dye Uptake | IC₅₀ ≥ 100 µg/mL OR fold-change ≥ 3× |
| Penciclovir | PRA | IC₅₀ ≥ 2.0 µg/mL OR fold-change ≥ 3× |
| Brivudine | PRA | Fold-change ≥ 3× (no universal absolute cutoff) |
| Pritelivir / Amenamevir | PRA | Fold-change ≥ 3× |
| Cidofovir | PRA | Fold-change ≥ 3× |

*Source: Safrin et al. AAC 1994, Sauerbrei et al. JAC 2016, Dähne et al. CMI 2025.*

---

*Instructions prepared by: Antigravity agent, 2026-11-22. Review and approve before distributing to co-annotator.*
