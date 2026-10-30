# Threshold Citation Checklist
## HSV Antiviral Resistance Benchmark — Track C Human Audit
**Purpose:** Human reviewer reads each cited source and confirms the fold-change / IC50 threshold is actually stated in the paper as written.  
**Governing protocol:** AUDIT_PROTOCOL.md Strict Version 4 (Frozen)  
**Prepared by:** Antigravity agent — 2026-10-31  
**Completed by:** [REVIEWER NAME] — [DATE]

---

> [!CAUTION]
> Do NOT accept a threshold from an abstract, secondary summary, or database entry alone.
> You must read the primary source table or methods and confirm the exact number.

---

## Instructions

For each row below:
1. Retrieve the paper via the PMID or DOI.
2. Find the specific table, supplementary table, or methods section that states the threshold.
3. Record the **exact text** from the paper that supports the threshold.
4. Mark **STATUS** as one of: `CONFIRMED`, `WRONG_NUMBER`, `NOT_IN_PAPER`, `ONLY_IN_SUPPLEMENT`, `NEEDS_REPLACEMENT`.
5. If `WRONG_NUMBER` or `NOT_IN_PAPER`: write the correct value or the correct citation in the Notes column.

---

## Table 1: Aciclovir — Plaque Reduction Assay (PRA)

| Field | Value |
|-------|-------|
| **Threshold in `thresholds.tsv`** | Fold-change ≥ 3.0; IC50 > 2.0 µg/mL (> 8.8 µM) |
| **Primary citation** | Sauerbrei et al. JAC 2016 **PMID: 26433780** |
| **Secondary citation** | Charles et al. 2024 (HerpesDRG) |
| **Tertiary citation** | Döhne et al. CMI 2025 Table S3 |
| **Current status in TSV** | CITED |

**Reviewer task:**
- [ ] Read Sauerbrei et al. 2016 (PMID 26433780). Confirm 3.0-fold or 2.0 µg/mL threshold appears explicitly.
- [ ] Read CMI 2025 Review Table S3. Confirm same threshold.
- [ ] Verify PMID 26433780 resolves to the correct paper (not a different Sauerbrei paper).

| Reviewer field | Entry |
|----------------|-------|
| Exact text from paper | |
| Page/Table/Supplementary ref | |
| STATUS | |
| Notes | |

---

## Table 2: Aciclovir — Dye Uptake Assay

| Field | Value |
|-------|-------|
| **Threshold** | Fold-change ≥ 3.0; IC50 > 2.0 µg/mL |
| **Primary citation** | Burrel et al. Antiviral Res 2010 **PMID: 20380849** |
| **Secondary citation** | Döhne et al. CMI 2025 Table S3 |
| **Current status** | CITED |

**Reviewer task:**
- [ ] Read Burrel et al. 2010 (PMID 20380849). Confirm threshold for dye uptake specifically.
- [ ] Note: This assay is methodologically different from PRA — confirm the paper distinguishes them.

| Reviewer field | Entry |
|----------------|-------|
| Exact text from paper | |
| Page/Table/Supplementary ref | |
| STATUS | |
| Notes | |

---

## Table 3: Aciclovir — TK Activity Assay

| Field | Value |
|-------|-------|
| **Threshold** | TK activity < 15% of wild-type |
| **Primary citation** | Sauerbrei et al. JAC 2016 **PMID: 26433780** |
| **Current status** | FLAG_HUMAN_NON_FC (no fold-change value — human decision required) |

> [!WARNING]
> This threshold is NOT a fold-change metric. It measures enzyme activity, not antiviral IC50.
> Decision required: Is this threshold comparable to PRA thresholds in the primary cohort?
> If TK-deficient mutants were labeled using this criterion, they may need a separate flag.

**Reviewer task:**
- [ ] Read Sauerbrei et al. 2016 (PMID 26433780). Confirm 15% WT activity threshold.
- [ ] Decide: Are any rows in the primary cohort labeled using TK activity alone? (Check `sources` column.)
- [ ] If yes: Flag those rows for sensitivity analysis.

| Reviewer field | Entry |
|----------------|-------|
| Exact text from paper | |
| TK-activity-only rows in primary cohort (count) | |
| Decision | INCLUDE_AS_IS / SENSITIVITY_FLAG / EXCLUDE |
| Notes | |

---

## Table 4: Foscarnet — PRA

| Field | Value |
|-------|-------|
| **Threshold** | Fold-change ≥ 3.0; IC50 > 100 µM (> 30 µg/mL) |
| **Primary citations** | Döhne et al. CMI 2025; Sauerbrei et al. JAC 2016 (PMID: 26433780); HerpesDRG (PMID: 39192205) |
| **Current status** | CITED |

**Reviewer task:**
- [ ] Confirm PMID 39192205 resolves to the HerpesDRG paper (Charles et al. 2024).
- [ ] Confirm 100 µM (not 100 µg/mL) is the correct molar threshold (FOS MW ≈ 126 g/mol → 30 µg/mL ≈ 238 µM, so check conversion).
- [ ] Confirm Sauerbrei et al. 2016 states this threshold for Foscarnet specifically.

| Reviewer field | Entry |
|----------------|-------|
| Exact text from paper | |
| Unit conversion verified? (µM vs µg/mL) | |
| STATUS | |
| Notes | |

---

## Table 5: Foscarnet — Dye Uptake

| Field | Value |
|-------|-------|
| **Threshold** | Fold-change ≥ 3.0; IC50 > 100 µM |
| **Primary citations** | Frobert et al. JAC 2014; Döhne et al. CMI 2025 Table S3 |
| **Current status** | CITED |
| **Note** | Frobert et al. JAC 2014 DOI not recorded in TSV — needs retrieval |

**Reviewer task:**
- [ ] Find Frobert et al. JAC 2014 (search: Frobert foscarnet dye uptake 2014 Journal of Antimicrobial Chemotherapy).
- [ ] Record DOI.
- [ ] Confirm 100 µM threshold appears for dye uptake specifically.

| Reviewer field | Entry |
|----------------|-------|
| DOI of Frobert et al. 2014 | |
| Exact text from paper | |
| STATUS | |
| Notes | |

---

## Table 6: Pritelivir — PRA

| Field | Value |
|-------|-------|
| **Threshold** | Fold-change ≥ 3.0; IC50 > 0.1 µM |
| **Primary citations** | Birkmann et al. AAC 2022 **PMID: 35862737**; Schalkwijk et al. JID 2024 |
| **Current status** | CITED |

**Reviewer task:**
- [ ] Confirm PMID 35862737 = Birkmann et al. AAC 2022. Read methods/table for 0.1 µM and 3-fold cutoffs.
- [ ] Find Schalkwijk et al. JID 2024 (no PMID recorded). Record PMID/DOI.

| Reviewer field | Entry |
|----------------|-------|
| Birkmann et al. 2022 threshold text | |
| Schalkwijk et al. JID 2024 PMID/DOI | |
| STATUS | |
| Notes | |

---

## Table 7: Pritelivir — CPE Assay

| Field | Value |
|-------|-------|
| **Threshold** | Fold-change ≥ 3.0 |
| **Primary citations** | Birkmann et al. AAC 2022; Döhne et al. CMI 2025 Table S3 |
| **Current status** | CITED |

**Reviewer task:**
- [ ] Confirm Birkmann et al. 2022 explicitly covers CPE in addition to PRA.

| Reviewer field | Entry |
|----------------|-------|
| CPE threshold text | |
| STATUS | |
| Notes | |

---

## Table 8: Amenamevir — PRA

| Field | Value |
|-------|-------|
| **Threshold** | Fold-change ≥ 3.0 (or > 4.0 in Vero cells) |
| **Primary citations** | Kawakami et al. AAC 2021 **PMID: 34228537**; Döhne et al. CMI 2025 Table S3 |
| **Current status** | CITED |

> [!NOTE]
> Two thresholds listed (3.0 and 4.0 in Vero). Decision: which applies to rows in our dataset?
> Check `sources` to see if any rows use Vero-cell assays with this drug.

**Reviewer task:**
- [ ] Read Kawakami et al. 2021 (PMID 34228537). Confirm 3.0-fold and 4.0-fold (Vero) thresholds.
- [ ] Check whether any primary cohort rows use amenamevir+Vero-cell data.

| Reviewer field | Entry |
|----------------|-------|
| Exact threshold text | |
| Vero-cell rows in primary cohort? | |
| STATUS | |
| Notes | |

---

## Table 9: Amenamevir — CPE

| Field | Value |
|-------|-------|
| **Threshold** | Fold-change ≥ 3.0 |
| **Primary citation** | Kawakami et al. AAC 2021 (PMID: 34228537) |
| **Current status** | CITED |

| Reviewer field | Entry |
|----------------|-------|
| STATUS | |
| Notes | |

---

## Table 10: Brivudine — Any Assay

| Field | Value |
|-------|-------|
| **Threshold** | NOT RECORDED |
| **Citation** | NOT RECORDED |
| **Current status** | BLANK_FLAG_HUMAN |

> [!CAUTION]
> No threshold or citation for Brivudine. Any Brivudine-labeled rows in the dataset
> cannot be verified. Decision required before these rows can be used in any analysis.

**Reviewer task:**
- [ ] Determine if any primary cohort rows use Brivudine phenotype data.
- [ ] If yes: Find an appropriate resistance threshold citation (e.g., CMI 2025 Review, Döhne et al.).
- [ ] If no primary cohort rows: mark as OUT_OF_SCOPE.

| Reviewer field | Entry |
|----------------|-------|
| Brivudine rows in primary cohort? | |
| Proposed threshold + citation | |
| Decision | ADD_THRESHOLD / OUT_OF_SCOPE |
| Notes | |

---

## Table 11: Cidofovir — Any Assay

| Field | Value |
|-------|-------|
| **Threshold** | NOT RECORDED |
| **Citation** | NOT RECORDED |
| **Current status** | BLANK_FLAG_HUMAN |

> [!CAUTION]
> Same issue as Brivudine. No threshold or citation.

**Reviewer task:**
- [ ] Determine if any primary cohort rows use Cidofovir phenotype data.
- [ ] If yes: Find threshold citation (suggested: CMI 2025, Sauerbrei 2016).
- [ ] If no primary cohort rows: mark as OUT_OF_SCOPE.

| Reviewer field | Entry |
|----------------|-------|
| Cidofovir rows in primary cohort? | |
| Proposed threshold + citation | |
| Decision | ADD_THRESHOLD / OUT_OF_SCOPE |
| Notes | |

---

## Summary Sign-Off

| Item | Status |
|------|--------|
| All CITED thresholds confirmed | [ ] Yes / [ ] Partially / [ ] No |
| All FLAG_HUMAN_NON_FC items decided | [ ] Yes / [ ] Partially |
| All BLANK_FLAG_HUMAN items decided | [ ] Yes / [ ] Partially |
| Any threshold changed → protocol amendment filed? | [ ] Yes / [ ] N/A |
| **Reviewer signature** | |
| **Date completed** | |
