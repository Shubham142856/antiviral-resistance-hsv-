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
| **Threshold in `thresholds.tsv`** | Fold-change ≥ 3.0; IC50 ≥ 2.0 µg/mL (predictive value for treatment failure 95%) |
| **Primary clinical citation** | Safrin et al. *Antimicrob Agents Chemother* 1994;38(6):1246–1250 **PMID: 8092821** [VERIFIED 2026-10-31] |
| **Database compilations** | Sauerbrei et al. JAC 2016 **PMID: 26433780**; Charles et al. 2024 (HerpesDRG); Döhne et al. CMI 2025 Table S3 |
| **Current status in TSV** | VERIFIED |

**Verification details (2026-10-31):**
- Primary text from Safrin et al. 1994: *"The predictive values of a susceptible result (50% effective concentrations, < 2 micrograms/ml for acyclovir and < 100 micrograms/ml for foscarnet) for complete healing of lesions were 62% for acyclovir and 82% for foscarnet; the predictive values of a resistant result for failure to heal were 95% for acyclovir and 88% for foscarnet."*
- Sauerbrei et al. 2016 is an authoritative database review compiling mutations from phenotyped strains, but did not independently define the 2.0 µg/mL cutoff.

---

## Table 2: Aciclovir — Dye Uptake Assay

| Field | Value |
|-------|-------|
| **Threshold** | Fold-change ≥ 3.0; IC50 > 2.0 µg/mL |
| **Primary citation** | Burrel et al. *Antimicrob Agents Chemother* 2010;54(11):4833–4842 **PMID: 20733037** [VERIFIED 2026-10-31] |
| **Secondary citation** | Döhne et al. CMI 2025 Table S3 |
| **Current status** | VERIFIED |

**Correction note (2026-10-31):**
- Previous entry incorrectly listed "Antiviral Res 2010 (PMID: 20380849)". PMID 20380849 is an unrelated nanoparticle paper (*Development of novel fibrinogen nanoparticles by two-step co-acervation method*).
- The true publication is Burrel et al., *Antimicrobial Agents and Chemotherapy* 2010;54(11):4833–4842 (PMID: 20733037, PMCID: PMC2976137).

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

---

## Table 4: Foscarnet — PRA

| Field | Value |
|-------|-------|
| **Threshold** | Fold-change ≥ 3.0; IC50 ≥ 100 µg/mL (predictive value for treatment failure 88%; no µM conversion asserted) |
| **Primary clinical citation** | Safrin et al. *Antimicrob Agents Chemother* 1994;38(6):1246–1250 **PMID: 8092821** [VERIFIED 2026-10-31] |
| **Reaffirmation** | Wang, *J Clin Microbiol* 2024;62(8):e00728-23 **PMID: 39051778** [VERIFIED 2026-10-31] |
| **Database compilations** | Sauerbrei et al. JAC 2016 (PMID: 26433780); Charles et al. 2024 (HerpesDRG); Döhne et al. CMI 2025 |
| **Current status** | VERIFIED |

**Verification details (2026-10-31):**
- Resolved ref 26 in Wang 2024: It is Safrin et al. 1994 (AAC 38:1246-1250; PMID 8092821), which established IC50 ≥ 100 µg/mL based on an 88% predictive value for lesion failure to heal in 115 HIV-infected patients.
- The chemical form was foscarnet sodium; no molar conversion was reported in the primary trial, so reporting strictly in µg/mL is authoritative.
- The phantom "Frobert et al. JAC 2014" and the spurious back-derived "30 µg/mL" conversion are permanently withdrawn.
- [ ] Decide whether the benchmark should cite any absolute cut-off at all, given that labels are assigned on fold-change ≥ 3.0.

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
| **Threshold** | **WITHDRAWN** — previously "Fold-change ≥ 3.0; IC50 > 100 µM" |
| **Primary citations** | "Frobert et al. JAC 2014" — **could not be located; appears not to exist** |
| **Current status** | PERMANENTLY_WITHDRAWN |

> [!CAUTION]
> **Agent pre-check, 2026-10-31 — this citation does not resolve.** A PubMed search for Frobert + foscarnet returns seven papers; **none is from 2014 and none is in *J Antimicrob Chemother***. The Frobert papers that exist are in *Antiviral Research* (2023, 2019), *Viruses* (2025), *Am J Ophthalmol* (2019) and French-language journals. No DOI for "Frobert et al. JAC 2014" was ever recorded in `thresholds.tsv`, and `logs/citation_check.tsv` already flagged that a DOI previously attributed to Frobert (10.1128/AAC.00977-15) is in fact Schmidt et al. AAC 2015.
> This was the **sole cited source for the 100 µM foscarnet cut-off**, which appears in the manuscript Methods §2.1. Until a real source is produced, the 100 µM figure has no provenance.

**Reviewer task:**
- [ ] Determine whether a real source for a 100 µM foscarnet cut-off exists. Candidate to check: Frobert et al., *Antiviral Research* 2014, "Resistance of herpes simplex viruses to acyclovir: an update from a ten-year survey in France" — **note this is an aciclovir paper and may state no foscarnet cut-off at all** (agent did not retrieve it).
- [ ] If no source is found: confirm removal of the 100 µM figure from `thresholds.tsv` and from the manuscript.
- [ ] Confirm whether any primary-cohort rows were labelled using a foscarnet dye-uptake absolute cut-off rather than fold-change. If none were, mark OUT_OF_SCOPE.

| Reviewer field | Entry |
|----------------|-------|
| Real source found? (cite or "none") | |
| Exact text from paper | |
| Rows affected in primary cohort | |
| Decision | REPLACE_CITATION / REMOVE_THRESHOLD / OUT_OF_SCOPE |
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

