# PI Sign-Off Form
## HSV Antiviral Resistance Benchmark Paper 2 — Formal Approval Checklist
**Project:** HSV-1 Antiviral Resistance Mutation Benchmark  
**Protocol version:** AUDIT_PROTOCOL.md Strict Version 4 (Frozen)  
**Prepared by:** Antigravity agent — 2026-10-31  
**Signed by:** [PI NAME], [INSTITUTION]  
**Date of sign-off:** [DATE]

---

> [!IMPORTANT]
> This form must be completed BEFORE any results are submitted to a journal or posted as a preprint.
> All items in Sections 1–5 must have a recorded decision.

---

## Section 1: Data Provenance & Scope

### 1.1 Primary Dataset
- Dataset: `data/processed/primary_position_folds.tsv`
- SHA-256: `7C4D0A3DE2CD171B07230757F44E0DFBD61C6A0A8570DC8F263195B354C66909`
- Rows: 929 (UL23-ACV: 286, UL30-ACV: 332, UL30-FOS: 311)
- Source databases: HerpesDRG-DB, NCBI ClinVar, published literature

| Checkpoint | Status |
|------------|--------|
| PI has reviewed the data provenance document (`docs/data_provenance.md` if exists) | [ ] Confirmed / [ ] Not yet |
| Elsevier-derived rows are excluded from any public release | [ ] Confirmed / [ ] Not applicable |
| Dataset composition is within the scope of [GRANT NUMBER/TITLE] | [ ] Confirmed / [ ] Out of scope → see Note |

**PI Note:**

---

### 1.2 Resistance Thresholds
- Threshold file: `thresholds.tsv`  
- Human threshold citation review: `docs/threshold_citation_checklist.md`

| Checkpoint | Status |
|------------|--------|
| Threshold citation checklist completed by human reviewer | [ ] Yes / [ ] Not yet |
| All CITED thresholds confirmed in primary source | [ ] Yes / [ ] Partially — see Note |
| FLAG_HUMAN_NON_FC items resolved | [ ] Yes / [ ] Partially |
| BLANK_FLAG_HUMAN items resolved (Brivudine, Cidofovir) | [ ] Yes / [ ] Out of scope |
| Any threshold changed → dated protocol amendment filed | [ ] Yes / [ ] N/A |

**PI Note:**

---

### 1.3 Co-Annotator Kappa
- Kappa sample: `logs/kappa_sample_50rows.tsv` (50 rows, stratified 25R/25S)
- Required: Second annotator (PI or designated co-author) independently labels each row

| Checkpoint | Status |
|------------|--------|
| Second annotator has completed `kappa_sample_50rows.tsv` | [ ] Yes / [ ] Not yet |
| Cohen's kappa computed and recorded here: κ = _______ | [ ] Computed / [ ] Not yet |
| κ ≥ 0.70 threshold met | [ ] Yes / [ ] No → see Action |
| If κ < 0.70: disagreement-resolution protocol applied | [ ] Yes / [ ] N/A |

**PI Note:**

**Action if κ < 0.70:**

---

## Section 2: Evaluation Protocol & Statistical Claims

### 2.1 Primary Metric
The primary metric is the **contiguous-block 5-fold mean AUROC**, with a **block-level Student's t 95% confidence interval** computed from per-block AUROC values. NOT pooled AUROC.

| Checkpoint | Status |
|------------|--------|
| PI confirms the primary metric as stated above | [ ] Confirmed |
| PI confirms pooled AUROC is used as secondary descriptive only | [ ] Confirmed |
| PI confirms protocol amendment dated 2026-10-29/30 (withdrawal of condition (a)) | [ ] Reviewed and approved |

**PI Note:**

---

### 2.2 Stopping Rule (amended)
- Original condition (a): "beat position-only RF under position-grouped CV" — **WITHDRAWN** by dated amendment 2026-10-29/30
- Condition (b) retained: block mean ≥ 0.65 with block t-interval above 0.50, on ≥ 2/3 pairs

| Checkpoint | Status |
|------------|--------|
| PI approves withdrawal of condition (a) | [ ] Approved / [ ] Disputed → see Note |
| PI confirms stopping rule result: 3/3 pairs pass for ESM-2 650M | [ ] Confirmed |
| PI confirms stopping rule result: 3/3 pairs pass for ProtBERT | [ ] Confirmed / [ ] Needs re-verification |

**PI Note:**

---

### 2.3 Model Selection
- Co-primary models: **ProtBERT** (Rostlab/prot_bert) and **ESM-2 650M** (facebook/esm2_t33_650M_UR50D)
- Both are zero-shot (no fine-tuning on the 929 rows)
- ESM-2 35M is a secondary baseline only

| Checkpoint | Status |
|------------|--------|
| PI confirms co-primary designation (not post-hoc winner selection) | [ ] Confirmed |
| No fine-tuning performed on primary cohort | [ ] Confirmed |

**PI Note:**

---

## Section 3: Novelty Claims

### 3.1 Novelty Log Status
- File: `logs/novelty_widened_raw.tsv`
- Human decision column: must be completed by PI or designated reviewer

| Checkpoint | Status |
|------------|--------|
| All 6 Track A novelty queries logged with top-10 results | [ ] Yes / [ ] Partial |
| Human decision column filled for all entries | [ ] Yes / [ ] Partial |
| Gate 0 decision recorded (proceed to Phase 2 or discussion-only) | [ ] Yes / [ ] Not yet |

**Gate 0 Decision:**
> [ ] Items 1–3 appear crowded or B1/B2 deflate the dissociation → mechanism story is Discussion only. Phase 2 skipped.  
> [ ] Items 1–3 appear open AND dissociation survives conservation baseline → Phase 2 is justified.

**PI Note:**

---

### 3.2 Novelty Claims Permitted
Based on the novelty log decisions, the following claims are permitted in the paper:

| Claim | Permitted? |
|-------|-----------|
| "Contiguous-block holdout reveals positional shortcut in HSV resistance CV" | [ ] Yes / [ ] Qualified |
| "Zero-shot PLMs outperform tabular baselines under strict structural holdout" | [ ] Yes / [ ] Qualified |
| "PLM signal on UL23 is tied with conservation baseline" | [ ] Yes (descriptive) |
| "Mechanism claim: [specific mechanism]" | [ ] Requires Gate 0 approval |

**PI Note:**

---

## Section 4: Grant Scope Confirmation

| Item | Confirmation |
|------|-------------|
| This work falls within the scope of grant [GRANT NUMBER/TITLE] | [ ] Yes / [ ] Partially |
| No data were used that are outside the grant data sharing agreement | [ ] Confirmed |
| Any Elsevier-source rows are excluded from public datasets | [ ] Confirmed |
| Authorship list agreed | [ ] Yes — list: |

**Grant Number:**  
**Grant Title:**  
**Authorship (in order):**  

**PI Note:**

---

## Section 5: Provisional Status Lift

All results are currently **PROVISIONAL** pending this sign-off. They become reportable only when all checkboxes in Sections 1–4 are filled and this form is signed.

| Final condition | Status |
|----------------|--------|
| κ ≥ 0.70 (or disagreement resolved) | [ ] Met |
| All threshold citations confirmed | [ ] Met |
| Novelty log human decisions complete | [ ] Met |
| Grant scope confirmed | [ ] Met |
| Protocol amendments reviewed and approved | [ ] Met |

---

## Signatures

| Role | Name | Signature | Date |
|------|------|-----------|------|
| Principal Investigator | | | |
| First Author | | | |
| Second Annotator (kappa) | | | |
| Data Custodian (if applicable) | | | |

---

*This form is part of the permanent audit record for this project.*  
*File this signed copy at: `docs/pi_signoff_form_SIGNED_[DATE].pdf`*
