# STATUS.md: Audit Rectification & Evidence Register

**Project:** HSV Antiviral-Resistance Mutation Prediction (Paper 2)  
**Date:** 2026-10-29  
**Audit Protocol Version:** `AUDIT_PROTOCOL.md` (Strict Version 4, Frozen)  
**Governing Standard:** Evidence-based audit; zero unverified assertions; strict provenance.  
**Git Baseline Commit:** `audit-v1` (`d647958`)

---

## 1. Audit Deliverables Matrix (Protocol Items 1–10 + Post-Audit Extensions)

| Item | Requirement | Status | Primary Evidence Log | Summary of Findings & Audit Action |
| :--- | :--- | :--- | :--- | :--- |
| **1. CITATIONS** | Resolve all DOIs/PMIDs, verify FDA NDAs & PRIME venue | **DONE** | [`logs/citation_check.tsv`](file:///d:/antiviral-resistance-hsv--1/logs/citation_check.tsv) | 40 identifiers audited via live Crossref, PubMed, and openFDA APIs. 26 confirmed matches, 14 flagged/fixed mismatches. Confirmed FDA NDA 020255 (Famvir/famciclovir, 1994) vs NDA 020068 (Foscavir/foscarnet, 1991). Confirmed PRIME publication venue is *BMC Genomics* (DOI: 10.1186/s12864-026-12976-5), not Nature Communications. Corrected 5 mismatched PMIDs (Coen & Schaffer: 6246531, Kawakami: 34228537, Biswas: 19403757, Field & Biswas: 21183396, ProtTrans: 34232869). |
| **2. LORA** | Measure ProtBERT-BFD architecture and LoRA trainable params | **DONE** | [`logs/lora_measurement.log`](file:///d:/antiviral-resistance-hsv--1/logs/lora_measurement.log) | Loaded config from HuggingFace Hub: $L=30$, $d=1024$, $d_{ff}=4096$, total BertModel parameters: 419,931,136. Measured LoRA ($r=8, \alpha=16$ on $W_q, W_v$): exactly 983,040 trainable parameters (0.2341% of base model). Replaced previous erroneous assertion of 1,179,648 (which assumed 36 layers). |
| **3. FIRST-REPORT YEAR** | State exact split basis; log publication years | **DONE** | [`logs/publication_year_lookup.tsv`](file:///d:/antiviral-resistance-hsv--1/logs/publication_year_lookup.tsv), [`logs/publication_year_methodology.log`](file:///d:/antiviral-resistance-hsv--1/logs/publication_year_methodology.log) | Documented that 'created_date' is database entry curation date. Split was recomputed strictly from resolved Crossref DOI publication years ('paper_pub_year') and labeled as 'Publication Year Split' (Pre-2016: 725 rows, 569 mutations; Post-2016: 226 rows, 189 mutations, 132 novel). Passes Option A threshold (>= 30 post-cutoff mutations). |
| **4. THRESHOLDS & CMI MERGE** | Merge Table S3 with HerpesDRG; build cited register; recount Table A | **DONE** | [`logs/table_a_merged_recount.tsv`](file:///d:/antiviral-resistance-hsv--1/logs/table_a_merged_recount.tsv), [`logs/cmi_herpesdrg_conflicts.log`](file:///d:/antiviral-resistance-hsv--1/logs/cmi_herpesdrg_conflicts.log), [`thresholds.tsv`](file:///d:/antiviral-resistance-hsv--1/thresholds.tsv) | Merged all 1,600 CMI 2025 Table S3 rows with 2,112 active HerpesDRG rows across HSV-1 and HSV-2. 1,024 exact overlapping tuples evaluated: 959 concordant (93.65%), 65 discordant/conflicting (6.35%). Recount confirms: UL23-ACV (160 resistant positions), UL30-ACV (84 resistant positions), UL30-FOS (57 resistant positions) pass Primary threshold (>= 15 pos). Re-testing Helicase-Primase confirms UL5 has 7 pritelivir and 4 amenamevir resistant positions, and UL52 has 2 and 3 resistant positions, strictly validating that UL5/UL52 must remain in secondary/pooled exploratory analysis. |
| **5. PURITY** | Recompute purity reporting both denominators | **DONE** | [`logs/label_purity_summary.tsv`](file:///d:/antiviral-resistance-hsv--1/logs/label_purity_summary.tsv), [`logs/label_purity.log`](file:///d:/antiviral-resistance-hsv--1/logs/label_purity.log) | Denominator 1 (all positions with >= 1 substitution): UL23-ACV 93.19% single-label, UL30-ACV 98.37%, UL30-FOS 98.72%. Denominator 2 (positions with >= 2 distinct substitutions): UL23-ACV 74.51% single-label, UL30-ACV 87.10%, UL30-FOS 90.00%. Proves overwhelming position identity memorization risk. |
| **6. C2 PROBES & PER-FOLD AUROC** | Benchmark models under 4 split regimes; compute per-fold AUROC | **DONE** | [`logs/c2_split_evaluation_per_fold.json`](file:///d:/antiviral-resistance-hsv--1/logs/c2_split_evaluation_per_fold.json), [`logs/c2_diagnostic_per_fold.log`](file:///d:/antiviral-resistance-hsv--1/logs/c2_diagnostic_per_fold.log) | Recomputed AUROC strictly per-fold and averaged across 5 folds. Proved mathematically that Position Prior under Position-Grouped splitting achieves Mean Per-Fold AUROC = 0.5000 +/- 0.0000 across all folds (every test sample is tied at constant fallback prior). Clarified that pooled AUROC was 0.4123 due to cross-fold prior baseline shifts. Explicitly documented that 'Physicochemical Random Forest' is a 12-feature tabular baseline, not a protein language model. |
| **7. CMI TABLE S3** | Live download attempt & SSL error capture | **DONE** | [`logs/cmi_table_s3_download.log`](file:///d:/antiviral-resistance-hsv--1/logs/cmi_table_s3_download.log) | Live download of mmc3.xlsx from Elsevier CDN succeeded (HTTP 200, 105,416 bytes, SHA-256: 760ab5f4e3185d5e22b6e37c92097c752e053a72ec274a20c1ed8daca81ceb52), matching local workspace file bit-for-bit. Elsevier SSL verify failures for main paper/corrigendum logged; local PDF verified (SHA-256: 621a93c7d51995eb00aed905635b25b575eaf57c12ef2474504305cce5978b90). Commercial supplements excluded from public git releases. |
| **8. CHARMD** | Save live query evidence that CHARMD holds no HSV data | **DONE** | [`logs/charmd_evidence_20261028.json`](file:///d:/antiviral-resistance-hsv--1/logs/charmd_evidence_20261028.json), [`logs/charmd_database_snapshot.html`](file:///d:/antiviral-resistance-hsv--1/logs/charmd_database_snapshot.html) | Live query of Limoges CODEXMV database confirmed it contains exclusively Cytomegalovirus (Human betaherpesvirus 5 genes UL27, UL51, UL54, UL56, UL89, UL97) and zero HSV-1/HSV-2 entries. |
| **9. WIDENED NOVELTY** | Systematic search across C4-C6, cross-viruses, and registries | **DONE** | [`logs/novelty_widened_raw.tsv`](file:///d:/antiviral-resistance-hsv--1/logs/novelty_widened_raw.tsv), [`logs/novelty_widened_execution.log`](file:///d:/antiviral-resistance-hsv--1/logs/novelty_widened_execution.log) | Widened search covering 18 queries across Europe PMC, PubMed, GitHub, Hugging Face Hub (Models & Datasets), and Zenodo. Screened HCMV (UL54/UL97) and VZV (ORF28/ORF36) literature. No existing sequence-based PLM or deep-learning model exists for HSV antiviral resistance prediction. Decision column preserved blank for human reviewer. |
| **10. LANGUAGE & SCOPE** | Replace unproven assertions, fix metric discrepancy, update status | **DONE** | [`AUDIT_REPORT.md`](file:///d:/antiviral-resistance-hsv--1/AUDIT_REPORT.md) | Replaced 'exactly 0 exist' with 'none found in the logged search'; marked C9 as 'pending'; reconciled 0.88 RF vs 0.8418 prior metric; removed 'locked-scope' and 'sealed' assertions; marked approval as 'Pending PI and Second Annotator Sign-Off'. |

---

## 2. Updated Merged Table A (HerpesDRG + CMI 2025 Review Table S3)

*Cutoff rule: Resistant if Fold-Change >= 3.0x or reported Resistant; Sensitive if Fold-Change < 3.0x or reported Sensitive.*

| Gene | Drug | Total Rows | HerpesDRG Rows | CMI 2025 Rows | Distinct Positions | Resistant Positions (Total) | Sensitive Positions | HSV-1 Res Positions | HSV-2 Res Positions | Audit Status / Verdict |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **UL23** | **Aciclovir** | 889 | 468 | 421 | 260 | **160** | 156 | 128 | 55 | **PASS (>= 15 pos) — Primary** |
| **UL23** | **Foscarnet** | 219 | 219 | 0 | 153 | *8\** | 149 | 7 | 2 | **Negative Control Pair (Not Primary)** |
| **UL23** | **Penciclovir**| 189 | 103 | 86 | 82 | **59** | 31 | 50 | 12 | **PASS (>= 15 pos) — Primary Eligible** |
| **UL23** | **Brivudine** | 135 | 77 | 58 | 64 | **48** | 24 | 37 | 12 | **PASS (>= 15 pos) — Primary Eligible** |
| **UL30** | **Aciclovir** | 865 | 432 | 433 | 373 | **84** | 314 | 62 | 26 | **PASS (>= 15 pos) — Primary** |
| **UL30** | **Foscarnet** | 740 | 402 | 338 | 358 | **57** | 325 | 38 | 21 | **PASS (>= 15 pos) — Primary** |
| **UL30** | **Penciclovir**| 140 | 64 | 76 | 75 | **28** | 62 | 17 | 12 | **PASS (>= 15 pos) — Primary Eligible** |
| **UL30** | **Cidofovir** | 224 | 113 | 111 | 112 | **22** | 105 | 18 | 5 | **PASS (>= 15 pos) — Primary Eligible** |
| **UL5** | **Pritelivir** | 24 | 20 | 4 | 8 | **7** | 1 | 5 | 3 | **FAIL (< 15 pos) — Secondary / Pooled** |
| **UL5** | **Amenamevir** | 25 | 14 | 11 | 8 | **4** | 4 | 3 | 2 | **FAIL (< 15 pos) — Secondary / Pooled** |
| **UL52** | **Pritelivir** | 4 | 3 | 1 | 2 | **2** | 0 | 1 | 1 | **FAIL (< 5 pos) — Relegated** |
| **UL52** | **Amenamevir** | 10 | 6 | 4 | 5 | **3** | 3 | 3 | 0 | **FAIL (< 5 pos) — Relegated** |
| **UL29** | **Aciclovir** | 1 | 1 | 0 | 1 | **0** | 1 | 0 | 0 | **FAIL (< 5 pos) — Relegated** |

*\*Note on UL23-Foscarnet: Apparent resistant positions are clinical isolate artifacts harboring co-mutations in UL30.*

---

## 3. Revised Modeling Architecture Plan (Addressing Review Feedback)

1. **Embedding Design (Preventing Signal Dilution):**
   - Discontinue uniform 512-token sliding window aggregation for UL30.
   - For any target mutation, extract token embeddings in a **centered window** of at most 510 residues (reserving 2 special tokens for `[CLS]` and `[SEP]`) centered on the mutated residue.
   - For UL23 (376 aa), extract the full sequence in a single pass.
   - Compute the **mutant-minus-wildtype embedding delta** ($\Delta \mathbf{h} = \mathbf{h}_{\text{mut}} - \mathbf{h}_{\text{wt}}$) and evaluate masked-marginal log-likelihood scores ($\Delta \text{LLR}$) to preserve local biochemical disruption signal.
2. **Claim C4 Revision (Multi-Task & Negative Control Routing):**
   - With primary evaluation focused on `UL23-ACV`, `UL30-ACV`, and `UL30-FOS`, biological pathway routing is redefined around two rigorous empirical tests:
     - **Multi-task drug conditioning:** Sharing representation across nucleoside analogues (`ACV`, `PCV`) vs pyrophosphate analogues (`FOS`) in UL30.
     - **Negative-control pathway gating:** Explicit rejection of `UL23-FOS` (TK does not metabolize foscarnet; predicted probability of TK-mediated resistance must collapse to zero).
3. **Realistic Performance Expectations under Position Holdouts:**
   - With 74% to 90% single-label positions among multi-substitution sites, held-out-position performance is fundamentally a position-level resistance propensity prediction task.
   - Baseline models without structural features will achieve modest AUROCs (~0.60–0.70). Incorporation of evolutionary conservation (ConSurf / Shannon entropy) and structure-derived solvent accessibility (DSSP from AlphaFold structures) is planned for the primary model.

---

## 4. Pending Human Actions Before Final Sign-Off

1. **Dual-Annotator Curation Sample (Cohen's $\kappa$):**
   - Sample 50 rows from the merged dataset and have two human reviewers independently assign Tier A (Direct Recombinant) vs Tier B (Clinical Isolate / Multi-sequenced).
2. **Threshold Citation Inspection:**
   - Reviewer to check the assay-specific cutoff definitions in `thresholds.tsv` against primary sources (CMI 2025 Table S3, Sauerbrei 2016, HerpesDRG 2024).
3. **PI Scope Approval:**
   - PI formal review and sign-off on the provisional scope statement in `AUDIT_REPORT.md`.
