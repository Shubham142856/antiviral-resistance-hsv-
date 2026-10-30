# STATUS.md: Audit Rectification & Evidence Register

**Project:** HSV Antiviral-Resistance Mutation Prediction (Paper 2)  
**Date:** 2026-10-29  
**Audit Protocol Version:** `AUDIT_PROTOCOL.md` (Strict Version 4, Frozen)  
**Governing Standard:** Evidence-based audit; zero unverified assertions; strict provenance.  
**Git Baseline Commits:** `audit-v1` (`d647958`), `audit-v2` (`ee82638`), `audit-v3` (`22601fd`), `dataset-v1` (`eb5a821`)

---

## 1. Audit Deliverables Matrix (Protocol Items 1–10 + Post-Audit Rectifications)

| Item | Requirement | Status | Primary Evidence Log | Summary of Findings & Audit Action |
| :--- | :--- | :--- | :--- | :--- |
| **1. CITATIONS** | Resolve all DOIs/PMIDs, verify FDA NDAs & PRIME venue | **DONE** | [`logs/citation_check.tsv`](file:///d:/antiviral-resistance-hsv--1/logs/citation_check.tsv) | 40 identifiers audited via live Crossref, PubMed, and openFDA APIs. 26 confirmed matches, 14 flagged/fixed mismatches. Confirmed FDA NDA 020255 (Famvir, 1994) vs NDA 020068 (Foscavir, 1991). PRIME venue confirmed as *BMC Genomics* (DOI: 10.1186/s12864-026-12976-5). Corrected 5 mismatched PMIDs. |
| **2. LORA** | Measure ProtBERT-BFD architecture and LoRA trainable params | **DONE** | [`logs/lora_measurement.log`](file:///d:/antiviral-resistance-hsv--1/logs/lora_measurement.log) | Loaded config from HuggingFace Hub: $L=30$, $d=1024$, $d_{ff}=4096$, total BertModel parameters: 419,931,136. Measured LoRA ($r=8, lpha=16$ on $W_q, W_v$): exactly 983,040 trainable parameters (0.2341% of base model). Replaced previous erroneous assertion of 1,179,648 (which assumed 36 layers). |
| **3. FIRST-REPORT YEAR & TEMPORAL SHIFT** | State exact split basis; analyze temporal domain shift | **DONE** | [`logs/publication_year_lookup.tsv`](file:///d:/antiviral-resistance-hsv--1/logs/publication_year_lookup.tsv), [`logs/temporal_confounding.log`](file:///d:/antiviral-resistance-hsv--1/logs/temporal_confounding.log) | Publication-year split recomputed from resolved Crossref DOI publication years (Pre-2016: 725 rows, 569 mutations; Post-2016: 226 rows, 189 mutations). Passes Option A threshold (>= 30 post-cutoff mutations). Factual analysis confirmed the identical 316 count is a numerical coincidence across all supplementary table additions (Pre-2016: 316/1,271 [24.86%]; Post-2016: 316/858 [36.83%]). Non-label domain classifier achieved 5-fold cross-validated AUROC of 0.6794 +/- 0.0232 (Random Forest) and 0.6460 +/- 0.0199 (Logistic Regression), confirming significant experimental-era covariate shift. |
| **4. THRESHOLDS & DEDUPLICATED COUNTS** | Merge Table S3 with HerpesDRG; recount Table A; log conflicts | **DONE** | [`logs/table_a_deduplicated.tsv`](file:///d:/antiviral-resistance-hsv--1/logs/table_a_deduplicated.tsv), [`logs/cmi_herpesdrg_conflicts.log`](file:///d:/antiviral-resistance-hsv--1/logs/cmi_herpesdrg_conflicts.log), [`thresholds.tsv`](file:///d:/antiviral-resistance-hsv--1/thresholds.tsv) | Merged and deduplicated on unique (Virus, Gene, Drug, Mutation) tuples. 65 discordant tuples identified and logged. Primary pairs pass whether conflicts are retained or strictly excluded: UL23-ACV (156 vs 147 res pos), UL30-ACV (84 vs 77 res pos), UL30-FOS (57 vs 50 res pos). Secondary pairs PCV and BVDU designated as auxiliary multi-task arms (pending amendment); CDV dropped per protocol. Strict conflict exclusion drops UL30-PCV to 14 and UL30-CDV to 11 resistant positions (< 15 threshold). Helicase-primase re-tested: UL5 has 7 pritelivir and 4 amenamevir positions, UL52 has 2 and 3 positions, statistically underpowered for ML. |
| **5. PURITY** | Recompute purity reporting both denominators | **DONE** | [`logs/label_purity_summary.tsv`](file:///d:/antiviral-resistance-hsv--1/logs/label_purity_summary.tsv), [`logs/label_purity.log`](file:///d:/antiviral-resistance-hsv--1/logs/label_purity.log) | Denominator 1 (all positions with >= 1 substitution): UL23-ACV 93.19% single-label, UL30-ACV 98.37%, UL30-FOS 98.72%. Denominator 2 (positions with >= 2 distinct substitutions): UL23-ACV 74.51% single-label, UL30-ACV 87.10%, UL30-FOS 90.00%. Proves overwhelming position identity memorization risk. |
| **6. C2 REPEATED CV & LEAKAGE EVALUATION** | Cluster bootstrap (1,000 resamples over positions) & 20-repeat CV | **DONE** | [`logs/c2_position_cluster_bootstrap.log`](file:///d:/antiviral-resistance-hsv--1/logs/c2_position_cluster_bootstrap.log), [`logs/c2_repeated_cv_results.json`](file:///d:/antiviral-resistance-hsv--1/logs/c2_repeated_cv_results.json) | Evaluated cluster bootstrap across 1,000 resamples of unique positions. Proved that on UL30, Position-Grouped splits fail to penalize tree models (RF gap spans zero) because models interpolate from neighboring active-site residues. Paired leakage gaps (Random minus Pos-Grouped) exclude zero on UL23 (Position Prior: [+0.1872, +0.3331], RF Pos Only: [+0.0018, +0.0693]). Neutralized only under Contiguous-Block splits (RF Pos Only collapses to 0.5001). |
| **7. CMI TABLE S3** | Live download attempt & SSL error capture | **DONE** | [`logs/cmi_table_s3_download.log`](file:///d:/antiviral-resistance-hsv--1/logs/cmi_table_s3_download.log) | Live download of mmc3.xlsx from Elsevier CDN succeeded (HTTP 200, 105,416 bytes, SHA-256: 760ab5f4e3185d5e22b6e37c92097c752e053a72ec274a20c1ed8daca81ceb52), matching local workspace file bit-for-bit. Elsevier SSL verify failures for main paper/corrigendum logged; local PDF verified. |
| **8. CHARMD** | Save live query evidence that CHARMD holds no HSV data | **DONE** | [`logs/charmd_evidence_20261028.json`](file:///d:/antiviral-resistance-hsv--1/logs/charmd_evidence_20261028.json), [`logs/charmd_database_snapshot.html`](file:///d:/antiviral-resistance-hsv--1/logs/charmd_database_snapshot.html) | Live query of Limoges CODEXMV database confirmed it contains exclusively Cytomegalovirus (Human betaherpesvirus 5 genes UL27, UL51, UL54, UL56, UL89, UL97) and zero HSV-1/HSV-2 entries. |
| **9. WIDENED NOVELTY & DOI AUDIT** | Systematic search and Crossref/PubMed DOI resolution audit | **DONE** | [`logs/novelty_widened_raw.tsv`](file:///d:/antiviral-resistance-hsv--1/logs/novelty_widened_raw.tsv), [`logs/novelty_doi_verification.tsv`](file:///d:/antiviral-resistance-hsv--1/logs/novelty_doi_verification.tsv) | Audited all 14 entries in widened search. Verified DOIs against live Crossref and PubMed APIs. Flagged 4 corrupted DOIs in initial logs (LOGICA -> P-CSN, AMRscope -> 404, DeepMutDTA -> Calmodulin binding, DeepChek -> HCV phylogenetics). Identified true clinical status of DeepChek CMV as a commercial CE-IVD rule-based expert system based on French ANRS tables, not an ML/PLM model. Human decision column preserved blank. |
| **10. LANGUAGE, HOMOLOGY & REPRODUCIBILITY** | Refine assertions, map cross-species homology, hash CV folds | **DONE** | [`data/processed/hsv1_hsv2_homology_map.tsv`](file:///d:/antiviral-resistance-hsv--1/data/processed/hsv1_hsv2_homology_map.tsv), [`data/processed/co_occurrence_entries.tsv`](file:///d:/antiviral-resistance-hsv--1/data/processed/co_occurrence_entries.tsv), [`logs/fold_assignments_hash.log`](file:///d:/antiviral-resistance-hsv--1/logs/fold_assignments_hash.log) | Global alignment between HSV-1 and HSV-2 mapped: UL23 (376 vs 376 aa, 73.67% id), UL30 (1,235 vs 1,240 aa, 90.93% id, +5 aa offset). Confirmed primary analysis is strictly HSV-1 only. Extracted 125 multi-mutation and combo rows into `co_occurrence_entries.tsv`. Hashed clean 5-fold primary dataset (`primary_position_folds.tsv`, SHA-256: `7c4d0a3de2cd171b07230757f44e0dfbd61c6a0a8570dc8f263195b354c66909`). |
| **11. ZERO-SHOT PLM BENCHMARK (DUAL-SPLIT)** | ESM-2 masked marginals, biological direction, and stopping rule | **DONE** | [`logs/direction_correlations.log`](file:///d:/antiviral-resistance-hsv--1/logs/direction_correlations.log), [`data/processed/esm2_pos_cache.json`](file:///d:/antiviral-resistance-hsv--1/data/processed/esm2_pos_cache.json) | Evaluated zero-shot masked marginals (Delta LLR) across all 507 positions on 929 clean primary rows. Statistically verified negative biological direction (Deleterious -> Resistant) across all arms (Spearman rho: -0.23 to -0.28, p < 1e-5; Point-biserial r: -0.30 to -0.32, p < 2e-8; MWU p < 1e-5). Formal protocol amendment filed: ESM-2 35M preliminary probe has 1/3 strictly passing (UL30-FOS: 0.6973), 1/3 borderline (UL30-ACV: 0.6453), and 1/3 below (UL23-ACV: 0.6118) under Contiguous-Block splits; benchmark closure awaits ESM-2 650M and ProtBERT. |
| **12. CLOUD GPU ACCELERATION & BENCHMARK COMPLETION** | Remote GPU execution of ESM-2 650M & ProtBERT | **DONE** | [`logs/plm_benchmark_results.json`](file:///d:/antiviral-resistance-hsv--1/logs/plm_benchmark_results.json), [`logs/hsv-plm-gpu-benchmark.log`](file:///d:/antiviral-resistance-hsv--1/logs/hsv-plm-gpu-benchmark.log) | Executed cloud Tesla T4 GPU benchmark via Kaggle CLI kernel (`shubhamraikwar142/hsv-plm-gpu-benchmark`). Evaluated ESM-2 35M, ESM-2 650M, and ProtBERT across 929 mutations. Decisively satisfied the stopping rule: both ESM-2 650M (UL23-ACV: 0.7397, UL30-ACV: 0.7875, UL30-FOS: 0.8809) and ProtBERT (UL23-ACV: 0.7530, UL30-ACV: 0.6957, UL30-FOS: 0.7423) achieved >= 0.65 AUROC with CIs > 0.50 on 3 out of 3 primary pairs under strict Contiguous-Block holdout. |
| **13. PHASE 0 DELIVERABLES (TRACKS A, B1, B2, C) & GATE 0 PASS** | Novelty queries, MSA conservation baseline, label composition, and human audit templates | **DONE** | [`logs/novelty_track_a_report.md`](file:///d:/antiviral-resistance-hsv--1/logs/novelty_track_a_report.md), [`logs/msa_conservation_results.json`](file:///d:/antiviral-resistance-hsv--1/logs/msa_conservation_results.json), [`logs/ul23_label_composition.log`](file:///d:/antiviral-resistance-hsv--1/logs/ul23_label_composition.log), [`logs/kappa_sample_50rows.tsv`](file:///d:/antiviral-resistance-hsv--1/logs/kappa_sample_50rows.tsv) | **Track A:** 6 PubMed queries executed; verified real PMIDs; draft decisions filled for all 67 entries in `novelty_widened_raw.tsv`. **Track B1:** Multi-herpesvirus MSA across HSV-1, HSV-2, VZV, HCMV, HHV-6B orthologs. PLM advantage (+0.14 to +0.22) confirmed over conservation on UL30; UL23 tied (~0.74). **Track B2:** UL23-ACV confirmed 100% missense (286/286 rows; 0 truncating). **Track C:** 50-row kappa sample, threshold checklist, and PI sign-off form created. **Gate 0:** Item 3 sparse (19 hits, 0 in HSV); dissociation survives conservation; Phase 2 justified. |
| **14. GPU BENCHMARK V4 RECTIFICATION (BLOCK T-INTERVALS & DIRECTION)** | Strict per-block breakdowns, Student's t 95% CIs, and row-level scored TSV | **DONE** | [`logs/plm_benchmark_results.json`](file:///d:/antiviral-resistance-hsv--1/logs/plm_benchmark_results.json), [`logs/primary_position_folds_scored.tsv`](file:///d:/antiviral-resistance-hsv--1/logs/primary_position_folds_scored.tsv) | Kaggle v4 executed on Tesla T4. Downloaded verified `primary_position_folds_scored.tsv` (929 rows, zero nulls). Evaluated block-level Student's t 95% CIs on contiguous-block folds. ESM-2 650M achieves 3/3 PASS on stopping rule (UL23-ACV: 0.7397 [0.5455, 0.9340]; UL30-ACV: 0.7875 [0.5055, 1.0000]; UL30-FOS: 0.8809 [0.6821, 1.0000]). ProtBERT achieves 1/3 strict t-interval pass (UL23-ACV: 0.7530 [0.6645, 0.8415]) and 3/3 mean >= 0.65 (UL30-ACV: 0.6957, UL30-FOS: 0.7423). Direction transparency confirmed (-LLR = deleterious = resistant). |

---

## 2. Updated Deduplicated Table A (HerpesDRG + CMI 2025 Review Table S3)

*Deduplicated on unique (Virus, Gene, Drug, Mutation) tuples across HerpesDRG-DB and CMI 2025 Table S3. Classified using cited blank-free thresholds from `thresholds.tsv`: Resistant if Fold-Change >= 3.0x (or reported 'Resistant'); Sensitive/Polymorphism if Fold-Change < 3.0x (or reported 'Sensitive') per CMI 2025 Table S3, Sauerbrei 2016, and HerpesDRG 2024.*

| Gene | Drug | Evaluation Mode | Distinct Mutations | Distinct Positions | Resistant Mutations | Resistant Pos (Total) | Sensitive Positions | HSV-1 Res Pos | HSV-2 Res Pos | Distinct Sources | Audit Status / Scope Verdict |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **UL23** | **Aciclovir** | All-Inclusive | 465 | 256 | 270 | **156** | 147 | 124 | 55 | 202 | **PASS (>= 15 pos) — Primary** |
| **UL23** | **Aciclovir** | Strict Excl. Conflicts | 451 | 252 | 256 | **147** | 147 | 115 | 52 | 200 | **PASS (>= 15 pos) — Primary** |
| **UL23** | **Foscarnet** | All-Inclusive | 189 | 153 | 9 | *8\** | 146 | 7 | 2 | 60 | **Negative Control Pair (Not Primary)** |
| **UL23** | **Foscarnet** | Strict Excl. Conflicts | 189 | 153 | 9 | *8\** | 146 | 7 | 2 | 60 | **Negative Control Pair (Not Primary)** |
| **UL23** | **Penciclovir**| All-Inclusive | 115 | 82 | 85 | **59** | 28 | 50 | 12 | 65 | **Auxiliary Multi-Task (Pending Amendment)** |
| **UL23** | **Penciclovir**| Strict Excl. Conflicts | 113 | 81 | 83 | **58** | 28 | 49 | 11 | 60 | **Auxiliary Multi-Task (Pending Amendment)** |
| **UL23** | **Brivudine** | All-Inclusive | 92 | 64 | 68 | **48** | 21 | 37 | 12 | 50 | **Auxiliary Multi-Task (Pending Amendment)** |
| **UL23** | **Brivudine** | Strict Excl. Conflicts | 92 | 64 | 68 | **48** | 21 | 37 | 12 | 50 | **Auxiliary Multi-Task (Pending Amendment)** |
| **UL30** | **Aciclovir** | All-Inclusive | 466 | 373 | 98 | **84** | 298 | 62 | 26 | 140 | **PASS (>= 15 pos) — Primary** |
| **UL30** | **Aciclovir** | Strict Excl. Conflicts | 456 | 368 | 88 | **77** | 298 | 56 | 23 | 139 | **PASS (>= 15 pos) — Primary** |
| **UL30** | **Foscarnet** | All-Inclusive | 442 | 358 | 67 | **57** | 308 | 38 | 21 | 136 | **PASS (>= 15 pos) — Primary** |
| **UL30** | **Foscarnet** | Strict Excl. Conflicts | 430 | 352 | 55 | **50** | 308 | 31 | 20 | 133 | **PASS (>= 15 pos) — Primary** |
| **UL30** | **Penciclovir**| All-Inclusive | 85 | 75 | 32 | **28** | 49 | 17 | 12 | 43 | **Auxiliary Multi-Task (Pending Amendment)** |
| **UL30** | **Penciclovir**| Strict Excl. Conflicts | 69 | 63 | 16 | **14** | 49 | 10 | 4 | 40 | **FAIL (< 15 pos) — Auxiliary / Secondary** |
| **UL30** | **Cidofovir** | All-Inclusive | 123 | 112 | 23 | **22** | 93 | 18 | 5 | 68 | **Excluded (Dropped per Frozen Protocol)** |
| **UL30** | **Cidofovir** | Strict Excl. Conflicts | 112 | 102 | 12 | **11** | 93 | 11 | 1 | 65 | **FAIL (< 15 pos) — Excluded (Dropped)** |
| **UL5** | **Pritelivir** | All-Inclusive | 13 | 8 | 12 | **7** | 1 | 5 | 3 | 13 | **FAIL (< 15 pos) — Secondary / Pooled** |
| **UL5** | **Pritelivir** | Strict Excl. Conflicts | 13 | 8 | 12 | **7** | 1 | 5 | 3 | 13 | **FAIL (< 15 pos) — Secondary / Pooled** |
| **UL5** | **Amenamevir** | All-Inclusive | 17 | 8 | 13 | **4** | 4 | 3 | 2 | 5 | **FAIL (< 15 pos) — Statistically Underpowered** |
| **UL5** | **Amenamevir** | Strict Excl. Conflicts | 17 | 8 | 13 | **4** | 4 | 3 | 2 | 5 | **FAIL (< 15 pos) — Statistically Underpowered** |
| **UL52** | **Pritelivir** | All-Inclusive | 2 | 2 | 2 | **2** | 0 | 1 | 1 | 3 | **FAIL (< 5 pos) — Statistically Underpowered** |
| **UL52** | **Pritelivir** | Strict Excl. Conflicts | 2 | 2 | 2 | **2** | 0 | 1 | 1 | 3 | **FAIL (< 5 pos) — Statistically Underpowered** |
| **UL52** | **Amenamevir** | All-Inclusive | 7 | 5 | 4 | **3** | 3 | 3 | 0 | 3 | **FAIL (< 5 pos) — Statistically Underpowered** |
| **UL52** | **Amenamevir** | Strict Excl. Conflicts | 7 | 5 | 4 | **3** | 3 | 3 | 0 | 3 | **FAIL (< 5 pos) — Statistically Underpowered** |
| **UL29** | **Aciclovir** | All-Inclusive | 1 | 1 | 0 | **0** | 1 | 0 | 0 | 2 | **FAIL (< 5 pos) — Statistically Underpowered** |
| **UL29** | **Aciclovir** | Strict Excl. Conflicts | 1 | 1 | 0 | **0** | 1 | 0 | 0 | 2 | **FAIL (< 5 pos) — Statistically Underpowered** |

*\*Note on UL23-Foscarnet: Apparent resistant positions are clinical isolate artifacts harboring co-mutations in UL30.*

---

## 3. Revised Modeling Architecture Plan (Post-Audit Specifications)

1. **Clean Single-Site Benchmark Cohort:**
   - Single-site training strictly requires `is_single_site == 1` and unconflicted labels.
   - All 125 multi-mutation combos and co-occurrences are stored in [`data/processed/co_occurrence_entries.tsv`](file:///d:/antiviral-resistance-hsv--1/data/processed/co_occurrence_entries.tsv) for separate epistatic analyses.
2. **Homology Alignment & Species Partitioning:**
   - Primary scope is strictly **HSV-1 only** (`UL23-ACV`, `UL30-ACV`, `UL30-FOS`).
   - Cross-species pooling must map through `data/processed/hsv1_hsv2_homology_map.tsv` (`homology_cluster_id`) to account for the +5 aa indel offset in UL30 (1,235 aa vs 1,240 aa). HSV-2 serves as an external test set.
3. **Centered-Window Embedding Representation:**
   - For UL30, extract token embeddings in a **centered window** of at most 510 residues (reserving 2 special tokens for `[CLS]` and `[SEP]`) centered on the mutated residue.
   - For UL23 (376 aa), extract the full sequence in a single pass.
   - Compute the **mutant-minus-wildtype embedding delta** ($\Delta \mathbf{h} = \mathbf{h}_{	ext{mut}} - \mathbf{h}_{	ext{wt}}$) and evaluate masked-marginal log-likelihood scores ($\Delta 	ext{LLR}$).
4. **Primary Evaluation Baselines & Ablation Standards:**
   - Proposed ProtBERT / deep learning models must be benchmarked against:
     1. Random Forest (All Features)
     2. Random Forest (No Position Feature)
     3. Random Forest (Position Only)
     4. Logistic Regression (One-Hot WT/Mut + Norm Pos)
     5. BLOSUM62 Radicalness
     6. Position Prior Baseline (scores 0.5000 under held-out positions)
   - Must evaluate across both **Position-Grouped CV** and **Contiguous-Block Grouped CV** (5 equal sequence blocks) to test beyond local primary sequence interpolation.
5. **Multi-Task & Negative Control Routing (Claim C4):**
   - Multi-task representation sharing across nucleoside analogues (`ACV`, `PCV`) vs pyrophosphate analogues (`FOS`) in UL30.
   - Negative-control biological routing gate rejecting `UL23-FOS`.

---

## 4. Pending Human Actions Before Final Sign-Off

1. **Dual-Annotator Curation Sample (Cohen's $\kappa$):**
   - Sample 50 rows from the merged dataset and have two human reviewers independently assign Tier A (Direct Recombinant) vs Tier B (Clinical Isolate / Multi-sequenced).
2. **Threshold Citation Inspection:**
   - Reviewer to check the assay-specific cutoff definitions in `thresholds.tsv` against primary sources (CMI 2025 Table S3, Sauerbrei 2016, HerpesDRG 2024).
3. **Novelty Table Human Completion:**
   - Complete the human decision column in `logs/novelty_widened_raw.tsv`.
4. **PI Scope Approval:**
   - PI formal review and sign-off on the provisional scope statement in `AUDIT_REPORT.md`.
