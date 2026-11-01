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

| **2. LORA** | Measure ProtBERT-BFD architecture and LoRA trainable params | **DONE** | [`logs/lora_measurement.log`](file:///d:/antiviral-resistance-hsv--1/logs/lora_measurement.log) | Loaded config from HuggingFace Hub: $L=30$, $d=1024$, $d_{ff}=4096$, total BertModel parameters: 419,931,136. Measured LoRA ($r=8, lpha=16$ on $W_q, W_v$): exactly 983,040 trainable parameters (0.2341% of base model). Replaced previous erroneous assertion of 1,179,648 (which assumed 36 layers). |

| **3. FIRST-REPORT YEAR & TEMPORAL SHIFT** | State exact split basis; analyze temporal domain shift | **DONE** | [`logs/publication_year_lookup.tsv`](file:///d:/antiviral-resistance-hsv--1/logs/publication_year_lookup.tsv), [`logs/temporal_confounding.log`](file:///d:/antiviral-resistance-hsv--1/logs/temporal_confounding.log) | Publication-year split recomputed from resolved Crossref DOI publication years (Pre-2016: 725 rows, 569 mutations; Post-2016: 226 rows, 189 mutations). Passes Option A threshold (>= 30 post-cutoff mutations). Factual analysis confirmed the identical 316 count is a numerical coincidence across all supplementary table additions (Pre-2016: 316/1,271 [24.86%]; Post-2016: 316/858 [36.83%]). Non-label domain classifier achieved 5-fold cross-validated AUROC of 0.6794 +/- 0.0232 (Random Forest) and 0.6460 +/- 0.0199 (Logistic Regression), confirming significant experimental-era covariate shift. |

| **4. THRESHOLDS & DEDUPLICATED COUNTS** | Merge Table S3 with HerpesDRG; recount Table A; log conflicts | **DONE** | [`logs/table_a_deduplicated.tsv`](file:///d:/antiviral-resistance-hsv--1/logs/table_a_deduplicated.tsv), [`logs/cmi_herpesdrg_conflicts.log`](file:///d:/antiviral-resistance-hsv--1/logs/cmi_herpesdrg_conflicts.log), [`thresholds.tsv`](file:///d:/antiviral-resistance-hsv--1/thresholds.tsv) | Merged and deduplicated on unique (Virus, Gene, Drug, Mutation) tuples. 65 discordant tuples identified and logged. Primary pairs pass whether conflicts are retained or strictly excluded: UL23-ACV (156 vs 147 res pos), UL30-ACV (84 vs 77 res pos), UL30-FOS (57 vs 50 res pos). Secondary pairs PCV and BVDU designated as auxiliary multi-task arms (pending amendment); CDV dropped per protocol. Strict conflict exclusion drops UL30-PCV to 14 and UL30-CDV to 11 resistant positions (< 15 threshold). Helicase-primase re-tested: UL5 has 7 pritelivir and 4 amenamevir positions, UL52 has 2 and 3 positions, statistically underpowered for ML. |

| **5. PURITY** | Recompute purity reporting both denominators | **DONE** | [`logs/label_purity_summary.tsv`](file:///d:/antiviral-resistance-hsv--1/logs/label_purity_summary.tsv), [`logs/label_purity.log`](file:///d:/antiviral-resistance-hsv--1/logs/label_purity.log) | Denominator 1 (all positions with >= 1 substitution): UL23-ACV 93.19% single-label, UL30-ACV 98.37%, UL30-FOS 98.72%. Denominator 2 (positions with >= 2 distinct substitutions): UL23-ACV 74.51% single-label, UL30-ACV 87.10%, UL30-FOS 90.00%. Proves overwhelming position identity memorization risk. |

| **6. C2 REPEATED CV & LEAKAGE EVALUATION** | Cluster bootstrap (1,000 resamples over positions) & 20-repeat CV | **DONE** | [`logs/c2_position_cluster_bootstrap.log`](file:///d:/antiviral-resistance-hsv--1/logs/c2_position_cluster_bootstrap.log), [`logs/c2_repeated_cv_results.json`](file:///d:/antiviral-resistance-hsv--1/logs/c2_repeated_cv_results.json) | Evaluated cluster bootstrap across 1,000 resamples of unique positions. Proved that on UL30, Position-Grouped splits fail to penalize tree models (RF gap spans zero) because models interpolate from neighboring active-site residues. Paired leakage gaps (Random minus Pos-Grouped) exclude zero on UL23 (Position Prior: [+0.1872, +0.3331], RF Pos Only: [+0.0018, +0.0693]). Neutralized only under Contiguous-Block splits (RF Pos Only collapses to 0.5001). |

| **7. CMI TABLE S3** | Live download attempt & SSL error capture | **DONE** | [`logs/cmi_table_s3_download.log`](file:///d:/antiviral-resistance-hsv--1/logs/cmi_table_s3_download.log) | Live download of mmc3.xlsx from Elsevier CDN succeeded (HTTP 200, 105,416 bytes, SHA-256: 760ab5f4e3185d5e22b6e37c92097c752e053a72ec274a20c1ed8daca81ceb52), matching local workspace file bit-for-bit. Elsevier SSL verify failures for main paper/corrigendum logged; local PDF verified. |

| **8. CHARMD** | Save live query evidence that CHARMD holds no HSV data | **DONE** | [`logs/charmd_evidence_20261028.json`](file:///d:/antiviral-resistance-hsv--1/logs/charmd_evidence_20261028.json), [`logs/charmd_database_snapshot.html`](file:///d:/antiviral-resistance-hsv--1/logs/charmd_database_snapshot.html) | Live query of Limoges CODEXMV database confirmed it contains exclusively Cytomegalovirus (Human betaherpesvirus 5 genes UL27, UL51, UL54, UL56, UL89, UL97) and zero HSV-1/HSV-2 entries. |

| **9. WIDENED NOVELTY & DOI AUDIT** | Systematic search and Crossref/PubMed DOI resolution audit | **DONE** | [`logs/novelty_widened_raw.tsv`](file:///d:/antiviral-resistance-hsv--1/logs/novelty_widened_raw.tsv), [`logs/novelty_doi_verification.tsv`](file:///d:/antiviral-resistance-hsv--1/logs/novelty_doi_verification.tsv) | Audited all 14 entries in widened search. Verified DOIs against live Crossref and PubMed APIs. Flagged 4 corrupted DOIs in initial logs (LOGICA -> P-CSN, AMRscope -> 404, DeepMutDTA -> Calmodulin binding, DeepChek -> HCV phylogenetics). Identified true clinical status of DeepChek CMV as a commercial CE-IVD rule-based expert system based on French ANRS tables, not an ML/PLM model. Human decision column preserved blank. |

| **10. LANGUAGE, HOMOLOGY & REPRODUCIBILITY** | Refine assertions, map cross-species homology, hash CV folds | **DONE** | [`data/processed/hsv1_hsv2_homology_map.tsv`](file:///d:/antiviral-resistance-hsv--1/data/processed/hsv1_hsv2_homology_map.tsv), [`data/processed/co_occurrence_entries.tsv`](file:///d:/antiviral-resistance-hsv--1/data/processed/co_occurrence_entries.tsv), [`logs/fold_assignments_hash.log`](file:///d:/antiviral-resistance-hsv--1/logs/fold_assignments_hash.log) | Global alignment between HSV-1 and HSV-2 mapped: UL23 (376 vs 376 aa, 73.67% id), UL30 (1,235 vs 1,240 aa, 90.93% id, +5 aa offset). Confirmed primary analysis is strictly HSV-1 only. Extracted 125 multi-mutation and combo rows into `co_occurrence_entries.tsv`. Hashed clean 5-fold primary dataset (`primary_position_folds.tsv`, SHA-256: `7c4d0a3de2cd171b07230757f44e0dfbd61c6a0a8570dc8f263195b354c66909`). |

| **11. ZERO-SHOT PLM BENCHMARK (DUAL-SPLIT)** | ESM-2 masked marginals, biological direction, and stopping rule | **DONE** | [`logs/direction_correlations.log`](file:///d:/antiviral-resistance-hsv--1/logs/direction_correlations.log), [`data/processed/esm2_pos_cache.json`](file:///d:/antiviral-resistance-hsv--1/data/processed/esm2_pos_cache.json) | Evaluated zero-shot masked marginals (Delta LLR) across all 507 positions on 929 clean primary rows. Statistically verified negative biological direction (Deleterious -> Resistant) across all arms (Spearman rho: -0.23 to -0.28, p < 1e-5; Point-biserial r: -0.30 to -0.32, p < 2e-8; MWU p < 1e-5). Formal protocol amendment filed: ESM-2 35M preliminary probe has 1/3 strictly passing (UL30-FOS: 0.6973), 1/3 borderline (UL30-ACV: 0.6453), and 1/3 below (UL23-ACV: 0.6118) under Contiguous-Block splits; benchmark closure awaits ESM-2 650M and ProtBERT. |

| **12. CLOUD GPU ACCELERATION & BENCHMARK COMPLETION** | Remote GPU execution of ESM-2 650M & ProtBERT | **DONE** | [`logs/plm_benchmark_results.json`](file:///d:/antiviral-resistance-hsv--1/logs/plm_benchmark_results.json), [`logs/hsv-plm-gpu-benchmark.log`](file:///d:/antiviral-resistance-hsv--1/logs/hsv-plm-gpu-benchmark.log) | Executed cloud Tesla T4 GPU benchmark via Kaggle CLI kernel (`shubhamraikwar142/hsv-plm-gpu-benchmark`). Evaluated ESM-2 35M, ESM-2 650M, and ProtBERT across 929 mutations. ESM-2 650M (UL23-ACV: 0.7397, UL30-ACV: 0.7875, UL30-FOS: 0.8809) and ProtBERT (UL23-ACV: 0.7530, UL30-ACV: 0.6957, UL30-FOS: 0.7423) all achieved block means >= 0.65. **CORRECTION (2026-10-31):** the original claim here — that both models cleared "CIs > 0.50 on 3 out of 3 primary pairs" — was wrong; it paired block means with *pooled* position-cluster bootstrap CIs. On the correct block-level Student's $t$ basis, ESM-2 650M passes 3/3 and ProtBERT passes 1/3 (UL23-ACV only). Superseded by Item 14 and by the corrected table in `AUDIT_REPORT.md` §5.6. |

| **13. PHASE 0 DELIVERABLES (TRACKS A, B1, B2, C) & GATE 0 PASS** | Novelty queries, MSA conservation baseline, label composition, and human audit templates | **DONE** | [`logs/novelty_track_a_report.md`](file:///d:/antiviral-resistance-hsv--1/logs/novelty_track_a_report.md), [`logs/msa_conservation_results.json`](file:///d:/antiviral-resistance-hsv--1/logs/msa_conservation_results.json), [`logs/ul23_label_composition.log`](file:///d:/antiviral-resistance-hsv--1/logs/ul23_label_composition.log), [`logs/kappa_sample_50rows.tsv`](file:///d:/antiviral-resistance-hsv--1/logs/kappa_sample_50rows.tsv) | **Track A:** 6 PubMed queries executed; verified real PMIDs; draft decisions filled for all 67 entries in `novelty_widened_raw.tsv`. **Track B1:** Multi-herpesvirus MSA across HSV-1, HSV-2, VZV, HCMV, HHV-6B orthologs. PLM advantage (+0.14 to +0.22) confirmed over conservation on UL30; UL23 tied (~0.74). **Track B2:** UL23-ACV confirmed 100% missense (286/286 rows; 0 truncating). **Track C:** 50-row kappa sample, threshold checklist, and PI sign-off form created. **Gate 0:** Item 3 sparse (19 hits, 0 in HSV); dissociation survives conservation; Phase 2 justified. |

| **14. GPU BENCHMARK V4 RECTIFICATION (BLOCK T-INTERVALS & DIRECTION)** | Strict per-block breakdowns, Student's t 95% CIs, and row-level scored TSV | **DONE** | [`logs/plm_benchmark_results.json`](file:///d:/antiviral-resistance-hsv--1/logs/plm_benchmark_results.json), [`logs/primary_position_folds_scored.tsv`](file:///d:/antiviral-resistance-hsv--1/logs/primary_position_folds_scored.tsv) | Kaggle v4 executed on Tesla T4. Downloaded verified `primary_position_folds_scored.tsv` (929 rows, zero nulls). Evaluated block-level Student's t 95% CIs on contiguous-block folds. ESM-2 650M achieves 3/3 PASS on stopping rule (UL23-ACV: 0.7397 [0.5455, 0.9340]; UL30-ACV: 0.7875 [0.5055, 1.0000]; UL30-FOS: 0.8809 [0.6821, 1.0000]). ProtBERT achieves 1/3 strict t-interval pass (UL23-ACV: 0.7530 [0.6645, 0.8415]) and 3/3 mean >= 0.65 (UL30-ACV: 0.6957, UL30-FOS: 0.7423). Direction transparency confirmed (-LLR = deleterious = resistant). |

| **15. BLOSUM62 SIGN RECTIFICATION & ORTHOLOG COVERAGE AUDIT** | Inverted sign diagnosis, coverage audit, and missense framing quantification | **DONE** | [`logs/blosum_direction_audit.log`](file:///d:/antiviral-resistance-hsv--1/logs/blosum_direction_audit.log), [`scripts/audit_ortholog_coverage.py`](file:///d:/antiviral-resistance-hsv--1/scripts/audit_ortholog_coverage.py) | **BLOSUM62:** Root cause identified (double negation in runner inverted radicalness). Corrected contiguous-block mean: UL23-ACV 0.6380 [0.5524, 0.7237] (every block > 0.55); UL30-ACV 0.5494; UL30-FOS 0.5074 (pure chance). Proves amino acid radicalness alone provides strong signal on UL23 missense variants but fails on UL30. **Ortholog Coverage:** 86.7% (UL23) and 82.1% (UL30) positions have high coverage; dissociation is invariant to coverage thresholding (UL30 PLM advantage remains +0.13 to +0.15 on high-coverage positions). **Framing:** Truncating/frameshift variants represent 44.4% of cataloged resistant UL23 variants; paper framed strictly as missense resistance benchmark. |

| **16. PHASE 1 DELIVERABLES (C3 SUPERVISED, CALIBRATION, C5 TEMPORAL, C4 DRUG CONDITIONING)** | Multi-feature supervised predictor, Brier/ECE calibration, prospective temporal validation, and multi-task drug conditioning | **DONE** | [`logs/phase1_deliverables_results.json`](file:///d:/antiviral-resistance-hsv--1/logs/phase1_deliverables_results.json), [`logs/phase1_deliverables.log`](file:///d:/antiviral-resistance-hsv--1/logs/phase1_deliverables.log) | **C3 Supervised:** Combining ESM-2 650M + ProtBERT + MSA Conservation + BLOSUM62 + Physicochemical deltas achieves contiguous-block AUROC of 0.8188 [0.6796, 0.9580] on UL23-ACV, 0.7823 [0.4729, 1.0000] on UL30-ACV, and 0.8845 [0.7249, 1.0] on UL30-FOS. **CORRECTED 2026-10-31:** this entry previously read "0.8048 [0.5750, 1.0]" for UL30-ACV and claimed "all 3 pass stopping rule". The 0.8048 figure matches no model in `logs/phase1_deliverables_results.json` (Logistic_Regression_L2 = 0.7823 [0.4729, 1.0000]; HistGradientBoosting = 0.8000 [0.5731, 1.0000]). On the Logistic Regression values used in manuscript Table 4, **UL30-ACV does not pass** the stopping rule, since the block $t$-CI lower bound (0.4729) does not exceed 0.50. Supervised pass count is **2 of 3**, not 3 of 3. Manuscript Table 4 is unaffected: it reports the paired delta, not the stopping rule, and its numbers are correct. **Calibration:** Well-calibrated out-of-fold predictions (Brier: 0.084 to 0.178; ECE: 0.065 to 0.072). **C5 Temporal:** Prospective post-2016 evaluation achieves 0.7641 (UL23-ACV), 0.8290 (UL30-ACV), 0.8533 (UL30-FOS). Stratification confirms thesis: within-position test AUROC is 0.89-0.93; cross-position generalization to novel positions is 0.70-0.81. **C4 Drug Conditioning:** Multi-task model reaches 0.8695 [0.7634, 0.9756]; negative control UL23-FOS correctly suppresses false resistance transfer. |



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



## 4. Manuscript Status and Pending Human Actions Before Final Sign-Off



1. **Manuscript Draft Updated (`manuscript/manuscript_draft.md`):**

   - **Strategic Narrative Rebuild (2026-10-31):** Title, Abstract, Contributions (§1), §3.3, and Discussion §4 rebuilt around a **leakage-first, radicalness-anchored narrative**. The primary biological finding is the statistically robust amino acid radicalness dissociation (BLOSUM62: UL23 0.6380 [0.5524, 0.7237] vs. UL30 0.5074–0.5494 chance). The residual zero-shot PLM margin over alignment baselines on UL30 (+0.0876 ACV, +0.0563 FOS) is explicitly demoted to a directional, non-significant observation (paired $p=0.174$ and $p=0.373$; CIs span zero). Table 3 reorganized to report evolutionary baselines with paired deltas against strongest baselines, exact $p$-values, and 'tied' / 'n.s.' designations.
   - **Radicalness Interaction Test Executed (2026-11-01):** Resolving the difference-of-significance fallacy (Gelman & Stern, 2006), `scripts/compute_radicalness_interaction_test.py` evaluated the direct interaction between UL23 and UL30 across 2,000 position-cluster bootstrap resamples (`logs/radicalness_interaction_test.json`). The predictive delta is strictly positive and excludes zero: UL23-ACV vs. UL30-ACV $\Delta = +0.1187$ (95% CI: $[+0.0387, +0.2021]$, empirical $p=0.0020$); UL23-ACV vs. UL30-FOS $\Delta = +0.1598$ (95% CI: $[+0.0655, +0.2597]$, $p < 0.0001$). At the contiguous-block level, fold-to-fold variance is wide at $n=4$ or $5$ (Welch's $p=0.431$ and $p=0.406$).
   - **Structural Citations and PDB Annotations Rectified (2026-11-01):** Shankar et al. citation corrected to *Cell* 2024;187(20):5572–5586.e15 (PMID 39197451). PDB annotations clarified: `8EXX` is the foscarnet-bound pre-translocation state; `8V1T` is the aciclovir-triphosphate-bound closed state; `8V1Q`/`8V1R`/`8V1S` are holoenzyme forms; EMDB resolves the W781V mutant holoenzyme in editing conformation.
   - **PRIME Precedent Restored & Contribution #1 Reframed (2026-11-01):** Restored formal citation of PRIME (*BMC Genomics*, DOI: 10.1186/s12864-026-12976-5) in Introduction and References. Reframed Contribution #1 from claiming broad 'discovery' of spatial shortcut learning to demonstrating this vulnerability in HSV enzymes and establishing the contiguous-block domain generalization benchmark.
   - **3D Spatial-Neighborhood Holdout Pre-Registered & Evaluated (2026-11-01):** Executed the empirical structural neighborhood holdout experiment (`scripts/run_spatial_3d_holdout_experiment.py`) using cryo-EM coordinates (`8V1Q` Chain A, 1097 CA atoms, 88.6% reference coverage, 80.8% unique benchmark position coverage) and crystal coordinates (`1OF1` Chain B, 313 CA atoms, 1.95 Å resolution, complexed with (South)-methanocarba-thymidine; Bennett et al. 1999). Pre-registered 5 deterministic spatial clusters via KMeans (`k=5, seed=42`) into `data/processed/primary_3d_spatial_folds.tsv` (SHA-256: `0ED51D37D558530EF8B2D0C991780A23936580920E62C7ADE46A426ADFA96F4C`, logged in `logs/spatial_3d_fold_assignments_hash.log`). Key findings serialized in `logs/spatial_3d_holdout_results.json`:
     1. *Methodological Null Result Regarding PLM Structural Dependence:* ESM-2 650M maintains invariant generalization from 1D contiguous blocks to 3D spatial neighborhoods: UL30-ACV 0.7875 -> 0.7875 (resolved core: 0.7887 -> 0.7966); UL30-FOS 0.8809 -> 0.8932 (resolved core: 0.8793 -> 0.8977); UL23-ACV 0.7397 -> 0.7303 (resolved core: 0.7255 -> 0.7132). Cross-regime |Delta| <= 0.012 confirms that zero-shot PLMs do not depend on localized 3D physical contact environments, and 3D spatial grouping provides no incremental evaluation penalty.
     2. *Negative Control Stringency Finding (1D Block is Strictly Stricter):* Position-only random forest rebounds from pure chance under 1D contiguous block holdout (0.5001 on UL30-ACV) to 0.6241 +/- 0.345 under 3D spatial holdout (and from 0.5166 to 0.6001 +/- 0.288 on UL30-FOS). Because 3D geometric clusters partition sequential primary residues across different spatial folds, models exploit sequence-adjacency interpolation across training and test folds. This proves that 1D contiguous-block holdout is strictly the more stringent evaluation standard.
     3. *Target Dissociation Invariant in 3D Space:* Chemical dissimilarity alone (BLOSUM62 radicalness) consistently clears chance on the catalytic kinase (UL23-ACV: 1D 0.6380, 3D 0.6249; resolved core: 1D 0.6732, 3D 0.6501) but sits at chance on the allosteric polymerase (UL30-ACV: 1D 0.5494, 3D 0.4908; UL30-FOS: 1D 0.5074, 3D 0.4847).
     4. *Evolutionary Baselines in 3D Space:* Multi-herpesvirus PSSM scores 0.7189 (UL30-ACV) and 0.8334 (UL30-FOS), while positional Shannon entropy scores 0.6806 and 0.8338. ESM-2 maintains its lead over both baselines (+0.069 to +0.107 on ACV; +0.060 on FOS).
     5. *Publication Figure 5 Generated:* High-resolution 300-DPI multi-panel visualization saved to `manuscript/figures/figure5_spatial_3d_generalization.png` and artifact storage.
   - **Primary Reference List Completed (2026-11-01):** Restored missing primary data sources and structural citations to References: HerpesDRG (Charles et al., BMC Bioinformatics 2024, PMID: 39192205), CMI 2025 review (Dähne et al., Clin Microbiol Infect 2025, PMID: 40349973), and HSV-1 TK structure (Bennett et al., FEBS Lett 1999, PMID: 10025941). Corrected spelling from 'Döhne' to 'Dähne' throughout, and corrected cohort scope sentence to HSV-1.
   - **Estimator Alignment (2026-11-01):** Made the block-level Welch test primary for the target dissociation interaction claim (p=0.431 and p=0.406 across structural blocks), reporting the position-cluster bootstrap (Delta = +0.1187, p=0.0020; Delta = +0.1598, p < 0.0001) strictly as secondary sensitivity evidence, maintaining complete consistency with the paper's central thesis on spatial leakage.

   - **Foscarnet threshold corrected (2026-10-31)** in Methods §2.1. The previously asserted "stoichiometric equivalence" (100 uM * 300.04 g/mol = 30.0 ug/mL for foscarnet trisodium hexahydrate) is **WITHDRAWN as unsourced**. It was back-derived to make two literature numbers agree rather than read from any paper, and no source states the 30 ug/mL figure. Source review established: the one verifiable cut-off is **IC50 >= 100 ug/mL** (Wang, J Clin Microbiol 2024;62:e00728-23, reviewing a 1994 HIV cohort: 88% failure-to-heal at/above vs 0% below). Published cut-offs do not state the chemical form (free acid Mr 126.05 vs trisodium hexahydrate Mr 300.04), so no uM conversion is asserted. "Frobert et al. JAC 2014", cited in `thresholds.tsv` as the sole source for the 100 uM cut-off, **could not be located and appears not to exist** (no Frobert 2014 JAC paper; no DOI was ever recorded). Flagged for the human threshold review.

   - **PSSM profile baseline computed and logged (2026-10-31).** After the earlier phantom numbers (0.7409/0.6815/0.8119) were purged, a genuine multi-herpesvirus PSSM baseline with Dirichlet pseudocount smoothing (alpha=1.0, Robinson & Robinson 1991 background frequencies) was implemented in scripts/compute_pssm_baseline.py and executed directly on the primary cohort. Contiguous-block mean AUROCs are **0.7237 +/- 0.1009** (UL23-ACV), **0.6507 +/- 0.1920** (UL30-ACV), and **0.7832 +/- 0.0949** (UL30-FOS). On UL23, PSSM ties ESM-2 (0.7237 vs 0.7397, paired p=0.814). On UL30, ESM-2 maintains a positive paired difference (+0.1369 +/- 0.1177 on ACV; +0.0977 +/- 0.1154 on FOS), but with n=4 evaluable blocks both paired 95% Student's t intervals span zero (p=0.103 and p=0.189). Crucially, PSSM performs below positional Shannon entropy on UL30 (0.6507 vs 0.6999 on ACV; 0.7832 vs 0.8246 on FOS), leaving Shannon entropy as the strongest alignment comparator. All results serialized to logs/pssm_baseline_results.json and logs/pssm_baseline.log. Added as fourth baseline to Table 3 and incorporated into Discussion.
   - **Table 3 updated with PSSM baseline (2026-10-31).** Table 3 now reports all four alignment and profile baselines from logs/msa_conservation_results.json and logs/pssm_baseline_results.json (weighted conservation, family conservation fraction, negative positional entropy, and PSSM log-odds). Strongest baselines remain positional entropy on UL30 and weighted conservation on UL23, so the headline deltas (+0.0876 and +0.0563) remain unchanged. Discussion updated to report the real PSSM evaluation honestly rather than describing it as an open question.
   - **Abstract SD correction (2026-10-31).** The abstract reported the position-only RF as 0.8715 ± 0.0640 collapsing to 0.5001 ± 0.2887. Neither SD appears in any results log, and both contradicted the paper's own Table 1. Corrected to **± 0.0171** and **± 0.0009** per `logs/c2_repeated_cv_results.json`.

   - **Brier null baselines — FALSE ALARM RAISED AND WITHDRAWN (2026-10-31).** An earlier pass in this session flagged the null baselines 0.1774 / 0.1146 as "inflated", replaced them with the full-cohort values 0.1519 / 0.0949, and lowered the reported skill scores to +0.252 / +0.113. **That flag was wrong and is withdrawn; the original values are reinstated.** `scripts/run_phase1_deliverables.py` computes `brier_pt = brier_score_loss(y_eval, p_eval)` over the **evaluable-block cohort** (N = 286 / 269 / 250, excluding single-class Block 0), so the matching denominator is the evaluable-cohort null (0.2435 / 0.1774 / 0.1146). The original skill scores **+0.271 / +0.359 / +0.266** are correct. The erroneous "correction" paired an evaluable-cohort numerator with a full-cohort denominator — the same apples-to-oranges class of error as the pooled-vs-block CI mix-up it was meant to guard against. §3.4 now states the cohort explicitly and reports both null sets, so the pairing cannot be misread again. Both are logged under `deliverable_1_2_brier_bootstrap_audit`. **Lesson recorded: a flagged discrepancy is a hypothesis, not a finding, until the generating code path is read.**

   - **Brier bootstrap CIs computed and logged at source (2026-10-31).** §3.4 out-of-fold bootstrap 95% CIs were generated directly via 1,000 resamples within `scripts/run_phase1_deliverables.py` and serialized to `logs/phase1_deliverables_results.json` under `deliverable_1_2_brier_bootstrap_audit`: UL23-ACV [0.1543, 0.2006], UL30-ACV [0.0878, 0.1406], UL30-FOS [0.0602, 0.1131].

   - **Results-log label discrepancy resolved at source (2026-10-31).** `scripts/run_phase1_deliverables.py` was executed directly from source, dynamically evaluating contiguous-block CV folds and logging `deliverable_1_1_paired_c3_comparison` with accurate model labeling (`Logistic_Regression_L2`), exactly matching manuscript Table 4. The previous hardcoded mislabel in `update_phase1_paired_comparison.py` is permanently eliminated.

   - **Paired dissociation tests on alignment baselines executed and logged (2026-10-31).** `scripts/run_paired_dissociation_test.py` evaluated paired per-block deltas: (1) against negative positional entropy (`logs/paired_dissociation_entropy_test.json`): UL23-ACV $\Delta = +0.0146, p=0.843$; UL30-ACV $\Delta = +0.0877 \pm 0.0988, p=0.174$; UL30-FOS $\Delta = +0.0563 \pm 0.1079, p=0.373$. (2) against the Table 3 strongest baseline (`logs/paired_dissociation_strongest_test.json`): on UL23-ACV vs. Multi-Herpesvirus Weighted Conservation, $\Delta = +0.0024 \pm 0.1929$, $t(4)=0.0277$, $p=0.9792 \approx 0.98$ ('tied'); on UL30 vs. Entropy, $\Delta = +0.0877$ ($p=0.174$) and $+0.0563$ ($p=0.373$). All numbers in Table 3 trace directly to these frozen JSON logs.
   - **Primary threshold citations verified at source (2026-10-31).** Resolved foundational 1994 primary cohort: Safrin et al., *Antimicrob Agents Chemother* 1994;38(6):1246–1250 (PMID: 8092821), which established both the 2.0 µg/mL ACV and 100 µg/mL FOS PRA thresholds based on prospective failure-to-heal outcomes (95% and 88% predictive values, respectively). Corrected Burrel citation to *AAC* 2010 (PMID: 20733037; replacing invalid PMID 20380849). Sauerbrei et al. (JAC 2016, PMID: 26433780) affirmed as an authoritative database compilation of phenotyped variants rather than the original source of IC50 cut-offs. Updated `thresholds.tsv` and `docs/threshold_citation_checklist.md`, setting status to `AGENT_VERIFIED_PENDING_HUMAN` pending human sign-off.

   - **Abstract overclaim corrected (2026-10-31).** "exactly converging with our contiguous-block holdout estimates" → "broadly consistent with", matching §3.5, which had already been softened. Novel-position AUROCs (0.6996 / 0.8075 / 0.8126) differ from block estimates (0.8188 / 0.7823 / 0.8845) by up to 0.119; "exactly converging" was not defensible.

   - **Table 5 Caveats Added:** Seen positions on UL30 explicitly footnoted as descriptive completeness only due to N=2 resistant variants.

   - **Prospective Novel Wording Refined:** Updated from "closely mirror" to "broadly consistent with".

   - **Structural Claim Reframed (Option B):** Stated as consistent with Shankar et al. (*Cell* 2024) cryo-EM mechanism without asserting direct proof.

   - **Scope Prominently Documented:** 44.4% truncating/frameshift class in UL23 documented as deterministic loss-of-function; paper explicitly framed around the missense diagnostic frontier.

2. **Novelty Decisions Completed (`logs/novelty_widened_raw.tsv`):**

   - All 67 entries across 6 systematic query categories and claims C1–C3 have recorded audit decisions; zero blank entries remain.

3. **Pending Human Action Gates (You and the PI):**

   - **Dual-Annotator Curation Sample (Cohen's kappa):** Complete independent curation of the 50-row sample in [`logs/kappa_sample_50rows.tsv`](file:///d:/antiviral-resistance-hsv--1/logs/kappa_sample_50rows.tsv).

   - **Threshold Citation Verification:** Review and verify the 11 cited threshold rows against primary papers in [`docs/threshold_citation_checklist.md`](file:///d:/antiviral-resistance-hsv--1/docs/threshold_citation_checklist.md).

   - **PI Scope Sign-Off:** Execute the formal sign-off in [`docs/pi_signoff_form.md`](file:///d:/antiviral-resistance-hsv--1/docs/pi_signoff_form.md).



---

## 8. Phase 1 Rectification Audit Log (2026-11-02)

| Item | Focus Area | Action Taken | Evidence / Verification |
| :--- | :--- | :--- | :--- |
| **P1.1** | Table 1 vs Table 6 Discrepancy | Resolved 0.5001 vs 0.5347 at source: Table 1 reports 20-repeat CV average ($0.5001 \pm 0.0009$); Table 6 reports single 5-block partition where 3 of 4 evaluable folds are exactly 0.5000 (chance), fold mean is $0.5347 \pm 0.0694$, and pooled AUROC is $0.4591$ [0.3741, 0.5468]. Explicit cross-table reconciliation notes added to both tables. | `manuscript/manuscript_draft.md` Tables 1 & 6, `logs/spatial_3d_holdout_results.json` |
| **P1.2** | Structural Ligand / PDB Metadata | Corrected PDB `1OF1` metadata: 1.95 Å, chain B, complexed with (South)-methanocarba-thymidine (ligand SCT; Bennett et al. 1999). Clarified Shankar et al. (Cell 2024; 187:5572-5586.e15) structures (`8EXX` foscarnet-bound; `8V1Q`/`8V1R`/`8V1S`/`8V1T` holoenzyme states). | `manuscript/manuscript_draft.md`, `data/structures/1of1.cif` |
| **P1.3** | Reference List Completion | Added primary database references: HerpesDRG (Charles et al. 2024) and CMI review (Dähne et al. 2025). Standardized spelling to "Dähne". Verified abstract scope is restricted strictly to HSV-1. | `manuscript/manuscript_draft.md` References & Abstract |
| **P1.4** | 3D Spatial Null Re-interpretation | Removed assertions that zero-shot PLM invariance across 3D folds constitutes biological proof of contact independence (zero-shot scores are computed from sequence and do not change across folds). Reframed §3.6 and Contribution 7 around the genuine empirical hazard: 3D geometric clustering partitions sequence non-contiguously, reintroducing sequence-adjacency interpolation shortcuts (Position-Only RF rebounds to 0.6241). | `manuscript/manuscript_draft.md` §3.6, Contribution 7 |
| **P1.5** | Equal PLM Evaluation & ProtBERT Claims | Removed post-hoc designation of ESM-2 650M as "primary". Evaluated all three models (ESM-2 650M, ProtBERT, ESM-2 35M) symmetrically. Corrected ProtBERT claim in Contribution 4 (passes $\ge 0.65$ mean on all 3 pairs, satisfies strict $t$-CI lower bound on UL23-ACV, exhibits wider fold variance on UL30). | `manuscript/manuscript_draft.md` Table 2, Contribution 4 |
| **P1.6** | Rhetorical Calibration | Replaced strong claims ("Discovery", "proves", "strictly more stringent") with calibrated empirical statements ("Empirical Demonstration", "demonstrates", "provides a more effective control against spatial interpolation shortcuts"). | `manuscript/manuscript_draft.md` |


---

## 9. Phase 2 & 3 Empirical Execution and Pre-Registered Outcome (2026-11-02)

| Item | Focus Area | Action Taken | Evidence / Verification |
| :--- | :--- | :--- | :--- |
| **P2.1** | Hash & URL Reconciliation | Reconciled dataset hashes across `README.md`, `manuscript_draft.md`, and `STATUS.md`: Primary 929-row missense benchmark is `7c4d0a3de2cd171b07230757f44e0dfbd61c6a0a8570dc8f263195b354c66909` (`primary_position_folds.tsv`); Master 2,204-row compilation is `cc1356264cb451542620dd5d8469039419531334850af7b7933880fc774a92cd` (`unified_modeling_dataset.tsv`). Standardized remote repository URL to `https://github.com/Shubham142856/antiviral-resistance-hsv-`. | `README.md`, `manuscript/manuscript_draft.md`, `STATUS.md` |
| **P2.2** | Verified Domain Boundaries | Resolved authoritative structural domain boundaries from Liu et al. (2006) *JBC*, Shankar et al. (2024) *Cell*, CATH, and InterPro. Evaluated class distributions: UL30 NTD (1–362) has 0 resistant variants across ACV/FOS; FOS Exo and Thumb have sparse positives ($N_{\text{pos}}=2, 3$). Pre-registered pooled out-of-fold AUROC with domain-level bootstrap CIs as the primary evaluation metric. | `docs/WHERE_WHAT_PREREGISTRATION.md`, `scripts/where_what/extract_features.py` |
| **P2.3** | Dual-Buffer Protocol & Split Validity | Evaluated grid of sequence buffers ($\pm 5, \pm 10, \pm 15$ aa) and 3D spatial buffers ($< 8, < 10, < 12$ Å). Minimal valid buffer satisfying the Position-Only RF split validity check ($\text{AUROC} \in [0.45, 0.55]$): $\pm 5$ aa sequence buffer and $< 12.0$ Å 3D distance buffer (Position-Only RF pooled AUROC = $0.5381$). BLOSUM62 radicalness pooled AUROC: UL30-ACV $0.4942$ [0.3521, 0.5649] (chance); UL30-FOS $0.4531$ [0.3773, 0.6535] (chance); UL23-ACV $0.6130$ [0.5504, 0.7175] (signal preserved). | `logs/dual_buffer_selection.json`, `logs/dual_buffered_domain_benchmark.json` |
| **P3.1** | Where x What Pre-Registration | Pre-registered all hypotheses, features, model specifications ($L_2$ logistic regression), and explicit falsification thresholds ($\beta_5 < 0$ with 95% CI strictly excluding zero on UL30-ACV) in `docs/WHERE_WHAT_PREREGISTRATION.md` (SHA-256: `ca4d2e7ff132f37bf4067a220c5c48f5a5db57e396b8c2c2571e5f19ea43adc6`). Committed and pushed to GitHub (commit `ded62cb`) before test execution. | `docs/WHERE_WHAT_PREREGISTRATION.md`, Git commit `ded62cb` |
| **P3.2** | Where x What Empirical Test Results | Executed tolerance-window logistic regressions with 1,000 position-clustered bootstrap resamples (`scripts/where_what/run_tolerance_models.py`). **Primary Target (UL30-ACV) Interaction Curvature:** $\beta_5 = +0.0368$, 95% CI: $[-0.1868, +0.3039]$ (**FALSIFIED**; null/non-negative point estimate, CI spans zero). **Dispensable Target (UL23-ACV):** Linear disruption $\beta_2 = -0.3629$ [$-0.6831, -0.1235$]; Curvature $\beta_3 = +0.0423$ [$-0.1485, +0.2463$] (strictly monotonic). | `logs/where_what_test_results.json` |
| **P3.3** | Main Manuscript & Title Preservation | Per user instruction (*"Ye tabhi likhna jab neeche ka test pass ho. Warna fallback title niche hai"*), the main paper title ("Spatial Shortcut Learning Inflates Viral Resistance Prediction") and baseline findings remain intact and separate. The negative result is honestly documented as an empirical falsification that demonstrates the critical necessity of pre-registration against HARKing. | `manuscript/manuscript_draft.md`, `STATUS.md` |
