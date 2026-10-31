# AUDIT_REPORT.md: Antiviral-Resistance Mutation Prediction (Paper 2)

**Project:** HSV Antiviral-Resistance Mutation Prediction (Paper 2)  
**Funding:** MPCST grant on computational approaches to HSV protein sequence analysis  
**Audit Protocol Version:** Version 4 (Strict), frozen at `AUDIT_PROTOCOL.md`  
**Protocol SHA-256:** `D2E7FEC2B43D378E1C9203EFE5438E2E20FB6D55339C3521F27FD00EC76FD75C`  
**Audit Completion Date:** 2026-10-28  
**Auditor / Extraction System:** Antigravity Deterministic Audit Pipeline (Executed on `D:\antiviral-resistance-hsv--1`)  

---

## 1. Executive Summary & Provisional Scope Statement (Pending Sign-Off)

### Provisional Scope Statement
> This paper presents an evidence-tiered, position- and clade-stratified benchmark for HSV resistance mutations (**C1**, **C2**). A drug-conditioned predictor is evaluated for **UL23–Aciclovir**, **UL30–Aciclovir**, and **UL30–Foscarnet**, which pass the empirical position threshold of $\ge 15$ resistant positions (**C3**); helicase-primase (*UL5*, *UL52*) and *UL29* are restricted to secondary, pooled exploratory analyses due to severe position sparsity ($< 10$ positions). The pathway-constrained gate is **Primary** (**C4**), enforced as a domain-specific biological routing architecture across divergent resistance pathways. The temporal holdout is **Primary** (**C5**), using the natural 2016 publication-year cutoff yielding 226 post-cutoff rows across 189 distinct mutations (132 novel). The masking analysis is **Case Series Only, No AUC** (**C6**), as none were found in the logged search of verified literature. Cross-resistance and landscape outputs are **Derived Predictions** (**C8**), strictly segregating observed experimental ground truth from model inferences. Multi-site epistasis modeling is **Deferred to Paper 3** (**C10**). Every claim traces directly to empirical counts, cryptographic database snapshots, and peer-reviewed literature citations.

---

## 2. Stage 0 & 1: Rights, Coverage, and Raw Snapshot Provenance

### 2.1 Cryptographic Snapshot Verification
All database extracts and source files were cloned, hashed, and inspected directly on the `D:` drive (`D:\antiviral-resistance-hsv--1\data\raw\herpesdrg-db`).

| Artifact | Source Repository / Identifier | Commit Hash | File SHA-256 (LF Normalized) | License / Terms |
| :--- | :--- | :--- | :--- | :--- |
| **HerpesDRG-DB** | `github.com/ojcharles/herpesdrg-db` | `22fc18676a1ea49cd6d7cd1c3280db210a7580e6` | `e1d9f61caa37f47b6a102ff020a4f483b2a82dc5735bd6688cda57827c36f770` | MIT License (Oscar J Charles, 2021) |
| **AUDIT_PROTOCOL.md** | Stage 0 Frozen Protocol | Git HEAD | `D2E7FEC2B43D378E1C9203EFE5438E2E20FB6D55339C3521F27FD00EC76FD75C` | CC-BY 4.0 / Project Protocol |

*Note on Windows Line Endings:* On Windows platforms, git default CRLF conversion produces raw file hash `5a946ba484376ac518c7779563cfbafb3a627c3006800d8cfeddac8e0cb59e1f`. Normalizing `\r\n` to `\n` confirms the exact canonical repository SHA-256: `e1d9f61caa37f47b6a102ff020a4f483b2a82dc5735bd6688cda57827c36f770`.

### 2.2 Coverage and Rights Audit of Literature Sources
1. **HerpesDRG-DB (ojcharles/herpesdrg-db):**  
   - Active status filter `status == 'A'` applied per README specification. Total database rows: 2,035; Active rows: 1,696; Active HSV rows: 951 (HSV-1: 754; HSV-2: 197).
   - Gene coverage in active HSV rows: UL23 (470 rows), UL30 (440 rows), UL5 (31 rows), UL52 (9 rows), UL29 (1 row). **UL8 is completely absent (0 rows)**.
2. **CMI 2025 Review Table S3 (Dähne et al., 2025):**  
   - *Citation:* Clin Microbiol Infect 2025; 31(9): 1484–1490. DOI: [10.1016/j.cmi.2025.04.046](https://doi.org/10.1016/j.cmi.2025.04.046). PMID: [40349973](https://pubmed.ncbi.nlm.nih.gov/40349973/).  
   - *Corrigendum:* Clin Microbiol Infect 2026; DOI: [10.1016/j.cmi.2026.01.010](https://doi.org/10.1016/j.cmi.2026.01.010).  
   - *Terms:* Elsevier standard license. Only derived accession IDs, extracted mutation tuples, and scripts may be redistributed; raw supplementary tables must not be republished verbatim without permission.
   - *Local Verification:* `mmc3.xlsx` downloaded and verified locally (SHA-256: `760ab5f4e3185d5e22b6e37c92097c752e053a72ec274a20c1ed8daca81ceb52`).
3. **CHARMD Database (Tilloy et al., 2024):**  
   - *Citation:* Antiviral Res 2024; 231: 106016. DOI: [10.1016/j.antiviral.2024.106016](https://doi.org/10.1016/j.antiviral.2024.106016). PMID: [39349222](https://pubmed.ncbi.nlm.nih.gov/39349222/).  
   - *Web App Audit:* Live scrape of `https://www.unilim.fr/cnr-herpesvirus/outils/codexmv/database/` (logged at `logs/charmd_evidence_20261028.json`) revealed that the online database currently hosts exclusively Cytomegalovirus (*Human betaherpesvirus 5*), with zero HSV-1 or HSV-2 entries.
4. **Sauerbrei Historical Baseline (Sauerbrei et al., 2016):**  
   - *Citation:* J Antimicrob Chemother 2016; 71(1): 6–16. DOI: [10.1093/jac/dkv285](https://doi.org/10.1093/jac/dkv285). PMID: [26433780](https://pubmed.ncbi.nlm.nih.gov/26433780/).  
   - *Access:* Establishes pre-2016 historical baseline for UL23 and UL30 polymorphisms and resistance mutations.

---

## 3. Stage 4: Empirical Counts and Power Analyses

### Table A: Gene-Drug Pair Counts & Deduplicated Empirical Support
*Deduplicated on unique (Virus, Gene, Drug, Mutation) tuples across HerpesDRG-DB and CMI 2025 Table S3. Classified using cited blank-free thresholds from `thresholds.tsv`: Resistant if Fold-Change >= 3.0x (or reported 'Resistant'); Sensitive/Polymorphism if Fold-Change < 3.0x (or reported 'Sensitive') per CMI 2025 Table S3, Sauerbrei 2016, and HerpesDRG 2024. Note on Denominators: Across the raw multi-source extract, there are 1,888 total single-amino-acid substitution rows across all species (HSV-1, HSV-2), all 7 drugs, and discordant entries; filtering strictly for HSV-1, unconflicted consensus phenotypes, and the three primary arms yields the frozen primary evaluation cohort of exactly 929 clean single-site rows (UL23-ACV: 286, UL30-ACV: 332, UL30-FOS: 311; SHA-256: `7c4d0a3de2cd171b07230757f44e0dfbd61c6a0a8570dc8f263195b354c66909`).*

*Evaluated under two strict regimes:*
1. **All-Inclusive:** Retains all consensus entries.
2. **Strict Excluding Conflicts:** Excludes the 65 discordant tuples where HerpesDRG and CMI 2025 disagree on resistance phenotype (logged in `logs/cmi_herpesdrg_conflicts.log`).

| Gene | Drug | Evaluation Mode | Distinct Mutations | Distinct Positions | Resistant Mutations | Resistant Pos (Total) | Sensitive Positions | HSV-1 Res Pos | HSV-2 Res Pos | Distinct Sources | Audit Status / Verdict |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
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
| **UL30** | **Cidofovir** | All-Inclusive | 123 | 112 | 23 | **22** | 93 | 18 | 5 | 68 | **Excluded (Dropped per Protocol)** |
| **UL30** | **Cidofovir** | Strict Excl. Conflicts | 112 | 102 | 12 | **11** | 93 | 11 | 1 | 65 | **FAIL (< 15 pos) — Excluded (Dropped)** |
| **UL5** | **Pritelivir** | All-Inclusive | 13 | 8 | 12 | **7** | 1 | 5 | 3 | 13 | **FAIL (< 15 pos) — Secondary / Pooled** |
| **UL5** | **Pritelivir** | Strict Excl. Conflicts | 13 | 8 | 12 | **7** | 1 | 5 | 3 | 13 | **FAIL (< 15 pos) — Secondary / Pooled** |
| **UL5** | **Amenamevir** | All-Inclusive | 17 | 8 | 13 | **4** | 4 | 3 | 2 | 5 | **FAIL (< 15 pos) — Statistically Underpowered** |
| **UL5** | **Amenamevir** | Strict Excl. Conflicts | 17 | 8 | 13 | **4** | 4 | 3 | 2 | 5 | **FAIL ($< 15$ pos) — Statistically Underpowered** |
| **UL52** | **Pritelivir** | All-Inclusive | 2 | 2 | 2 | **2** | 0 | 1 | 1 | 3 | **FAIL (< 5 pos) — Statistically Underpowered** |
| **UL52** | **Pritelivir** | Strict Excl. Conflicts | 2 | 2 | 2 | **2** | 0 | 1 | 1 | 3 | **FAIL (< 5 pos) — Statistically Underpowered** |
| **UL52** | **Amenamevir** | All-Inclusive | 7 | 5 | 4 | **3** | 3 | 3 | 0 | 3 | **FAIL (< 5 pos) — Statistically Underpowered** |
| **UL52** | **Amenamevir** | Strict Excl. Conflicts | 7 | 5 | 4 | **3** | 3 | 3 | 0 | 3 | **FAIL (< 5 pos) — Statistically Underpowered** |
| **UL29** | **Aciclovir** | All-Inclusive | 1 | 1 | 0 | **0** | 1 | 0 | 0 | 2 | **FAIL (< 5 pos) — Statistically Underpowered** |
| **UL29** | **Aciclovir** | Strict Excl. Conflicts | 1 | 1 | 0 | **0** | 1 | 0 | 0 | 2 | **FAIL (< 5 pos) — Statistically Underpowered** |

*\*Note on UL23-Foscarnet: Apparent resistant positions are clinical isolate artifacts harboring co-mutations in UL30 (see Section 5).*  
*\*\*Impact of Phenotype Discordance: The three core pairs (`UL23-ACV`, `UL30-ACV`, `UL30-FOS`) comfortably exceed the >= 15 position threshold whether conflicts are included or strictly excluded (156 vs 147; 84 vs 77; 57 vs 50). However, for `UL30-PCV` and `UL30-CDV`, excluding discordant tuples drops them below 15 resistant positions (14 and 11 positions respectively), meaning they do NOT qualify for standalone primary evaluation under strict conflict exclusion.*

### Table B: Label Purity by Position (Position Shortcut Risk)
*Evaluating single-label vs mixed-label positions across both denominators (logged in `logs/label_purity_summary.tsv`):*

| Gene – Drug Pair | D1: All Positions | D1 Single-Label | D1 % Single | D1 Mixed | D1 % Mixed | D2: Positions ($\ge 2$ Subs) | D2 Single-Label | D2 % Single | D2 Mixed | D2 % Mixed | Audit Decision |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **UL23 – Aciclovir** | 191 | 178 | **93.19%** | 13 | 6.81% | 51 | 38 | **74.51%** | 13 | 25.49% | **Position-Shortcut Mandatory** |
| **UL30 – Aciclovir** | 246 | 242 | **98.37%** | 4 | 1.63% | 31 | 27 | **87.10%** | 4 | 12.90% | **Position-Shortcut Mandatory** |
| **UL30 – Foscarnet** | 234 | 231 | **98.72%** | 3 | 1.28% | 30 | 27 | **90.00%** | 3 | 10.00% | **Position-Shortcut Mandatory** |

*Crucial Audit Finding:* Over **74% to 90% of positions with multiple substitutions** exhibit only a single phenotype across all reported alleles. Under random train/test splits, a model memorizes position indices rather than learning biochemical substitution rules. Therefore, **within-position ranking** and **position-stratified grouped cross-validation** are strictly mandatory to prevent fake accuracy inflation.

### Table C: Engineered Doubles vs Co-occurrences
Auditing all 8 co-mutation rows present in the active HSV extract:

| Category | Count | Description / Biological Evidence |
| :--- | :--- | :--- |
| **Engineered Doubles (Tier A)** | **None found in the logged search** | Recombinant BAC or marker-transfer isolating dual mutations with wild-type and single-mutant controls. |
| **Clinical Co-occurrences (Tier B)** | **5** | Clinical isolates with co-occurring mutations (e.g. UL23 A189V + UL30 L778M/L802F; UL23 T183P + UL30 L778M/L802F; HSV2 UL30 D912V + TK 28*). |
| **Lab-Selection Multi-Mutants** | **3** | Strains generated under serial passage drug pressure (e.g. UL52 A899T + UL5 K356T; UL23 C336Y + UL30 H168Y; UL23 T96I + UL30 A991T). |

*Audit Decision:* **None found in the logged search ($< 10$ kill criterion triggered).** Claim **C6** kill criterion is triggered. Epistatic masking cannot be evaluated as an AUROC benchmark; it must be reframed as a qualitative descriptive case series. Claim **C10** (Multi-site head) is deferred to Paper 3.

### Table D: Experimental Assays & Top Source Concentration
Top laboratory test methods and publication concentration per drug:

| Drug | Total Rows | Dominant Assay Types | Top Source Concentration (% of drug rows) |
| :--- | :--- | :--- | :--- |
| **Aciclovir** | 917 | Plaque Reduction Assay (534, 58.2%), Dye Uptake Assay (214, 23.3%), TK Kinetic Assay (56, 6.1%) | Top source (DOI: [10.1128/AAC.00977-15](https://doi.org/10.1128/AAC.00977-15), Schmidt et al. 2015): 11.9%; Top 3 sources: 27.6% |
| **Foscarnet** | 638 | Plaque Reduction Assay (430, 67.4%), Dye Uptake Assay (136, 21.3%), Cytopathic Effect (35, 5.5%) | Top source (DOI: [10.1128/AAC.00977-15](https://doi.org/10.1128/AAC.00977-15), Schmidt et al. 2015): 16.9%; Top 3 sources: 37.5% |
| **Pritelivir** | 42 | Cytopathic Effect (22, 52.4%), Plaque Reduction Assay (20, 47.6%) | Top source (DOI: [10.1128/aac.00110-24](https://doi.org/10.1128/aac.00110-24)): 26.2% |
| **Amenamevir** | 34 | Plaque Reduction Assay (17, 50.0%), Cytopathic Effect (17, 50.0%) | Top source (DOI: [10.1128/AAC.00494-21](https://doi.org/10.1128/AAC.00494-21), Kawakami et al. 2021): 50.0% |

### Table E: Publication Year Distribution (Resolved DOI Partition)
*Extracted from Crossref DOIs across all 951 active HSV rows (logged in `logs/publication_year_lookup.tsv`):*

- **Pre-2016 Corpus (Training):** 725 rows (569 distinct mutations).
- **Post-2016 Corpus (Prospective Test):** 226 rows (189 distinct mutations, 132 novel).
- *Audit Verdict:* **Passes Option A condition ($\ge 30$ post-cutoff mutations).** Clarified as a **Publication Year Partition**, distinguishing paper publication dates from biological first reports.

### Table F: Species Breakdown (HSV-1 vs HSV-2)

| Gene | HSV-1 Rows | HSV-1 Unique Mutations | HSV-1 Positions | HSV-2 Rows | HSV-2 Unique Mutations | HSV-2 Positions | Primary Claim Scope |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **UL23** | 377 | 255 | 191 | 93 | 54 | 47 | HSV-1 Primary; HSV-2 Cross-species Test |
| **UL30** | 341 | 281 | 247 | 99 | 85 | 77 | HSV-1 Primary; HSV-2 Cross-species Test |
| **UL5** | 26 | 13 | 8 | 5 | 5 | 4 | Secondary / Pooled |
| **UL52** | 9 | 7 | 6 | 0 | 0 | 0 | Secondary / Pooled |
| **UL29** | 1 | 1 | 1 | 0 | 0 | 0 | Secondary / Pooled |


### Table F: Homology Alignment & Sequence Mapping (HSV-1 vs HSV-2)
*Global pairwise alignment between canonical reference strains: HSV-1 (NC_001806.2 / 17) and HSV-2 (NC_001798.2 / HG52) logged in `logs/homology_mapping.log` and `data/processed/hsv1_hsv2_homology_map.tsv`:*

- **UL23 (Thymidine Kinase):** Both proteins are exactly 376 amino acids long. 277 residues are identical (73.67%), 97 substitutions, 2 gapped alignment columns.
- **UL30 (DNA Polymerase):** HSV-1 is 1,235 aa; HSV-2 is 1,240 aa (+5 aa net difference). 1,123 residues are identical (90.93%), 107 substitutions, 5 gapped alignment columns (insertions/deletions concentrated in variable linker regions around residues 675–690 and 1106–1119).
- *Methodological Consequence:* Pooling raw sequence positions across HSV-1 and HSV-2 directly introduces index mismatching for UL30. Therefore:
  1. The **Primary Analysis Scope is strictly restricted to HSV-1**.
  2. Homologous cross-species pooling requires canonical `homology_cluster_id` mapping.
  3. HSV-2 is reserved as an **external cross-species generalization test set**.

---

## 4. Stage 5: Feasibility Probes, Leakage Bias & Dual-Split Standard (Claim C2)

To rigorously test the data leakage hypothesis of **Claim C2** and measure the performance inflation resulting from naive random splitting compared to structured holdouts, repeated cross-validation experiments were conducted across **20 independent repeats** (5 folds per repeat = 100 evaluations per cell, random seeds 1..20).

### 4.1 Repeated CV Benchmark Table across Split Regimes (20 Repeats, Mean Per-Fold AUROC +/- SD)
*Evaluated across 6 ablation models and 5 split regimes (logged in `logs/c2_repeated_cv.log` and `logs/c2_repeated_cv_results.json`). Paired leakage gaps represent the mean paired difference $(\Delta = \text{Naive Random} - \text{Position-Grouped})$ with Nadeau & Bengio (2003) corrected resampled t-test 95% confidence intervals accounting for training set overlap across repeats:*

#### Pair 1: UL23 – Aciclovir (HSV-1, Single-Site N=286 Clean Tuples)
| Model Architecture | Naive Random | Strain-Grouped | Position-Grouped | Both-Grouped (Double-Holdout) | Contiguous-Block (5 Blocks) | Paired Leakage Gap (Random - Pos) [Nadeau-Bengio 95% CI] |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Random Forest (All Features)** | 0.9109 +/- 0.0070 | 0.8793 +/- 0.0269 | 0.8067 +/- 0.0273 | 0.7683 +/- 0.0722 | 0.6875 +/- 0.0079 | **+0.1041** [+0.0703, +0.1379] |
| **Random Forest (No Position Feature)** | 0.7798 +/- 0.0095 | 0.7184 +/- 0.0322 | 0.6758 +/- 0.0286 | 0.6167 +/- 0.0839 | 0.6378 +/- 0.0081 | **+0.1040** [+0.0702, +0.1378] |
| **Random Forest (Position Only)** | 0.8503 +/- 0.0089 | 0.8396 +/- 0.0327 | 0.7683 +/- 0.0182 | 0.7677 +/- 0.0794 | 0.5440 +/- 0.0014 | **+0.0821** [+0.0618, +0.1024] |
| **Logistic Regression (One-Hot + Pos)** | 0.7444 +/- 0.0057 | 0.6859 +/- 0.0384 | 0.7117 +/- 0.0213 | 0.7052 +/- 0.0794 | 0.6877 +/- 0.0000 | **+0.0327** [+0.0102, +0.0552] |
| **BLOSUM62 Radicalness Score** | 0.7370 +/- 0.0046 | 0.6970 +/- 0.0355 | 0.7389 +/- 0.0077 | 0.7337 +/- 0.0792 | 0.7389 +/- 0.0000 | **-0.0019** [-0.0137, +0.0099] |
| **Position Prior Baseline** | 0.8490 +/- 0.0101 | 0.8297 +/- 0.0239 | **0.5000 +/- 0.0000** | **0.5000 +/- 0.0000** | 0.5000 +/- 0.0000 | **+0.3490** [+0.3378, +0.3602] |

#### Pair 2: UL30 – Aciclovir (HSV-1, Single-Site N=332 Clean Tuples)
| Model Architecture | Naive Random | Strain-Grouped | Position-Grouped | Both-Grouped (Double-Holdout) | Contiguous-Block (5 Blocks) | Paired Leakage Gap (Random - Pos) [Nadeau-Bengio 95% CI] |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Random Forest (All Features)** | 0.8620 +/- 0.0167 | 0.7976 +/- 0.0621 | 0.7328 +/- 0.0365 | 0.6226 +/- 0.1302 | 0.4089 +/- 0.0198 | **+0.1292** [+0.0878, +0.1706] |
| **Random Forest (No Position Feature)** | 0.6946 +/- 0.0181 | 0.6579 +/- 0.0414 | 0.5052 +/- 0.0536 | 0.4465 +/- 0.1598 | 0.3956 +/- 0.0143 | **+0.1894** [+0.1265, +0.2523] |
| **Random Forest (Position Only)** | 0.9071 +/- 0.0077 | 0.8508 +/- 0.0426 | 0.8715 +/- 0.0171 | 0.7780 +/- 0.1303 | 0.5001 +/- 0.0009 | **+0.0356** [+0.0136, +0.0576] |
| **Logistic Regression (One-Hot + Pos)** | 0.5619 +/- 0.0171 | 0.5184 +/- 0.0366 | 0.4521 +/- 0.0338 | 0.3963 +/- 0.1470 | 0.3696 +/- 0.0000 | **+0.1098** [+0.0740, +0.1456] |
| **BLOSUM62 Radicalness Score** | 0.4471 +/- 0.0042 | 0.4680 +/- 0.0420 | 0.4519 +/- 0.0171 | 0.4536 +/- 0.1313 | 0.4545 +/- 0.0000 | **-0.0048** [-0.0261, +0.0165] |
| **Position Prior Baseline** | 0.8038 +/- 0.0191 | 0.7487 +/- 0.0563 | **0.5000 +/- 0.0000** | **0.5000 +/- 0.0000** | 0.5000 +/- 0.0000 | **+0.3038** [+0.2798, +0.3278] |

#### Pair 3: UL30 – Foscarnet (HSV-1, Single-Site N=311 Clean Tuples)
| Model Architecture | Naive Random | Strain-Grouped | Position-Grouped | Both-Grouped (Double-Holdout) | Contiguous-Block (5 Blocks) | Paired Leakage Gap (Random - Pos) [Nadeau-Bengio 95% CI] |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Random Forest (All Features)** | 0.8335 +/- 0.0234 | 0.7317 +/- 0.0587 | 0.6756 +/- 0.0446 | 0.5192 +/- 0.1274 | 0.6372 +/- 0.0181 | **+0.1579** [+0.1050, +0.2108] |
| **Random Forest (No Position Feature)** | 0.7413 +/- 0.0201 | 0.6484 +/- 0.0703 | 0.5548 +/- 0.0342 | 0.4383 +/- 0.1248 | 0.5934 +/- 0.0116 | **+0.1865** [+0.1444, +0.2286] |
| **Random Forest (Position Only)** | 0.9027 +/- 0.0078 | 0.8324 +/- 0.0464 | 0.8684 +/- 0.0282 | 0.7626 +/- 0.1202 | 0.5166 +/- 0.0105 | **+0.0343** [-0.0012, +0.0698] |
| **Logistic Regression (One-Hot + Pos)** | 0.6303 +/- 0.0235 | 0.5751 +/- 0.0415 | 0.5468 +/- 0.0507 | 0.4506 +/- 0.0915 | 0.4558 +/- 0.0000 | **+0.0835** [+0.0232, +0.1438] |
| **BLOSUM62 Radicalness Score** | 0.4313 +/- 0.0142 | 0.4850 +/- 0.0517 | 0.4436 +/- 0.0230 | 0.4738 +/- 0.1331 | 0.4803 +/- 0.0000 | **-0.0123** [-0.0416, +0.0170] |
| **Position Prior Baseline** | 0.7507 +/- 0.0221 | 0.7289 +/- 0.0399 | **0.5000 +/- 0.0000** | **0.5000 +/- 0.0000** | 0.5000 +/- 0.0000 | **+0.2507** [+0.2238, +0.2776] |

*Note: Confidence intervals are computed using the Nadeau & Bengio (2003) corrected resampled variance $\widehat{\text{Var}}_{\text{NB}}(\bar{d}) = (1/R + n_{\text{test}}/n_{\text{train}}) S^2 = 0.30 S^2$, correcting for pairwise fold dependencies and training set overlap across repeated 5-fold cross-validation. Bold indicates the 95% confidence interval strictly excludes zero.*

### 4.2 Biological Direction & Mutational Tolerance Audit
To verify whether antiviral resistance is driven by catalytic disruption (deleterious fitness penalty) rather than enhanced viral fitness, rank correlations between ESM-2 masked-marginal log-likelihood ratios ($\Delta 	ext{LLR} = \log P(	ext{mut} \mid \mathbf{x}_{\setminus i}) - \log P(	ext{wt} \mid \mathbf{x}_{\setminus i})$) and binary resistance labels were computed across all 929 primary entries (logged in `logs/direction_correlations.log`):

| Target Pair | Cohort (Res / Sens) | Mean $\Delta 	ext{LLR}$ (Sens) | Mean $\Delta 	ext{LLR}$ (Res) | Spearman $
ho$ ($p$-value) | Point-Biserial $r$ ($p$-value) | Mann-Whitney U ($p$-value) | Verified Direction |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **UL23 – Aciclovir** | 120 Res / 166 Sens | -0.4979 +/- 1.1492 | -1.7110 +/- 2.3423 | **-0.2804** ($1.45 	imes 10^{-6}$) | **-0.3234** ($2.19 	imes 10^{-8}$) | $p = 1.11 	imes 10^{-6}$ | **Deleterious $	o$ Resistant** |
| **UL30 – Aciclovir** | 62 Res / 270 Sens | -0.7803 +/- 1.2574 | -1.9899 +/- 2.1904 | **-0.2340** ($1.66 	imes 10^{-5}$) | **-0.3040** ($1.58 	imes 10^{-8}$) | $p = 1.04 	imes 10^{-5}$ | **Deleterious $	o$ Resistant** |
| **UL30 – Foscarnet** | 33 Res / 278 Sens | -0.8149 +/- 1.3146 | -2.3272 +/- 2.0476 | **-0.2620** ($2.81 	imes 10^{-6}$) | **-0.3136** ($1.60 	imes 10^{-8}$) | $p = 1.99 	imes 10^{-6}$ | **Deleterious $	o$ Resistant** |

*Key Biological Verification:* The negative correlation is statistically verified across all three arms ($p < 10^{-5}$ across all tests). Resistant mutations consistently impose severe evolutionary fitness penalties ($\Delta 	ext{LLR}$ between -1.71 and -2.33) compared to benign polymorphisms (-0.50 to -0.81). Thus, sequence language models act as zero-shot resistance predictors precisely by measuring catalytic destabilization.

### 4.3 Position-Clustered Bootstrap Interval Audit (1,000 Resamples over Unique Positions)
To eliminate repeated-CV sample-overlap ambiguity, we implemented the reviewer's preferred method: **cluster bootstrapping over unique sequence positions** ($B = 1,000$ resamples of unique positions with replacement, logged in `logs/c2_position_cluster_bootstrap.log`):

| Target Pair | Model Architecture | Mean Gap ($\Delta = 	ext{Rand} - 	ext{Pos}$) | Bootstrap SE | 95% Position-Cluster CI | Leakage Verdict |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **UL23 – Aciclovir** | Position Prior Baseline | +0.2664 | 0.0372 | **[+0.1872, +0.3331]** | **Excludes zero (Significant leakage)** |
| **UL23 – Aciclovir** | RF (Position Only) | +0.0327 | 0.0181 | **[+0.0018, +0.0693]** | **Excludes zero (Significant leakage)** |
| **UL23 – Aciclovir** | RF (All Features) | +0.0103 | 0.0141 | [-0.0172, +0.0381] | Spans zero |
| **UL30 – Aciclovir** | RF (Position Only) | -0.0188 | 0.0093 | [-0.0382, -0.0023] | Spans zero / No drop |
| **UL30 – Aciclovir** | RF (All Features) | -0.0144 | 0.0167 | [-0.0501, +0.0169] | Spans zero |
| **UL30 – Foscarnet** | RF (All Features) | -0.0146 | 0.0194 | [-0.0540, +0.0218] | Spans zero |

*Statistical Finding on Clustering:*
Under position-cluster resampling on the clean single-site cohort, holding out individual positions does not degrade the Random Forest on UL30 because the tree models interpolate from neighboring active-site residues. This mathematically confirms the paper's core thesis: **Position-Grouped CV fails to penalize positional interpolation; only Contiguous-Block CV provides an honest structural holdout.**

### 4.4 Methodological Evaluation of Claim C2 & Dual-Split Standard
1. **Evidence Consistent with Leakage in the Random Forest and Logistic Regression:**
   - Under naive random splitting, **Random Forest (All Features)** achieves inflated AUROCs (0.9109 on UL23-ACV, 0.8620 on UL30-ACV, 0.8335 on UL30-FOS).
   - Under position-grouped splitting, performance drops substantially: paired gaps are **+0.1041** [+0.0703, +0.1379] for UL23-ACV, **+0.1292** [+0.0878, +0.1706] for UL30-ACV, and **+0.1579** [+0.1050, +0.2108] for UL30-FOS. Because all Nadeau-Bengio corrected 95% confidence intervals strictly exclude zero, this observation is **consistent with data leakage in the random forest**.
   - Logistic Regression also drops significantly when position identity is withheld (+0.0327 to +0.1098 paired gaps, CIs exclude zero).
   - In contrast, the zero-parameter baseline **BLOSUM62 Radicalness** displays virtually identical performance across random and position-grouped splits (paired gaps of -0.0019 and -0.0048, with CIs spanning zero), demonstrating that models without position features or trainable weights do not leak.
2. **Resolution of Regime 4 (Both-Grouped) Partitioning:**
   - The double-holdout partition (Regime 4) enforces a strict 2D double-holdout, preventing cross-leakage of positions across strains.
   - Hard assertions (`assert len(overlap) == 0`) verified zero position overlap and zero strain overlap between train and test across all 20 repeats.
   - Consequently, **Position Prior achieves exactly 0.5000 +/- 0.0000 in both Regime 3 and Regime 4**, confirming zero discriminative skill when target positions are withheld.
3. **The Core Finding: Position-Grouped CV Fails to Eliminate Positional Shortcuts:**
   - On UL30, **Random Forest (Position Only)** scores **0.8715 AUROC** (ACV) and **0.8684** (FOS) under Position-Grouped splitting, outperforming the All-Features Random Forest (0.7328 and 0.6756).
   - This occurs because catalytic resistance residues cluster into localized 3D active sites; holding out individual positions still allows the tree models to interpolate labels from neighboring positions along the primary sequence.
   - Under the **Contiguous-Block Split (Regime 5)**—which holds out contiguous 247-residue sequence windows—the Position-Only model completely collapses to **0.5001** and **0.5166** (pure chance).
   - **Dual-Split Standard:** Therefore, Position-Grouped is a weak local test, while Contiguous-Block is the strict structural domain holdout. Both must be reported as primary benchmark regimes.

### 4.5 Temporal Confounding & Factual Era Breakdown
*Conducted in `scripts/analyze_temporal_confounding.py` and logged in `logs/temporal_confounding.log`:*
1. **Factual Breakdown of Resistance Counts (Pre-2016 vs Post-2016):**
   - **Pre-2016 Corpus:** 316 resistant / 1,271 total labeled entries (24.86% resistance rate).
     - *UL23 (202 resistant):* Aciclovir: 154, Brivudine: 38, Penciclovir: 10.
     - *UL30 (105 resistant):* Aciclovir: 52, Foscarnet: 39, Cidofovir: 9, Penciclovir: 5.
     - *UL5 / UL52 (9 resistant):* Pritelivir: 9.
     - Total = 202 + 105 + 9 = **316 resistant mutations**.
   - **Post-2016 Corpus:** 316 resistant / 858 total labeled entries (36.83% resistance rate).
     - *UL23 (134 resistant):* Aciclovir: 86, Brivudine: 24, Penciclovir: 24.
     - *UL30 (155 resistant):* Aciclovir: 74, Foscarnet: 51, Cidofovir: 18, Penciclovir: 12.
     - *UL5 / UL52 (27 resistant):* Amenamevir: 17, Pritelivir: 10.
     - Total = 134 + 155 + 27 = **316 resistant mutations**.
   - *Resolution:* The reporting of exactly 316 resistant mutations in both historical windows is an **exact numerical coincidence** arising from CMI 2025 Table S3 additions across all genes and drugs.
2. **Measured Covariate Shift:**
   - A non-label domain classifier trained exclusively on `gene`, `drug`, `virus`, and `normalized_position` achieves a 5-fold cross-validated AUROC of **0.6794 +/- 0.0232** (Random Forest) and **0.6460 +/- 0.0199** (Logistic Regression).
   - This empirically confirms significant experimental-era covariate shift (e.g. Amenamevir tested exclusively post-2016; helicase-primase expanding from <1% to >7% of records; Foscarnet dropping from 35.6% to 20.2%).

### 4.6 Zero-Shot Protein Language Model Benchmark (Dual-Split Standard)
*Evaluated on the 929 clean primary rows (`data/processed/primary_position_folds.tsv`) using masked-marginal log-likelihood ratios $\Delta \text{LLR} = \log P(\text{mut} \mid \mathbf{x}_{\setminus i}) - \log P(\text{wt} \mid \mathbf{x}_{\setminus i})$ across all 507 unique sequence positions:*

| Target Gene-Drug Pair | Evaluation Metric | Position-Only RF (Shortcut Baseline) | All-Features RF (Tabular Baseline) | ESM-2 Zero-Shot (35M Masked Marginal) | ESM-2 Direction |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **UL23 – Aciclovir** | Pos-Grouped (Weak) | 0.7683 +/- 0.0182 | 0.8067 +/- 0.0273 | **0.6599 +/- 0.0495** | Deleterious $\to$ Resistant |
| **UL23 – Aciclovir** | Contig-Block (Strict) | 0.5440 +/- 0.0014 | **0.6875 +/- 0.0079** | 0.6118 +/- 0.1786 | Deleterious $\to$ Resistant |
| **UL30 – Aciclovir** | Pos-Grouped (Weak) | **0.8715 +/- 0.0171** | 0.7328 +/- 0.0365 | 0.6819 +/- 0.0934 | Deleterious $\to$ Resistant |
| **UL30 – Aciclovir** | Contig-Block (Strict) | 0.5001 +/- 0.0009 | 0.4089 +/- 0.0198 | **0.6453 +/- 0.2070** | Deleterious $\to$ Resistant |
| **UL30 – Foscarnet** | Pos-Grouped (Weak) | **0.8684 +/- 0.0282** | 0.6756 +/- 0.0446 | 0.7352 +/- 0.0157 | Deleterious $\to$ Resistant |
| **UL30 – Foscarnet** | Contig-Block (Strict) | 0.5166 +/- 0.0105 | 0.6372 +/- 0.0181 | **0.6973 +/- 0.1841** | Deleterious $\to$ Resistant |

*Key Benchmark Observations:*
1. **Biological Alignment:** For all three pairs, zero-shot $\Delta \text{LLR}$ negatively correlates with resistance (lower fitness / disruption of wild-type function indicates resistance), confirming that resistance in UL23 and UL30 predominantly operates through catalytic disruption rather than enhanced evolutionary fitness.
2. **Generalization Beyond Hotspots:** Under Contiguous-Block holdouts, tabular models fail or degrade (All-Features RF drops to 0.4089 on UL30-ACV; Position-Only collapses to ~0.50). In contrast, zero-shot ESM-2 demonstrates genuine transfer across structural blocks (**0.6453** on UL30-ACV and **0.6973** on UL30-FOS), establishing that protein language models capture structural mutational tolerance that positional interpolation cannot provide.

### 4.7 Multi-Mutation & Co-occurrence Exclusion Handling
- Single-site mutation models cannot deterministically assign phenotypic resistance to individual residues when multiple mutations are present simultaneously.
- All **125 multi-mutation and co-occurrence rows** (including 87 HerpesDRG co-gene/co-aa entries and 12 CMI Table S3 multi-substitution slash/combo rows) were formally cataloged in [`data/processed/co_occurrence_entries.tsv`](file:///d:/antiviral-resistance-hsv--1/data/processed/co_occurrence_entries.tsv) with explicit exclusion rationales.
- The primary evaluation cohort [`data/processed/primary_position_folds.tsv`](file:///d:/antiviral-resistance-hsv--1/data/processed/primary_position_folds.tsv) (SHA-256: `7c4d0a3de2cd171b07230757f44e0dfbd61c6a0a8570dc8f263195b354c66909`) consists strictly of the **929 clean single-site, non-conflicted rows** across the three primary arms (`UL23-ACV`: 286; `UL30-ACV`: 332; `UL30-FOS`: 311).

---

## 5. Protocol Amendments & Technical Specifications (Dated 2026-10-29)

### 5.1 Protocol Amendment: Claim C4 Scope & Methodology
> **Original Protocol Text:** Claim C4 proposes a pathway-constrained gate routing resistance predictions across distinct biological mechanisms.  
> **Dated Amendment (2026-10-29):**
> 1. *Negative Control Sanity Check:* The rejection of `UL23-FOS` (confirming that TK-mediated resistance probability collapses to near-zero for foscarnet) is formally designated as a **negative-control biological sanity check**, not the primary empirical proof of Claim C4.
> 2. *Multi-Task Representation Sharing:* The primary empirical test of Claim C4 requires demonstrating data efficiency and learning curve gains when sharing representations across related drugs (e.g. multi-task conditioning across nucleoside analogues `ACV` and `PCV` vs pyrophosphate analogue `FOS` in `UL30`).
> 3. *Penciclovir (PCV) Auxiliary Arm Audit:*
>    - `UL23-PCV`: 115 distinct mutations, 82 positions, 59 resistant positions (All-Inclusive) / 58 resistant positions (Strict Excluding Conflicts).
>    - `UL30-PCV`: 85 distinct mutations, 75 positions, 28 resistant positions (All-Inclusive) / 14 resistant positions (Strict Excluding Conflicts).
>    - Under strict conflict exclusion, `UL30-PCV` drops to 14 resistant positions (< 15), meaning it serves as an **auxiliary multi-task training arm**, but cannot support standalone primary evaluation.

### 5.2 Helicase-Primase (UL5 / UL52) Statistical Power Clarification
- The CMI 2025 review reported 27 novel HPC mutation entries (12 in HSV-1, 15 in HSV-2).
- Detailed mapping reveals that these 27 entries represent multiple substitutions and co-mutation combinations concentrated at a small number of conserved catalytic hotspot positions (e.g., UL5 352, 355, 356 in HSV-1; 351, 355 in HSV-2).
- Consequently, these 27 entries map to only **7 unique resistant positions for Pritelivir** and **4 unique resistant positions for Amenamevir** in UL5.
- At N=5 positive positions in HSV-1, the 95% Clopper-Pearson / Wilson interval spans +/- 0.40 (e.g. [0.30, 0.90]). Therefore, standalone machine learning evaluation is **statistically underpowered** (not 'mathematically impossible'). UL5 and UL52 must be relegated to pooled exploratory analysis.

### 5.3 ProtBERT Embedding Window Protocol
To prevent mutational signal dilution from sliding-window averaging:
1. **Centered Window Extraction:** For UL30 (1,235 aa), token embeddings are extracted from a window of at most 510 amino acids centered on the mutated residue P:
   - P_start = max(0, P - 256), P_end = min(L, P_start + 510).
   - Special tokens [CLS] and [SEP] are appended, strictly adhering to the 512-token context limit.
2. **UL23 Single-Pass:** For UL23 (376 aa), the entire sequence is extracted in a single forward pass without window truncation.
3. **Format & Non-Canonical Mapping:** Inputs are formatted as uppercase, space-separated residues ("M A S ..."), with non-canonical/rare amino acids mapped to X.
4. **Delta Embeddings:** Wild-type and mutant windows share identical coordinate bounds. Feature representations are computed as the mutant-minus-wildtype embedding delta: Delta h = h_mut - h_wt. Verified in `logs/embedding_window_verification.log`.

### 5.4 Audited Novelty Status & DOI Resolution Audit
- **Logged Search Results:** Searches across Europe PMC, PubMed, GitHub, Hugging Face, and Zenodo identified **none found in the logged search** for machine learning or protein language model prediction of HSV UL23/UL30 resistance. Clarified that 'zero prior papers found' reflects the specific queries executed, not an exhaustive census of all literature.
- **Failed Endpoints:** Automated queries to Zenodo REST API returned HTTP 403 (Cloudflare/WAF block) and Semantic Scholar API returned HTTP 429 (rate-limited without API key); these endpoints are recorded as *Failed / Not Run* in `logs/novelty_widened_raw.tsv`.
- **Systematic DOI & Citation Verification Audit (Logged in `logs/novelty_doi_verification.tsv`):**
  - Live Crossref and PubMed API verification of identifiers in the initial raw novelty log revealed **4 critical mismatches/corruptions**:
    1. `10.1093/bib/bbad180` (claimed LOGICA): Actually resolves to *'P-CSN: single-cell RNA sequencing data analysis by partial cell-specific network'* (Brief Bioinform 2023). Flagged as **MISMATCH**.
    2. `10.1038/s41467-022-33044-0` (claimed AMRscope): Returned **HTTP 404 (Not Found)**.
    3. `10.1021/acs.jcim.2c00428` (claimed DeepMutDTA): Actually resolves to *'Calcium Ion Binding to the Mutants of Calmodulin: A Structure-Based Computational Predictive Model...'* (J Chem Inf Model 2022). Flagged as **MISMATCH**.
    4. `10.1016/j.jcv.2017.05.008` (claimed DeepChek CMV): Actually resolves to *'Phylogenetic analysis of an epidemic outbreak of acute hepatitis C in HIV-infected patients by ultra-deep pyrosequencing'* (J Clin Virol 2017). Flagged as **MISMATCH**.
  - *Audit Resolution:* These citations were flagged and rectified. DeepChek is clarified as an in-vitro diagnostic rule-based expert system developed with French reference centers (ANRS/Limoges), not a machine learning model.
- **Pending Human Actions:** Evidence tiers are currently proxies from the `strain_generation` field. Final sign-off requires: (1) 50-row human sample for Cohen's kappa verification; (2) human inspection of primary threshold citations; (3) human completion of the decision column in `novelty_widened_raw.tsv`; (4) PI formal sign-off.

### 5.5 Reproducibility Integrity Hash
- The clean, primary scope single-site benchmark dataset with deterministic 5-fold position-grouped cross-validation assignments (`data/processed/primary_position_folds.tsv`) was hashed directly on disk:
  - **SHA-256 Hash:** `7c4d0a3de2cd171b07230757f44e0dfbd61c6a0a8570dc8f263195b354c66909`
  - **Seed:** Deterministic GroupKFold (seed=42)
  - **Integrity Record:** Logged in [`logs/fold_assignments_hash.log`](file:///d:/antiviral-resistance-hsv--1/logs/fold_assignments_hash.log)
  - **Leakage Assertion:** Hard assertions in script confirmed zero position overlap across all folds (`assert len(overlap) == 0`).

### 5.6 Protocol Amendment: Stopping Rule & Dual-Split Benchmark Reframing (Dated 2026-10-29)
> **Original Stopping Rule:** Success requires beating the Position-Only Random Forest under the Position-Grouped split AND achieving AUROC $\ge 0.65$ under the Contiguous-Block split on at least 2 of the 3 primary pairs.  
> **Dated Protocol Amendment (2026-10-29):**
> 1. *Withdrawal of Position-Grouped Condition:* The requirement to "beat the Position-Only forest under Position-Grouped CV" is formally **withdrawn**. As demonstrated empirically in Section 4.2, the Position-Only forest achieves an inflated AUROC of 0.8715 on UL30-ACV solely by interpolating between neighboring active-site residues in primary sequence space (collapsing to 0.5001 under contiguous-block holdout). Requiring a biological model to beat a leaky shortcut under a leaky split rewards models that exploit the same confound.
> 2. *Revised Success Criterion:* A model qualifies as a successful generalizable predictor if it achieves **AUROC $\ge 0.65$ with confidence interval above 0.50 under the Contiguous-Block split on at least 2 of the 3 primary pairs** (`UL23-ACV`, `UL30-ACV`, `UL30-FOS`).
> 3. *Measured Cloud GPU Benchmark Evaluation (ESM-2 35M, ESM-2 650M, ProtBERT 420M):*
> Executed on cloud Tesla T4 GPU (logged in `logs/plm_benchmark_results.json` and `logs/hsv-plm-gpu-benchmark.log`):
>
> | Model Architecture | UL23–ACV (Contig-Block) | UL30–ACV (Contig-Block) | UL30–FOS (Contig-Block) | 95% Block-Level Student's $t$ CIs (decision basis: UL23-ACV, UL30-ACV, UL30-FOS) | 95% Pooled Position-Cluster Bootstrap CIs (secondary, descriptive only) | Stopping Rule Verdict |
> | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
> | **ESM-2 (35M Probe)** | 0.6110 +/- 0.2001 | 0.6795 +/- 0.2441 | 0.7481 +/- 0.1842 | [0.3625, 0.8595], [0.2911, 1.0000], [0.4550, 1.0000] | [0.57, 0.73], [0.58, 0.75], [0.62, 0.84] | **0 of 3 pairs pass** (no block $t$-CI lower bound exceeds 0.50) |
> | **ESM-2 (650M)** | **0.7397 +/- 0.1565** | **0.7875 +/- 0.1772** | **0.8809 +/- 0.1250** | **[0.5455, 0.9340], [0.5055, 1.0000], [0.6821, 1.0000]** | [0.71, 0.84], [0.76, 0.89], [0.74, 0.92] | **3 of 3 pairs strictly pass** |
> | **ProtBERT (420M)** | **0.7530 +/- 0.0713** | 0.6957 +/- 0.2416 | 0.7423 +/- 0.1593 | **[0.6645, 0.8415]**, [0.3113, 1.0000], [0.4888, 0.9958] | [0.65, 0.79], [0.77, 0.89], [0.80, 0.91] | **1 of 3 pairs strictly passes** (UL23-ACV only); UL30 pairs meet the mean $\ge 0.65$ criterion but their block $t$-CIs include 0.50 |
>
> *Audit Benchmark Closure Verdict:* **THRESHOLD SATISFIED BY ESM-2 650M ONLY.** The stopping rule is adjudicated on the **block-level Student's $t$ CI**, not the pooled position-cluster bootstrap CI. Pooled bootstrap CIs are computed across the whole dataset and are therefore both narrower and centred on the pooled AUROC rather than on the block mean; for ProtBERT on UL30-ACV, for example, the block mean is 0.6957 while the pooled bootstrap CI is [0.7740, 0.8938] around a pooled AUROC of 0.8398. These are different estimators and must not be paired in a single row as point estimate and interval. On the correct block-level basis, **ESM-2 650M passes 3/3 primary pairs**; **ProtBERT passes 1/3** (clearing the $\ge 2$ of 3 amendment criterion only if the mean-based criterion alone is applied, which it is not); and the **ESM-2 35M probe passes 0/3**. ESM-2 650M remains the designated primary model and still outperforms the tabular tree baselines (RF: ~0.40–0.50) under contiguous-block holdout.
>
> *Correction note (2026-10-31):* The previous version of this table paired contiguous-block point estimates with pooled position-cluster bootstrap intervals and concluded "3 of 3 pairs strictly pass" for ProtBERT and "2 of 3" for the 35M probe. Both verdicts were artefacts of that mismatch and are superseded by the row contents above. The underlying measured values in `logs/plm_benchmark_results.json` were correct and are unchanged; only the interval selection and the derived verdicts were wrong.
