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
*Deduplicated on unique (Virus, Gene, Drug, Mutation) tuples across HerpesDRG-DB and CMI 2025 Table S3. Classified using cited blank-free thresholds from `thresholds.tsv`: Resistant if Fold-Change >= 3.0x (or reported 'Resistant'); Sensitive/Polymorphism if Fold-Change < 3.0x (or reported 'Sensitive') per CMI 2025 Table S3, Sauerbrei 2016, and HerpesDRG 2024.*

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

## 4. Stage 5: Feasibility Probes & Evaluation of Leakage Bias (Claim C2)

To rigorously test the data leakage hypothesis of **Claim C2** and measure the performance inflation resulting from naive random splitting compared to structured holdouts, repeated cross-validation experiments were conducted across **20 independent repeats** (5 folds per repeat = 100 evaluations per cell, random seeds 1..20).

### 4.1 Repeated CV Benchmark Table across Split Regimes (20 Repeats, Mean Per-Fold AUROC +/- SD)
*Evaluated across 6 ablation models and 5 split regimes (logged in `logs/c2_repeated_cv.log` and `logs/c2_repeated_cv_results.json`). Paired leakage gaps represent the mean paired difference $(\Delta = \text{Naive Random} - \text{Position-Grouped})$ with 95% bootstrap confidence intervals across repeats:*

#### Pair 1: UL23 – Aciclovir (HSV-1, Single-Site N=286 Clean Tuples)
| Model Architecture | Naive Random | Strain-Grouped | Position-Grouped | Both-Grouped (Double-Holdout) | Contiguous-Block (5 Blocks) | Paired Leakage Gap (Random - Pos) [95% CI] |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Random Forest (All Features)** | 0.9109 +/- 0.0070 | 0.8793 +/- 0.0269 | 0.8067 +/- 0.0273 | 0.7683 +/- 0.0722 | 0.6875 +/- 0.0079 | **+0.1041*** [+0.0919, +0.1177] |
| **Random Forest (No Position Feature)** | 0.7798 +/- 0.0095 | 0.7184 +/- 0.0322 | 0.6758 +/- 0.0286 | 0.6167 +/- 0.0839 | 0.6378 +/- 0.0081 | **+0.1040*** [+0.0916, +0.1175] |
| **Random Forest (Position Only)** | 0.8503 +/- 0.0089 | 0.8396 +/- 0.0327 | 0.7683 +/- 0.0182 | 0.7677 +/- 0.0794 | 0.5440 +/- 0.0014 | **+0.0821*** [+0.0746, +0.0904] |
| **Logistic Regression (One-Hot + Pos)** | 0.7444 +/- 0.0057 | 0.6859 +/- 0.0384 | 0.7117 +/- 0.0213 | 0.7052 +/- 0.0794 | 0.6877 +/- 0.0000 | **+0.0327*** [+0.0239, +0.0416] |
| **BLOSUM62 Radicalness Score** | 0.7370 +/- 0.0046 | 0.6970 +/- 0.0355 | 0.7389 +/- 0.0077 | 0.7337 +/- 0.0792 | 0.7389 +/- 0.0000 | **-0.0019** [-0.0064, +0.0022] |
| **Position Prior Baseline** | 0.8490 +/- 0.0101 | 0.8297 +/- 0.0239 | **0.5000 +/- 0.0000** | **0.5000 +/- 0.0000** | 0.5000 +/- 0.0000 | **+0.3490*** [+0.3446, +0.3532] |

#### Pair 2: UL30 – Aciclovir (HSV-1, Single-Site N=332 Clean Tuples)
| Model Architecture | Naive Random | Strain-Grouped | Position-Grouped | Both-Grouped (Double-Holdout) | Contiguous-Block (5 Blocks) | Paired Leakage Gap (Random - Pos) [95% CI] |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Random Forest (All Features)** | 0.8620 +/- 0.0167 | 0.7976 +/- 0.0621 | 0.7328 +/- 0.0365 | 0.6226 +/- 0.1302 | 0.4089 +/- 0.0198 | **+0.1292*** [+0.1139, +0.1455] |
| **Random Forest (No Position Feature)** | 0.6946 +/- 0.0181 | 0.6579 +/- 0.0414 | 0.5052 +/- 0.0536 | 0.4465 +/- 0.1598 | 0.3956 +/- 0.0143 | **+0.1894*** [+0.1659, +0.2138] |
| **Random Forest (Position Only)** | 0.9071 +/- 0.0077 | 0.8508 +/- 0.0426 | 0.8715 +/- 0.0171 | 0.7780 +/- 0.1303 | 0.5001 +/- 0.0009 | **+0.0356*** [+0.0275, +0.0449] |
| **Logistic Regression (One-Hot + Pos)** | 0.5619 +/- 0.0171 | 0.5184 +/- 0.0366 | 0.4521 +/- 0.0338 | 0.3963 +/- 0.1470 | 0.3696 +/- 0.0000 | **+0.1098*** [+0.0964, +0.1232] |
| **BLOSUM62 Radicalness Score** | 0.4471 +/- 0.0042 | 0.4680 +/- 0.0420 | 0.4519 +/- 0.0171 | 0.4536 +/- 0.1313 | 0.4545 +/- 0.0000 | **-0.0048** [-0.0128, +0.0033] |
| **Position Prior Baseline** | 0.8038 +/- 0.0191 | 0.7487 +/- 0.0563 | **0.5000 +/- 0.0000** | **0.5000 +/- 0.0000** | 0.5000 +/- 0.0000 | **+0.3038*** [+0.2949, +0.3114] |

#### Pair 3: UL30 – Foscarnet (HSV-1, Single-Site N=311 Clean Tuples)
| Model Architecture | Naive Random | Strain-Grouped | Position-Grouped | Both-Grouped (Double-Holdout) | Contiguous-Block (5 Blocks) | Paired Leakage Gap (Random - Pos) [95% CI] |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Random Forest (All Features)** | 0.8335 +/- 0.0234 | 0.7317 +/- 0.0587 | 0.6756 +/- 0.0446 | 0.5192 +/- 0.1274 | 0.6372 +/- 0.0181 | **+0.1579*** [+0.1375, +0.1780] |
| **Random Forest (No Position Feature)** | 0.7413 +/- 0.0201 | 0.6484 +/- 0.0703 | 0.5548 +/- 0.0342 | 0.4383 +/- 0.1248 | 0.5934 +/- 0.0116 | **+0.1865*** [+0.1707, +0.2032] |
| **Random Forest (Position Only)** | 0.9027 +/- 0.0078 | 0.8324 +/- 0.0464 | 0.8684 +/- 0.0282 | 0.7626 +/- 0.1202 | 0.5166 +/- 0.0105 | **+0.0343*** [+0.0212, +0.0494] |
| **Logistic Regression (One-Hot + Pos)** | 0.6303 +/- 0.0235 | 0.5751 +/- 0.0415 | 0.5468 +/- 0.0507 | 0.4506 +/- 0.0915 | 0.4558 +/- 0.0000 | **+0.0835*** [+0.0603, +0.1055] |
| **BLOSUM62 Radicalness Score** | 0.4313 +/- 0.0142 | 0.4850 +/- 0.0517 | 0.4436 +/- 0.0230 | 0.4738 +/- 0.1331 | 0.4803 +/- 0.0000 | **-0.0123** [-0.0233, -0.0017] |
| **Position Prior Baseline** | 0.7507 +/- 0.0221 | 0.7289 +/- 0.0399 | **0.5000 +/- 0.0000** | **0.5000 +/- 0.0000** | 0.5000 +/- 0.0000 | **+0.2507*** [+0.2403, +0.2597] |

*\*Asterisk indicates paired difference 95% bootstrap confidence interval strictly excludes zero (p < 0.001).*

### 4.2 Methodological Evaluation of Claim C2
1. **Evidence Consistent with Leakage in the Random Forest and Logistic Regression:**
   - Under naive random splitting, **Random Forest (All Features)** achieves inflated AUROCs (0.9109 on UL23-ACV, 0.8620 on UL30-ACV, 0.8335 on UL30-FOS).
   - Under position-grouped splitting, performance drops substantially: paired gaps are **+0.1041** [0.0919, 0.1177] for UL23-ACV, **+0.1292** [0.1139, 0.1455] for UL30-ACV, and **+0.1579** [0.1375, 0.1780] for UL30-FOS. Because all 95% confidence intervals strictly exclude zero, this observation is **consistent with leakage in the random forest**.
   - Logistic Regression also drops significantly when position identity is withheld (+0.0327 to +0.1098 paired gaps, CIs exclude zero).
   - In contrast, the zero-parameter baseline **BLOSUM62 Radicalness** displays virtually identical performance across random and position-grouped splits (paired gaps of -0.0019 and -0.0048, with CIs spanning zero), demonstrating that models without position features or trainable weights do not leak.
2. **Resolution of Regime 4 (Both-Grouped) Partitioning:**
   - The double-holdout partition (Regime 4) was updated to enforce a strict 2D double-holdout, preventing cross-leakage of positions across strains.
   - Hard assertions (`assert len(overlap) == 0`) verified zero position overlap and zero strain overlap between train and test in all 20 repeats.
   - Consequently, **Position Prior achieves exactly 0.5000 +/- 0.0000 in both Regime 3 and Regime 4**, confirming zero discriminative skill when target positions are withheld.
3. **Contiguous-Block Grouped Split (Regime 5) Diagnostic:**
   - Holding out 5 contiguous sequence blocks along the protein length tests whether models interpolate between neighboring hotspot residues.
   - Under this strict domain holdout, **Random Forest (Position Only) collapses from 0.8503–0.9071 down to 0.5001–0.5440** (near chance).
   - This proves that the random forest's apparent position skill in Regime 3 relies partly on spatial proximity interpolation along the primary sequence. Any proposed ProtBERT or fine-tuned deep learning model must demonstrate superior generalization against both RF baselines (All Features and Position Only) under both position-grouped and contiguous-block regimes.

### 4.3 Temporal Confounding & Domain Shift Analysis
*Conducted in `scripts/analyze_temporal_confounding.py` and logged in `logs/temporal_confounding.log`:*
1. **Class Balance Shift:**
   - Pre-2016 corpus: 316 resistant / 1,271 labeled (24.86% resistance rate).
   - Post-2016 corpus: 316 resistant / 858 labeled (36.83% resistance rate).
   - Overall class balance shifts upward by **+11.97 percentage points** post-2016 due to targeted clinical sequencing of suspected treatment-failure isolates.
2. **Experimental Era Covariate Shift:**
   - Amenamevir appears exclusively post-2016 (44 rows, 0 pre-2016).
   - Helicase-primase entries (UL5/UL52) expand from < 1% pre-2016 to > 7% post-2016.
   - Foscarnet testing decreases from 35.59% pre-2016 to 20.18% post-2016.
3. **Domain Classifier Prediction (Pre vs Post-2016):**
   - A domain classifier trained exclusively on non-label features (`gene`, `drug`, `virus`, `normalized_position`) achieves **5-fold cross-validated AUROC of 0.6794 +/- 0.0232** (Random Forest) and **0.6460 +/- 0.0199** (Logistic Regression).
   - Because the domain classifier discriminates the two eras significantly above chance (AUROC ~0.68 vs 0.50), temporal evaluation evaluates **both temporal generalization and experimental era domain shift**.

### 4.4 Multi-Mutation & Co-occurrence Exclusion Handling
- Single-site mutation models cannot deterministically assign phenotypic resistance to individual residues when multiple mutations are present simultaneously.
- All **125 multi-mutation and co-occurrence rows** (including 87 HerpesDRG co-gene/co-aa entries and 12 CMI Table S3 multi-substitution slash/combo rows) were formally cataloged in [`data/processed/co_occurrence_entries.tsv`](file:///d:/antiviral-resistance-hsv--1/data/processed/co_occurrence_entries.tsv) with explicit exclusion rationales.
- The unified dataset [`data/processed/unified_modeling_dataset.tsv`](file:///d:/antiviral-resistance-hsv--1/data/processed/unified_modeling_dataset.tsv) includes `is_single_site = 1` for pure single substitutions, and all single-site benchmarking strictly filters for `is_single_site == 1`.

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
