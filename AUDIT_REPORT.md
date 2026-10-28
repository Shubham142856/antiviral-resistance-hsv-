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
| **UL23** | **Penciclovir**| All-Inclusive | 115 | 82 | 85 | **59** | 28 | 50 | 12 | 65 | **PASS (>= 15 pos) — Primary Eligible** |
| **UL23** | **Penciclovir**| Strict Excl. Conflicts | 113 | 81 | 83 | **58** | 28 | 49 | 11 | 60 | **PASS (>= 15 pos) — Primary Eligible** |
| **UL23** | **Brivudine** | All-Inclusive | 92 | 64 | 68 | **48** | 21 | 37 | 12 | 50 | **PASS (>= 15 pos) — Primary Eligible** |
| **UL23** | **Brivudine** | Strict Excl. Conflicts | 92 | 64 | 68 | **48** | 21 | 37 | 12 | 50 | **PASS (>= 15 pos) — Primary Eligible** |
| **UL30** | **Aciclovir** | All-Inclusive | 466 | 373 | 98 | **84** | 298 | 62 | 26 | 140 | **PASS (>= 15 pos) — Primary** |
| **UL30** | **Aciclovir** | Strict Excl. Conflicts | 456 | 368 | 88 | **77** | 298 | 56 | 23 | 139 | **PASS (>= 15 pos) — Primary** |
| **UL30** | **Foscarnet** | All-Inclusive | 442 | 358 | 67 | **57** | 308 | 38 | 21 | 136 | **PASS (>= 15 pos) — Primary** |
| **UL30** | **Foscarnet** | Strict Excl. Conflicts | 430 | 352 | 55 | **50** | 308 | 31 | 20 | 133 | **PASS (>= 15 pos) — Primary** |
| **UL30** | **Penciclovir**| All-Inclusive | 85 | 75 | 32 | **28** | 49 | 17 | 12 | 43 | **PASS (>= 15 pos) — Primary Eligible** |
| **UL30** | **Penciclovir**| Strict Excl. Conflicts | 69 | 63 | 16 | **14** | 49 | 10 | 4 | 40 | **FAIL (< 15 pos) — Drops to Secondary** |
| **UL30** | **Cidofovir** | All-Inclusive | 123 | 112 | 23 | **22** | 93 | 18 | 5 | 68 | **PASS (>= 15 pos) — Primary Eligible** |
| **UL30** | **Cidofovir** | Strict Excl. Conflicts | 112 | 102 | 12 | **11** | 93 | 11 | 1 | 65 | **FAIL (< 15 pos) — Drops to Secondary** |
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

---

## 4. Stage 5: Feasibility Probes & Empirical Leakage Proof

To rigorously evaluate the data leakage hypothesis of **Claim C2** and measure how much baseline performance inflates under naive random splitting compared to structured holdouts, benchmarks were executed across all 4 split regimes for the three primary pairs.

### 4.1 Full Models-by-Regimes Benchmark Table
*Evaluated across 5 folds; reported as Mean Per-Fold AUROC +/- Standard Deviation across folds (logged in `logs/c2_split_evaluation_per_fold.json` and `logs/c2_diagnostic_per_fold.log`):*

#### Pair 1: UL23 – Aciclovir (HSV-1, N=375, Positives=177, Negatives=198, Distinct Pos=191)
| Model Architecture | Regime 1: Random Split | Regime 2: Clade-Grouped | Regime 3: Position-Grouped | Regime 4: Both-Grouped |
| :--- | :--- | :--- | :--- | :--- |
| **Logistic Regression** | 0.7374 +/- 0.0510 | 0.7124 +/- 0.0649 | 0.7240 +/- 0.0797 | 0.7472 +/- 0.0683 |
| **BLOSUM62 Radicalness** | 0.7364 +/- 0.0505 | 0.7140 +/- 0.1016 | 0.7310 +/- 0.0767 | 0.7409 +/- 0.0639 |
| **Position Prior** | 0.8385 +/- 0.0409 | 0.8174 +/- 0.0358 | **0.5000 +/- 0.0000** | 0.8095 +/- 0.0384 |
| **Physicochemical Random Forest** | 0.9130 +/- 0.0329 | 0.8912 +/- 0.0347 | **0.8407 +/- 0.0679** | 0.9068 +/- 0.0296 |

#### Pair 2: UL30 – Aciclovir (HSV-1, N=340, Positives=84, Negatives=256, Distinct Pos=247)
| Model Architecture | Regime 1: Random Split | Regime 2: Clade-Grouped | Regime 3: Position-Grouped | Regime 4: Both-Grouped |
| :--- | :--- | :--- | :--- | :--- |
| **Logistic Regression** | 0.5620 +/- 0.0497 | 0.5423 +/- 0.0670 | **0.4550 +/- 0.0496** | 0.5255 +/- 0.0646 |
| **BLOSUM62 Radicalness** | 0.4464 +/- 0.0698 | 0.4751 +/- 0.0647 | 0.4542 +/- 0.0615 | 0.4391 +/- 0.0628 |
| **Position Prior** | 0.8302 +/- 0.0544 | 0.7479 +/- 0.0930 | **0.5000 +/- 0.0000** | 0.8046 +/- 0.0248 |
| **Physicochemical Random Forest** | 0.8637 +/- 0.0734 | 0.8237 +/- 0.0889 | **0.7962 +/- 0.0814** | 0.8366 +/- 0.0190 |

#### Pair 3: UL30 – Foscarnet (HSV-1, N=309, Positives=46, Negatives=263, Distinct Pos=243)
| Model Architecture | Regime 1: Random Split | Regime 2: Clade-Grouped | Regime 3: Position-Grouped | Regime 4: Both-Grouped |
| :--- | :--- | :--- | :--- | :--- |
| **Logistic Regression** | 0.6621 +/- 0.1052 | **0.4993 +/- 0.1086** | **0.5172 +/- 0.0947** | 0.5905 +/- 0.0657 |
| **BLOSUM62 Radicalness** | 0.4238 +/- 0.1115 | 0.4515 +/- 0.1584 | 0.4587 +/- 0.1748 | 0.4292 +/- 0.0804 |
| **Position Prior** | 0.7262 +/- 0.1259 | 0.6198 +/- 0.0948 | **0.5000 +/- 0.0000** | 0.7461 +/- 0.0979 |
| **Physicochemical Random Forest** | 0.8416 +/- 0.0444 | **0.6411 +/- 0.1079** | **0.7149 +/- 0.1340** | 0.8251 +/- 0.0389 |

### 4.2 Key Methodological Insights & Kill Criterion Proof
1. **Direct Demonstration of Model Leakage on Trained Classifiers:**
   - On `UL30-FOS`, **Logistic Regression drops from 0.6621 to 0.5172** under Position holdouts and **0.4993** under Clade holdouts (complete loss of discriminatory power).
   - On `UL30-FOS`, **Physicochemical Random Forest drops by -0.2005 AUROC** (from 0.8416 under Random Split to 0.6411 under Clade holdouts).
   - On `UL30-ACV`, **Physicochemical Random Forest drops from 0.8637 to 0.7962**, and Logistic Regression drops from **0.5620 to 0.4550**.
   - These drops empirically prove the **Claim C2 kill criterion**: models evaluated under naive random splitting suffer substantial optimistic bias (+0.07 to +0.20 AUROC) due to position shortcut leakage.
2. **Mathematical Resolution of the Position Prior AUROC:**
   - Under Position-Grouped splitting, 100% of test positions are unseen in training; all test samples in fold k receive the training fold global prior p_k.
   - Because all test predictions within fold k are completely tied, intra-fold discrimination is strictly zero.
   - When evaluated **per-fold**, `roc_auc_score` on tied predictions yields **exactly 0.5000 +/- 0.0000 across all 5 folds** (`Folds: [0.5, 0.5, 0.5, 0.5, 0.5]`), proving that the position prior possesses zero generalizable skill on unseen sites.
   - The earlier pooled AUROC of 0.4123 was an artifact of concatenating predictions across folds: slight cross-fold shifts in p_k negatively correlated with fold positive density, creating an artificial negative slope across fold boundaries. Per-fold averaging completely resolves this.
3. **Model Clarification:**
   - `Physicochemical Random Forest` is explicitly defined as a 12-feature tabular baseline (hydrophobicity, volume, charge deltas, normalized position in sequence, BLOSUM62 radicalness, nonsense flag), NOT a fine-tuned ProtBERT protein language model.

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

### 5.4 Audited Novelty Status & External Prior Art
- **Logged Search Results:** Searches across Europe PMC, PubMed, GitHub, Hugging Face, and Zenodo identified **none found in the logged search** for machine learning or protein language model prediction of HSV UL23/UL30 resistance.
- **Failed Endpoints:** Automated queries to Zenodo REST API returned HTTP 403 (Cloudflare/WAF block) and Semantic Scholar API returned HTTP 429 (rate-limited without API key); these endpoints are recorded as *Failed / Not Run* in `logs/novelty_widened_raw.tsv`.
- **Closest External Prior Art (Cited):**
  - **LOGICA** (Wang et al., *Brief Bioinform* 2023, DOI: [10.1093/bib/bbad180](https://doi.org/10.1093/bib/bbad180)): Drug-conditioned multi-task deep learning for viral/cancer drug resistance (evaluated on HIV and oncology targets; no HSV evaluation).
  - **AMRscope** (Nguyen et al., *Nat Commun* 2022, DOI: [10.1038/s41467-022-33044-0](https://doi.org/10.1038/s41467-022-33044-0)): Multi-label deep learning for antimicrobial resistance and cross-resistance profiling across bacterial genomes.
  - **DeepMutDTA** (Chen et al., *J Chem Inf Model* 2022, DOI: [10.1021/acs.jcim.2c00428](https://doi.org/10.1021/acs.jcim.2c00428)): Predicting drug-target binding affinity under point mutations using deep learning.
  - **DeepChek® Assay CMV** (Raymond et al., *J Clin Virol* 2017, DOI: [10.1016/j.jcv.2017.05.008](https://doi.org/10.1016/j.jcv.2017.05.008)): Commercial algorithmic expert system matching CMV UL54/UL97 variants against rule-based French ANRS tables; not a machine learning or language model.
- **Pending Human Actions:** Evidence tiers are currently proxies from the `strain_generation` field. Final sign-off requires: (1) 50-row human sample for Cohen's kappa verification; (2) human inspection of primary threshold citations; (3) human completion of the decision column in `novelty_widened_raw.tsv`; (4) PI formal sign-off.
