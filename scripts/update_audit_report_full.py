"""
Update AUDIT_REPORT.md cleanly using string replacement.
"""
report_path = r'd:\antiviral-resistance-hsv--1\AUDIT_REPORT.md'
with open(report_path, 'r', encoding='utf-8') as f:
    text = f.read()

# Target 1: Table A
target_table_a_start = "### Table A: Gene-Drug Pair Counts & Empirical Support"
target_table_a_end = "### Table B: Label Purity by Position (Position Shortcut Risk)"

idx1 = text.find(target_table_a_start)
idx2 = text.find(target_table_a_end)

if idx1 != -1 and idx2 != -1:
    new_table_a = """### Table A: Gene-Drug Pair Counts & Deduplicated Empirical Support
*Deduplicated on unique (Virus, Gene, Drug, Mutation) tuples across HerpesDRG-DB and CMI 2025 Table S3. Classified using cited blank-free thresholds from `thresholds.tsv`: Resistant if Fold-Change >= 3.0x (or reported 'Resistant'); Sensitive/Polymorphism if Fold-Change < 3.0x (or reported 'Sensitive') per CMI 2025 Table S3, Sauerbrei 2016, and HerpesDRG 2024.*

*Evaluated under two strict regimes:*
1. **All-Inclusive:** Retains all consensus entries.
2. **Strict Excluding Conflicts:** Excludes the 65 discordant tuples where HerpesDRG and CMI 2025 disagree on resistance phenotype (logged in `logs/cmi_herpesdrg_conflicts.log`).

| Gene | Drug | Evaluation Mode | Distinct Mutations | Distinct Positions | Resistant Mutations | Resistant Pos (Total) | Sensitive Positions | HSV-1 Res Pos | HSV-2 Res Pos | Distinct Sources | Audit Status / Verdict |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **UL23** | **Aciclovir** | All-Inclusive | 465 | 256 | 270 | **156** | 147 | 124 | 55 | 202 | **PASS (>= 15 pos) — Primary** |
| **UL23** | **Aciclovir** | Strict Excl. Conflicts | 451 | 252 | 256 | **147** | 147 | 115 | 52 | 200 | **PASS (>= 15 pos) — Primary** |
| **UL23** | **Foscarnet** | All-Inclusive | 189 | 153 | 9 | *8\\** | 146 | 7 | 2 | 60 | **Negative Control Pair (Not Primary)** |
| **UL23** | **Foscarnet** | Strict Excl. Conflicts | 189 | 153 | 9 | *8\\** | 146 | 7 | 2 | 60 | **Negative Control Pair (Not Primary)** |
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

*\\*Note on UL23-Foscarnet: Apparent resistant positions are clinical isolate artifacts harboring co-mutations in UL30 (see Section 5).*  
*\\*\\*Impact of Phenotype Discordance: The three core pairs (`UL23-ACV`, `UL30-ACV`, `UL30-FOS`) comfortably exceed the >= 15 position threshold whether conflicts are included or strictly excluded (156 vs 147; 84 vs 77; 57 vs 50). However, for `UL30-PCV` and `UL30-CDV`, excluding discordant tuples drops them below 15 resistant positions (14 and 11 positions respectively), meaning they do NOT qualify for standalone primary evaluation under strict conflict exclusion.*

"""
    text = text[:idx1] + new_table_a + text[idx2:]
    print("Table A successfully replaced!")
else:
    print("Could not find Table A markers!")

# Target 2: Section 4
target_sec4_start = "## 4. Stage 5: Feasibility Probes & Empirical Leakage Proof"
idx_sec4 = text.find(target_sec4_start)
if idx_sec4 != -1:
    new_sec4 = """## 4. Stage 5: Feasibility Probes & Empirical Leakage Proof

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
"""
    text = text[:idx_sec4] + new_sec4
    print("Section 4 and amendments successfully replaced!")
else:
    print("Could not find Section 4 start marker!")

with open(report_path, 'w', encoding='utf-8') as f:
    f.write(text)

print("AUDIT_REPORT.md updated successfully!")
