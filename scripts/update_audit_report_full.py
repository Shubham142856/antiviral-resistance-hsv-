import pathlib
"""
Update AUDIT_REPORT.md with:
- 20-repeat CV results with paired leakage gaps and 95% CIs
- Contiguous-block split and ablation models (All, No Pos, Pos Only)
- Fixed Regime 4 double-holdout (Position Prior = 0.5000)
- Homology alignment HSV-1 vs HSV-2 (UL23 376/376, UL30 1235/1240)
- Revert PCV, BVDU to auxiliary (pending amendment) and drop CDV
- Multi-mutation exclusions (125 rows in co_occurrence_entries.tsv, is_single_site==1)
- Temporal confounding & domain classifier results (RF AUROC 0.6794)
- Novelty DOI audit findings (4 corrupt DOIs flagged in novelty_doi_verification.tsv)
- Primary folds SHA-256 hash
"""
import re

def main():
    report_path = str(pathlib.Path(__file__).resolve().parent.parent / 'AUDIT_REPORT.md')
    with open(report_path, 'r', encoding='utf-8') as f:
        content = f.read()

    # 1. Update Table A in Section 3
    # Replace "Primary Eligible" for PCV, BVDU, CDV with "Auxiliary, Pending Amendment" / "Dropped"
    content = content.replace(
        "| **UL23** | **Penciclovir**| All-Inclusive | 115 | 82 | 85 | **59** | 28 | 50 | 12 | 65 | **PASS (>= 15 pos) — Primary Eligible** |",
        "| **UL23** | **Penciclovir**| All-Inclusive | 115 | 82 | 85 | **59** | 28 | 50 | 12 | 65 | **Auxiliary Multi-Task (Pending Amendment)** |"
    )
    content = content.replace(
        "| **UL23** | **Penciclovir**| Strict Excl. Conflicts | 113 | 81 | 83 | **58** | 28 | 49 | 11 | 60 | **PASS (>= 15 pos) — Primary Eligible** |",
        "| **UL23** | **Penciclovir**| Strict Excl. Conflicts | 113 | 81 | 83 | **58** | 28 | 49 | 11 | 60 | **Auxiliary Multi-Task (Pending Amendment)** |"
    )
    content = content.replace(
        "| **UL23** | **Brivudine** | All-Inclusive | 92 | 64 | 68 | **48** | 21 | 37 | 12 | 50 | **PASS (>= 15 pos) — Primary Eligible** |",
        "| **UL23** | **Brivudine** | All-Inclusive | 92 | 64 | 68 | **48** | 21 | 37 | 12 | 50 | **Auxiliary Multi-Task (Pending Amendment)** |"
    )
    content = content.replace(
        "| **UL23** | **Brivudine** | Strict Excl. Conflicts | 92 | 64 | 68 | **48** | 21 | 37 | 12 | 50 | **PASS (>= 15 pos) — Primary Eligible** |",
        "| **UL23** | **Brivudine** | Strict Excl. Conflicts | 92 | 64 | 68 | **48** | 21 | 37 | 12 | 50 | **Auxiliary Multi-Task (Pending Amendment)** |"
    )
    content = content.replace(
        "| **UL30** | **Penciclovir**| All-Inclusive | 85 | 75 | 32 | **28** | 49 | 17 | 12 | 43 | **PASS (>= 15 pos) — Primary Eligible** |",
        "| **UL30** | **Penciclovir**| All-Inclusive | 85 | 75 | 32 | **28** | 49 | 17 | 12 | 43 | **Auxiliary Multi-Task (Pending Amendment)** |"
    )
    content = content.replace(
        "| **UL30** | **Penciclovir**| Strict Excl. Conflicts | 69 | 63 | 16 | **14** | 49 | 10 | 4 | 40 | **FAIL (< 15 pos) — Drops to Secondary** |",
        "| **UL30** | **Penciclovir**| Strict Excl. Conflicts | 69 | 63 | 16 | **14** | 49 | 10 | 4 | 40 | **FAIL (< 15 pos) — Auxiliary / Secondary** |"
    )
    content = content.replace(
        "| **UL30** | **Cidofovir** | All-Inclusive | 123 | 112 | 23 | **22** | 93 | 18 | 5 | 68 | **PASS (>= 15 pos) — Primary Eligible** |",
        "| **UL30** | **Cidofovir** | All-Inclusive | 123 | 112 | 23 | **22** | 93 | 18 | 5 | 68 | **Excluded (Dropped per Protocol)** |"
    )
    content = content.replace(
        "| **UL30** | **Cidofovir** | Strict Excl. Conflicts | 112 | 102 | 12 | **11** | 93 | 11 | 1 | 65 | **FAIL (< 15 pos) — Drops to Secondary** |",
        "| **UL30** | **Cidofovir** | Strict Excl. Conflicts | 112 | 102 | 12 | **11** | 93 | 11 | 1 | 65 | **FAIL (< 15 pos) — Excluded (Dropped)** |"
    )

    # 2. Add Homology Mapping Subsection under Section 3
    homology_section = """
### Table F: Homology Alignment & Sequence Mapping (HSV-1 vs HSV-2)
*Global pairwise alignment between canonical reference strains: HSV-1 (NC_001806.2 / 17) and HSV-2 (NC_001798.2 / HG52) logged in `logs/homology_mapping.log` and `data/processed/hsv1_hsv2_homology_map.tsv`:*

- **UL23 (Thymidine Kinase):** Both proteins are exactly 376 amino acids long. 277 residues are identical (73.67%), 97 substitutions, 2 gapped alignment columns.
- **UL30 (DNA Polymerase):** HSV-1 is 1,235 aa; HSV-2 is 1,240 aa (+5 aa net difference). 1,123 residues are identical (90.93%), 107 substitutions, 5 gapped alignment columns (insertions/deletions concentrated in variable linker regions around residues 675–690 and 1106–1119).
- *Methodological Consequence:* Pooling raw sequence positions across HSV-1 and HSV-2 directly introduces index mismatching for UL30. Therefore:
  1. The **Primary Analysis Scope is strictly restricted to HSV-1**.
  2. Homologous cross-species pooling requires canonical `homology_cluster_id` mapping.
  3. HSV-2 is reserved as an **external cross-species generalization test set**.
"""
    if "### Table F: Homology Alignment" not in content:
        # insert before Section 4
        content = content.replace("--- \n\n## 4. Stage 5:", homology_section + "\n---\n\n## 4. Stage 5:")
        content = content.replace("---\n\n## 4. Stage 5:", homology_section + "\n---\n\n## 4. Stage 5:")

    # 3. Replace Section 4 with 20-repeat CV results
    new_stage5 = """## 4. Stage 5: Feasibility Probes & Evaluation of Leakage Bias (Claim C2)

To rigorously test the data leakage hypothesis of **Claim C2** and measure the performance inflation resulting from naive random splitting compared to structured holdouts, repeated cross-validation experiments were conducted across **20 independent repeats** (5 folds per repeat = 100 evaluations per cell, random seeds 1..20).

### 4.1 Repeated CV Benchmark Table across Split Regimes (20 Repeats, Mean Per-Fold AUROC +/- SD)
*Evaluated across 6 ablation models and 5 split regimes (logged in `logs/c2_repeated_cv.log` and `logs/c2_repeated_cv_results.json`). Paired leakage gaps represent the mean paired difference $(\\Delta = \\text{Naive Random} - \\text{Position-Grouped})$ with 95% bootstrap confidence intervals across repeats:*

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
"""

    # Replace Section 4 in content
    # Find start of Section 4 and start of Section 5
    s4_idx = content.find("## 4. Stage 5:")
    s5_idx = content.find("## 5. Protocol Amendments")
    if s4_idx != -1 and s5_idx != -1:
        content = content[:s4_idx] + new_stage5 + "\n---\n\n" + content[s5_idx:]

    # 4. Update Section 5.4 with DOI Verification Findings & Hash
    new_s54 = """### 5.4 Audited Novelty Status & DOI Resolution Audit
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
"""
    # Replace Section 5.4 to end of document
    s54_idx = content.find("### 5.4 Audited Novelty Status")
    if s54_idx != -1:
        content = content[:s54_idx] + new_s54

    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(content)
        
    print(f"Successfully updated {report_path}!")

if __name__ == '__main__':
    main()
