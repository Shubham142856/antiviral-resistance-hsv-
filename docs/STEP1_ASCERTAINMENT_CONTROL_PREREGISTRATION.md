# Step 1 Pre-Registration: Within-Study Ascertainment Control and Source Diagnostics

- **Protocol Version:** 1.0 (Audit-Grade Frozen Specification)
- **Date:** November 2, 2026
- **Status:** FROZEN PRIOR TO BENCHMARK EXECUTION
- **Targets Evaluated:** HSV-1 DNA Polymerase (UL30; Aciclovir, Foscarnet) and Thymidine Kinase (UL23; Aciclovir)
- **Primary Data Source:** Deduplicated benchmark dataset (`data/processed/where_what_features.tsv`) with literature citations from `data/raw/HerpesDRG-DB` and `data/raw/CMI_2025_review`.

---

## 1. Scientific Rationale & Scope

Retrospective curation databases suffer from a fundamental ascertainment vulnerability: broad surveillance studies often sequence viral genomes and deposit uncharacterized polymorphisms (labeled sensitive), while specialized clinical resistance studies specifically test suspected failure isolates and clone resistant mutations. On UL30, this produces a substantial curation-source disparity ($P(\text{resistant} \mid \text{in CMI}) = 67.5\%$ vs. $3.2\%$ outside CMI on ACV; $55.2\%$ vs. $0.4\%$ on FOS), whereas UL23 shows no such association ($41.9\%$ vs. $42.9\%$).

This pre-registration defines **Step 1: Within-Study Ascertainment Control**, a protocol to test whether the predictive signal of the zero-shot structural/evolutionary `Where` composite on UL30-ACV and the cross-target radicalness dissociation survive when comparing resistant and sensitive variants **tested within the exact same laboratory study**.

---

## 2. Pre-Registered Protocol Specification

### A. Positivity Check (Computed Before Any AUROC and Reported First)

Before computing any predictive performance metrics, each target arm must undergo a deterministic positivity check to verify that within-study pairs are supported by adequate sample size and not concentrated in a single laboratory cohort:

1. **Metrics Computed per Target:**
   - Mixed-study count ($M$): Number of distinct studies reporting $\ge 1$ resistant AND $\ge 1$ sensitive variant.
   - Total same-study Resistant-Sensitive (R-S) pair count:
     $$N_{\text{pairs}} = \sum_{s \in \text{mixed}} (|R_s| \times |S_s|)$$
   - Top 5 studies by pair count with absolute count and percentage share of total same-study pairs.
2. **Positivity Decision Rules (Prerequisite for Hypothesis Testing):**
   - Condition 1: Total same-study pairs $N_{\text{pairs}} \ge 100$.
   - Condition 2: Mixed-study count $M \ge 5$.
   - Condition 3: Concentration rule: No single study accounts for $> 50.0\%$ of total same-study pairs.
   - **Verdict:** If any of these three conditions is violated, the target arm is declared **INCONCLUSIVE due to positivity failure** and cannot serve for binary hypothesis pass/fail decisions.
3. **Target Arm Roles:**
   - **UL30-Aciclovir:** Primary target arm for the within-study survival test.
   - **UL30-Foscarnet:** Secondary only. Expected thin ($N_{\text{pairs}} < 100$ anticipated). Furthermore, UL30-ACV and UL30-FOS evaluate the identical enzyme with overlapping substitutions, rendering them non-independent.
   - **UL23-Aciclovir:** Counter-arm for the cross-gene dissociation contrast $D$.

---

### B. S0 Stratum: Primary Within-Study Evaluator

- **Cohort Definition:** Restricted strictly to mixed studies ($|R_s| \ge 1$ and $|S_s| \ge 1$). Single-class studies (all-resistant or all-sensitive) contribute zero concordant/discordant pairs to within-study comparisons and are excluded. Study ID is extracted from the `reference` / `PMID` / `DOI` column.
- **Primary Estimator (Pair-Weighted AUROC):**
  Pairs every resistant variant $i \in R_s$ with every sensitive variant $j \in S_s$ within study $s$:
  $$\text{AUROC}_{\text{pair-weighted}} = \frac{1}{\sum_{s=1}^M |R_s||S_s|} \sum_{s=1}^M \sum_{i \in R_s} \sum_{j \in S_s} \left[ \mathbb{I}(score_i > score_j) + 0.5 \cdot \mathbb{I}(score_i = score_j) \right]$$
- **Sensitivity Estimator (Study-Equal-Weighted Mean AUROC):**
  Computes the within-study AUROC $\text{AUROC}_s$ for each mixed study and averages unweighted across all $M$ mixed studies:
  $$\text{AUROC}_{\text{study-mean}} = \frac{1}{M} \sum_{s=1}^M \text{AUROC}_s$$
  This protects against quadratic pair domination by any single large mixed study.
- **Uncertainty Quantification:**
  - Study-clustered bootstrap with $B = 2000$ resamples (sampling entire studies with replacement).
  - For the cross-gene contrast $D$, a **single shared study-cluster resample** is drawn across the unified study universe to preserve inter-gene study covariance.

---

### C. Secondary Strata (Descriptive Only, No Pass/Fail)

The following strata are evaluated descriptively for contextual comparison with prior literature:
1. **S1 (CMI 2025 Cell-Culture Confirmed = Yes):** Variants directly verified by cell-culture phenotyping.
2. **S2 (HerpesDRG Numeric Fold-Change):** Variants with quantitative fold-change determinations.
3. **S3 (Recombinant / Marker Transfer + BAC):** Variants engineered into viral reference backbones. Evaluated as a single descriptive row on UL30-ACV and explicitly flagged as underpowered ($N=37$, 22 Res / 15 Sens).
- **Rule:** No S4 (lab-induced). No binary pass/fail decisions are made on S1, S2, or S3; all are reported descriptively.

---

### D. Source Diagnostics (Descriptive, No Pass/Fail, WITHIN Label Class)

Because curation source correlates strongly with phenotype on UL30, predicting source across mixed labels is circular. Ascertainment bias must be tested **within fixed phenotypic classes**:

1. **Within Sensitive Rows ($y = 0$):**
   - Test whether zero-shot scores (`Where` and `BLOSUM62`) discriminate between survey-only-study rows versus mixed-study rows (AUROC of score predicting source).
   - Test whether zero-shot scores discriminate between in-CMI rows versus not-in-CMI rows.
2. **Within Resistant Rows ($y = 1$):**
   - Compute identical source discrimination AUROCs where both source groups have $n \ge 10$ variants; otherwise report as underpowered/omitted.

---

### E. Scores Evaluated (Frozen, Zero-Shot, No Fitting)

All evaluations use frozen zero-shot metrics without any parameter fitting or retraining:
- **`Where` Score:** $W = -Z(\text{Evolutionary Constraint}) - Z(d_{\text{active-site}})$ (spec-exact formula from `docs/WHERE_WHAT_PREREGISTRATION.md`).
- **`BLOSUM62` Score:** Radicalness matrix substitution score.

---

### F. Pre-Registered Decision Rules & Outcomes (Frozen Before Any Analysis)

All outcomes are determined strictly by pre-registered criteria:

#### 1. UL30-ACV Where Survival Test:
- **PASS:** S0 pair-weighted AUROC point estimate $\ge 0.65$ AND 95% study-clustered bootstrap CI lower bound $> 0.55$.
- **ARTIFACT:** S0 pair-weighted AUROC point estimate $< 0.55$.
- **INCONCLUSIVE:** Neither PASS nor ARTIFACT criteria met (e.g., point estimate in $[0.55, 0.65)$ or CI lower bound $\le 0.55$).

#### 2. Cross-Target Dissociation Survival Test:
- Contrast metric:
  $$D = (\text{Where} - \text{BLOSUM})_{\text{UL30}} - (\text{Where} - \text{BLOSUM})_{\text{UL23}}$$
  evaluated across same-study pairs in S0.
- **PASS:** Contrast point estimate $D > 0$ AND 95% study-clustered bootstrap CI strictly excludes zero (CI lower bound $> 0$).
- **INCONCLUSIVE:** Contrast 95% CI spans zero or point estimate $D \le 0$.

---

### G. Interpretation & Residual Confound Acknowledgment

- **What PASS Shows:** Demonstrates that the predictive superiority of `Where` on UL30-ACV and the cross-target dissociation survive when variants are tested under identical laboratory protocols, eliminating study-to-study curation artifacts as the primary explanatory driver.
- **What PASS Cannot Show (Residual Confound):**
  - **Position-Selection Confound:** Authors selectively choose which sequence positions to clone and assay, heavily concentrating efforts near active-site hotspots. Observational data cannot remove this investigator selection bias.
- **What ARTIFACT Shows:** Confirms that `Where` discriminability was an artifact of pooling distal benign polymorphisms from surveillance studies against proximal resistant mutations from resistance studies.

---

### H. Single-Run Protocol Gate

A single deterministic execution script (`scripts/run_step1_ascertainment_control.py`) will execute this protocol. Output metrics will be serialized to `logs/step1_ascertainment_control_results.json`. Following execution, all results including nulls will be printed in a final numbers table, and execution will **STOP immediately** without further extensions or iterative re-analyses.

---

## 3. Pre-Run Amendment (November 2, 2026)

*This amendment was committed and pushed prior to the execution of any Step 1 benchmark script. It adds descriptive reporting diagnostics and refines interpretation without altering any pre-registered thresholds, estimators, or pass/fail decision rules.*

### Amendment Item 1: Additional Descriptive Diagnostics (No Pass/Fail)
1. **Contributing Resistant Variant Count:** Report, for each target arm, the exact number of distinct resistant variants contributing to same-study pairs. (Rationale: As observed in the positivity check, large mixed studies may have high pair counts driven by only one or two resistant variants; quantifying distinct resistant variants clarifies the effective diversity of resistant observations).
2. **Leave-One-Study-Out S0 Sensitivity:** Report the leave-one-study-out (LOSO) pair-weighted AUROC for S0, systematically excluding each mixed study one at a time. This identifies whether any single laboratory cohort exerts outsized leverage on the S0 estimate (analogous to the position 927 diagnostic).

### Amendment Item 2: Interpretation Refinement and Scope of S0 vs. S2
1. **Scope of S0 Within-Study Survival:** A PASS outcome for the Where composite on S0 demonstrates that the structural and evolutionary prior separates resistant variants from tolerated variants within the same laboratory study. However, because many sensitive variants in mixed studies are natural surveillance polymorphisms that occur at evolutionarily tolerant positions, this within-study separation cannot prove *resistance-specificity* per se (i.e. distinguishing resistant mutations from benign mutations at the same constrained sites).
2. **S2 as Nearest Resistance-Specificity Benchmark:** Stratum S2 (HerpesDRG numeric fold-change, where both resistant and sensitive variants were subjected to direct quantitative phenotypic assay) is designated as the descriptive stratum closest to true resistance-specificity. S2 results will be highlighted alongside S0 to contextualize whether positional constraint discriminates when all evaluated variants are functionally characterized.

