# Confirmatory Pre-Registration: Feature-Family × Split-Scheme Degradation Benchmark & Split-Separation Diagnostic

**Document Version:** 1.0 (Confirmatory Protocol)  
**Status:** CONFIRMATORY PRE-REGISTRATION (Strictly Zero Supervised Degradation Runs Executed)  
**Repository Anchor Commit:** Prior to execution  
**Target Proteins:** HSV-1 UL30 (DNA Polymerase, 1,235 aa) and HSV-1 UL23 (Thymidine Kinase, 376 aa)  
**Evaluated Cohorts:**
- `UL30-Aciclovir` ($N=332$, 62 R / 270 S)
- `UL30-Foscarnet` ($N=311$, 33 R / 278 S)
- `UL23-Aciclovir` ($N=286$, 120 R / 166 S)

---

## 1. Scientific Motivation & Core Novelty

Variant effect predictors in clinical virology are frequently evaluated using "position-grouped" cross-validation under the assumption that grouping variants by amino acid position prevents data leakage. 

However, our label-free spatial separation diagnostic (`logs/split_separation_diagnostic.json`) revealed that:
1. **The sequence-neighbourhood information channel remains wide open under position-grouped splits:** In UL23, **96.5%** of held-out test positions lie within $\le 2$ residues of a training position (median sequence distance = 1.0 residue). In UL30, **63.1%** lie within $\le 2$ residues (median distance = 2.0 residues).
2. Consequently, models with local sequence sensitivity can interpolate across position-grouped splits without learning generalizable biophysical rules.
3. Conversely, comparing unfitted zero-shot scores (like raw BLOSUM62) to supervised tabular models across splits is conceptually flawed: zero-shot scores are invariant to split schemes by construction ($\Delta_{\text{degrade}} \equiv 0.0000$).

To provide a rigorous, apples-to-apples evaluation, this protocol pre-registers a formal **Feature-Family $\times$ Split-Scheme Degradation Benchmark** in which all feature families are trained under identical supervised cross-validation across three split schemes on the exact same folds.

---

## 2. Experimental Design & Split Schemes

Evaluated using 5-fold cross-validation pre-computed in `data/processed/primary_3d_spatial_folds.tsv`:

1. **Split Scheme 1: Position-Grouped 5-Fold (`cv_fold_id`)**
   Variants at the exact same residue position are assigned to the same fold, but neighbouring positions remain in train folds (median train distance = 1–2 residues).
2. **Split Scheme 2: 3D k-Means Spatial 5-Fold (`spatial_cluster_id`)**
   K-means clustering on C$\alpha$ 3D coordinates.
3. **Split Scheme 3: Contiguous Sequence Block 5-Fold (`seq_block_id`)**
   Contiguous amino acid blocks held out (median train distance = 22–81 residues), strictly closing the local sequence-neighbourhood channel.

---

## 3. Supervised Feature Families & Negative Control

All models are trained strictly on training folds and evaluated on held-out test folds (out-of-fold pooled AUROC):

### Family A: Positional Identity Models
- Input: Normalized sequence position $p_i / L_{\text{protein}}$.
- Mechanism: Tests pure spatial/positional memorization and local smoothing.

### Family B: Evolutionary Conservation & Substitution Severity
- Inputs: `score_blosum62`, `score_neg_entropy`, `score_weighted_cons`, `score_pssm`.
- Mechanism: Classical alignment-derived conservation independent of target coordinates.

### Family C: Protein Language Model Representations
- Inputs: Zero-shot masked marginals and embeddings from `score_esm2_650m`, `score_protbert`, `score_esm2_35m`.
- Mechanism: Tests whether deep pre-trained contextual representations generalize across unseen sequence blocks.

### Family D: Multi-Modal Combined Representation
- Inputs: All features combined (Position + Conservation + PLMs).

### Negative Control: Fixed Zero-Shot Baseline
- Input: Unfitted zero-shot `score_blosum62` evaluated directly across folds without model fitting.
- Serves as an explicit negative control: proves mathematically that $\Delta_{\text{degradation}} \equiv 0.0000$ when models are not trained on training sets.

---

## 4. Supervised Model Architectures

To ensure conclusions do not depend on tree-based vs. linear inductive biases, two model families will be evaluated:
1. **Primary Model:** `RandomForestClassifier(n_estimators=100, max_depth=5, random_state=42)` (representing tree-based non-linear tabular models).
2. **Secondary Model:** `LogisticRegression(penalty='l2', C=1.0, random_state=42)` (representing regularized linear models).

---

## 5. Statistical Estimator & Paired Bootstrap Contrasts

1. **Metric:** Out-of-fold (OOF) AUROC pooled across all 5 test folds.
2. **Degradation Delta:**
   $$\Delta_{\text{degrade}} = \text{AUROC}(\text{position-grouped}) - \text{AUROC}(\text{contiguous-block})$$
   $$\Delta_{\text{spatial}} = \text{AUROC}(\text{position-grouped}) - \text{AUROC}(\text{3D k-means})$$
3. **Paired Bootstrap:** $B = 2,000$ bootstrap resamples of variants with replacement (seed `42`), computing paired differences in AUROC within each resample.

---

## 6. Confirmatory Hypotheses & Decision Rules

1. **Hypothesis 1 (Catastrophic Positional Collapse):**
   - **CONFIRMED:** Positional Identity models exhibit $\Delta_{\text{degrade}} \ge 0.30$ on UL30-ACV AND contiguous-block AUROC drops to near chance ($\le 0.550$).
   - **FAIL:** Contiguous-block AUROC remains $> 0.600$.

2. **Hypothesis 2 (PLM Representation Resilience):**
   - **CONFIRMED:** PLM representations (Family C) exhibit significantly attenuated degradation ($\Delta_{\text{degrade}} < 0.150$) AND retain strong discriminative power under contiguous holdout (contiguous-block AUROC $\ge 0.650$).
   - **FAIL:** PLMs collapse to near chance ($\le 0.550$).

3. **Hypothesis 3 (Target-Size Nuance on 3D Holdout):**
   - On UL30 (1,235 aa), 3D spatial holdout significantly degrades positional models ($\Delta_{\text{spatial}} \ge 0.200$).
   - On UL23 (376 aa), 3D spatial holdout retains substantial positional inflation due to sequence proximity leakage ($\Delta_{\text{spatial}} < \Delta_{\text{degrade}}$).

4. **Hypothesis 4 (Negative Control Invariance):**
   - Unfitted zero-shot BLOSUM exhibits $\Delta_{\text{degrade}} \equiv 0.0000$ ($p = 1.000$).

---

## 7. Execution & Traceability

- **Execution Script:** `scripts/run_split_degradation_benchmark.py`
- **Output Log:** `logs/split_degradation_benchmark_results.json`
- **Integration:** Integrated into `scripts/make_tables_from_logs.py` as Table 11.

