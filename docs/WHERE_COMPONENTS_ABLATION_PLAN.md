# Where Decomposition Ablation Protocol: Constraint-Only vs. Distance-Only vs. Composite

**Document Version:** 1.2 (Exploratory Protocol Draft)  
**Status:** EXPLORATORY (Motivated by Post Hoc UL54 Observation; Zero Analysis Executed)  
**Primary Target:** HSV-1 UL30-Aciclovir ($S_0$ Pair-Weighted AUROC)  
**Descriptive Targets:** HSV-1 UL30-Foscarnet, HSV-1 UL23-Aciclovir  
**Rule Framework:** Purely descriptive/exploratory decomposition; NO confirmatory pass/fail thresholds applied.  

---

## 1. Scientific Motivation & Core Question

In the Step 2 cross-virus replication benchmark on HCMV UL54, evaluating components of the positional prior separately revealed an unexpected dissociation:
1. **Active-site proximity alone (Option A: Swiss-Model)** failed replication on the pre-registered primary arm ($S_0 = 0.4135$).
2. **Positional evolutionary constraint alone (ESM-2 650M masked marginals $\mu_{\text{pLM-19}}$)** retained strong discriminative power ($S_0 = 0.7619$).
3. **Structural Representation Invariance:** A label-free geometric audit (`logs/option_ab_distance_agreement.json`) confirmed that Swiss-Model coordinates and experimental cryo-EM coordinates (PDB 8EXX) agree at Spearman $\rho = 0.989$ and median $|\Delta| = 0.39\text{ \AA}$ across common residues, demonstrating that homology modeling did not corrupt distance geometry.
4. **Both coordinate options must be cited, not only Option A.** The frozen Step 2 log (`logs/step2_hcmv_composite_results.json`) contains two coordinate specifications on the same UL54 cohort, and on the Ganciclovir arm they disagree sharply: distance-only $S_0 = 0.4135$ under Option A against $0.7043$ under the original Option B. Quoting Option A alone would make the motivation selectively sourced. The divergence has since been resolved and is *not* a coordinate-provenance effect: the original Option B silently altered the reference atom set (ligand plus three catalytic residues in place of the 19-residue shell) and inadvertently included the distal exonuclease magnesium MG 1303, which sits 21.7 \AA\ from the polymerase shell (`logs/option_b_harmonization_and_apportionment_results.json`). With the shell harmonized, Option A and Option B agree to $\le 0.0092$ AUROC on every drug.
5. **The dissociation is drug-mechanism-dependent, not provenance-dependent.** With a single validated distance vector, $S_0$ is $0.8822$ on Foscarnet, $0.4135$ on Ganciclovir and $0.2507$ on Cidofovir, tracking the fraction of each drug's resistant variants that lie in the exonuclease region (10.0%, 59.2%, 77.2%; `logs/ul54_drug_localization_exploratory.json`). Any claim in this ablation about distance-term value must therefore be conditioned on inhibitor mechanism rather than on coordinate source.

This exploratory observation raises an essential architectural question for the foundational HSV-1 benchmark:
> **How much of the predictive signal in the zero-parameter positional prior (`Where`) on HSV-1 targets is driven by positional evolutionary constraint alone versus 3D active-site proximity alone, and does their combination yield a genuine super-additive or additive excess over constraint alone?**

---

## 2. Experimental Design & Target Prioritization

The ablation evaluates all three HSV-1 targets with strict hierarchical prioritization:

| Target System | Benchmark Role | Evaluation Hierarchy | Cohort Size ($N$) | Resistant / Sensitive | Mixed Studies ($M$) | Same-Study Pairs ($N_{\text{pairs}}$) |
|---|---|---|---|---|---|---|
| **HSV-1 UL30-Aciclovir** | Primary Essential Target | **Primary Arm** ($S_0$ Pair-Weighted) | 332 | 62 R / 270 S | 14 | 267 |
| **HSV-1 UL30-Foscarnet** | Secondary Essential Target | **Descriptive Only** | 311 | 33 R / 278 S | 8 | 82 |
| **HSV-1 UL23-Aciclovir** | Comparative Non-Essential Kinase | **Descriptive Only** | 286 | 120 R / 166 S | 17 | 331 |

*Methodological Note:* Because this analysis is an exploratory decomposition motivated by observations in HCMV UL54, **no confirmatory pass/fail thresholds** are imposed. All point estimates, paired deltas, and confidence intervals are reported descriptively.

---

## 3. Score Formulations Tested

For each variant $i$ at sequence position $p_i$:

### Score 1: Evolutionary Positional Constraint Only
$$Z(\text{Constraint}) = - \frac{\mu_{\text{pLM-19}}(p_i) - \bar{\mu}}{\sigma_\mu}$$
Where $\mu_{\text{pLM-19}}(p_i)$ is the mean zero-shot masked marginal log-likelihood ratio across all 19 alternative amino acids from ESM-2 650M at position $p_i$. Positions with greater intolerance to substitution receive higher constraint scores.

### Score 2: 3D Active-Site Proximity Only
Evaluated under both linear standardized distance and rank-normalized distance:
1. **Standardized Distance Formulation:**
   $$Z(d) = \frac{d(i) - \bar{d}}{\sigma_d}, \quad \text{Score}_{\text{dist}} = - Z(d)$$
2. **Rank-Normalized Distance Formulation:**
   $$\text{RankNorm}(d) = 1 - \frac{\text{rank}(d(i)) - 1}{N - 1}$$
   Where $d(i)$ is the minimum Euclidean distance in Angstroms from the C$\alpha$ of residue $p_i$ to any heavy atom of the holoenzyme active-site contact residues (PDB 8V1Q / 8EXX for UL30; 1OF1 for UL23).

### Score 3: Composite Where (Current Benchmark Specification)
$$W_{\text{composite}} = \frac{Z(\text{Constraint}) - Z(d)}{\text{std}(Z(\text{Constraint}) - Z(d))}$$

### Reference Baselines:
- **BLOSUM62 Radicalness:** Substitution severity independent of position.
- **Position-Only RF:** Evaluated under D1 pooled OOF and D2 dual-buffer holdout.

### Pre-Registered Coordinate-Quality Covariate

Fixed here before any ablation AUROC is computed, so that it cannot be selected after the fact. Should the distance term underperform on any arm, the following label-free statistic is the pre-specified explanatory covariate, and no alternative explanation of coordinate quality may be substituted later:

$$Q_{\text{shell}} = \text{median}_{j \in \mathcal{A}} \min_{k} \| \text{C}\alpha_j^{\text{model}} - \text{C}\alpha_k^{\text{experimental}} \|_2$$

the median C$\alpha$ deviation over the active-site shell residues $\mathcal{A}$ between the coordinate set used and the experimental reference. For HSV-1 both UL30 (`8V1Q`/`8EXX`) and UL23 (`1OF1`) coordinates are experimental, so $Q_{\text{shell}} \equiv 0$ by definition and the covariate cannot explain any HSV-1 result. This is stated explicitly to preclude a post hoc appeal to coordinate quality on HSV-1 arms: on this benchmark, a weak distance term must be explained by inhibitor mechanism or by genuine absence of spatial signal, not by geometry error.

**Distance rank-normalization.** UL23 distances span 0.0–30.9 \AA\ and UL30 distances span 0.0–68.3 \AA, so an unnormalized UL30-versus-UL23 distance comparison partly measures protein size rather than target dispensability. All cross-target distance contrasts are therefore computed on within-target rank-normalized distances, and the raw-scale values are reported alongside for transparency.

---

## 4. Statistical Estimators, Strata & Paired Clustered Bootstrap

1. **Primary Estimator: Within-Study $S_0$ (Pair-Weighted):**
   - Calculated strictly over same-study discordant pairs $(i \in \mathcal{R}_s, j \in \mathcal{S}_s)$:
     $$S_0 = \frac{\sum_{s=1}^M \sum_{i \in \mathcal{R}_s} \sum_{j \in \mathcal{S}_s} \mathbf{1}[S_i > S_j] + 0.5 \cdot \mathbf{1}[S_i = S_j]}{\sum_{s=1}^M |\mathcal{R}_s| \cdot |\mathcal{S}_s|}$$
   - Alongside unweighted study-mean AUROC: $\bar{A}_{\text{study}} = \frac{1}{M}\sum_{s=1}^M \text{AUROC}_s$.
2. **Paired Study-Clustered Bootstrap Engine:**
   - $B = 2,000$ study-clustered bootstrap resamples (fixed seed `42`).
   - In each resample iteration $b$, the identical resampled studies are evaluated simultaneously across all scores.
   - Paired deltas computed directly within each resample:
     $$\Delta_{\text{Composite - Constraint}}^{(b)} = S_0(W_{\text{composite}})^{(b)} - S_0(\text{Constraint})^{(b)}$$
     $$\Delta_{\text{Composite - Distance}}^{(b)} = S_0(W_{\text{composite}})^{(b)} - S_0(\text{Distance})^{(b)}$$
     $$\Delta_{\text{Distance - Constraint}}^{(b)} = S_0(\text{Distance})^{(b)} - S_0(\text{Constraint})^{(b)}$$
   - 95% bootstrap confidence intervals extracted as [2.5th, 97.5th] percentiles.
3. **Quality Covariates & Evaluated Strata:**
   - **S0 Full Cohort (Primary):** Same-study pairs across all mixed studies ($N=332$ on UL30-ACV).
   - **S2 Matched Quantitative Stratum:** High-confidence subset with numeric $\text{EC}_{50}$ fold-change measurements ($N=78$ on UL30-ACV, $N=52$ on UL30-FOS, $N=66$ on UL23-ACV).
   - **D1 Pooled Cohort:** Uncontrolled naive baseline for ascertaining study-clustering bias.

---

## 5. Descriptive Interpretation Matrix

| Outcome Pattern on UL30-ACV | Empirical Signature | Scientific Interpretation |
|---|---|---|
| **Additive / Complementary** | $S_0(W) > S_0(\text{Constraint})$ and $S_0(W) > S_0(\text{Distance})$; 95% CI of $\Delta > 0$ | Spatial 3D distance and evolutionary constraint provide non-redundant predictive signal. |
| **Constraint-Dominated** | $S_0(W) \approx S_0(\text{Constraint}) \gg S_0(\text{Distance})$; $\Delta_{\text{Comp - Const}} \approx 0$ | Predictive signal of `Where` is driven almost entirely by evolutionary conservation; active-site proximity adds negligible marginal value. |
| **Distance-Dominated** | $S_0(W) \approx S_0(\text{Distance}) \gg S_0(\text{Constraint})$; $\Delta_{\text{Comp - Dist}} \approx 0$ | Active-site pocket geometry drives resistance; evolutionary constraint is largely a proxy for proximity. |
| **Target Dissociation** | Distance contributes on UL30 (essential polymerase), but degrades on UL23 (dispensable kinase) | Structural proximity is functional only for targets where resistance requires maintaining catalytic competence. |

---

## 6. Execution Plan & Traceability

- **Script:** `scripts/run_where_components_ablation.py`.
- **Log Destination:** `logs/where_components_ablation_results.json`.
- **Generator Integration:** Summary table integrated into `scripts/make_tables_from_logs.py`.
- **Deterministic Verification:** Enforced by `scripts/test_table_generator_determinism.py`.
- **STOP AND VERIFY:** Protocol locked before running ablation calculations.
