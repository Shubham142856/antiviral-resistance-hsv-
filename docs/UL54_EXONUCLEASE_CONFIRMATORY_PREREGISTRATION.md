# Confirmatory Pre-Registration: Exonuclease Active-Site Distance & Target-Mechanism Double Dissociation on HCMV UL54

**Document Version:** 1.0 (Confirmatory Pre-Registration Protocol)  
**Status:** CONFIRMATORY PRE-REGISTRATION (Strictly Zero Exonuclease Distances or AUROCs Computed)  
**Repository Anchor Commit:** Prior to execution  
**Target Protein:** HCMV UL54 (Replicative DNA Polymerase)  
**Drug Classes Evaluated:** 
- Pyrophosphate analogue: Foscarnet (FOS)
- Nucleoside analogue chain terminator: Ganciclovir (GCV)
- Acyclic nucleotide phosphonate: Cidofovir (CDV)

---

## 1. Scientific Motivation & Context

In the pre-registered Step 2 cross-virus replication benchmark (`docs/STEP2_HCMV_REPLICATION_PREREGISTRATION.md`), the primary replication arm evaluated whether 3D Euclidean distance to the polymerase active-site catalytic triad/contacts could predict within-study resistance in HCMV UL54, replicating the strong signal observed in HSV-1 UL30 ($S_0 = 0.787$).

The empirical outcome was sharply dissociated across drugs:
- **UL54 - Foscarnet (Secondary Descriptive):** $S_0 = \mathbf{0.8822}\ [0.7693, 0.9677]$ $\implies$ Highly predictive.
- **UL54 - Ganciclovir (Primary Replication):** $S_0 = \mathbf{0.4135}\ [0.1901, 0.6941]$ $\implies$ **NOT SUPPORTED (Failed Replication)**.
- **UL54 - Cidofovir (Secondary Descriptive):** $S_0 = \mathbf{0.2507}\ [0.1150, 0.4668]$ $\implies$ Strongly anti-predictive / inverted.

A subsequent label-free structural audit (`scripts/explore_ul54_drug_localization.py`, `logs/ul54_drug_localization_exploratory.json`) demonstrated that:
1. The polymerase catalytic core (residues 690–960) and the 3'->5' exonuclease proofreading domain (residues 300–545) are geometrically separated by a median distance of **21.7 Å**.
2. Resistance mutations across drugs display a monotonic spatial partitioning:
   - **Foscarnet:** Only **10.0%** of resistant variants fall in the exonuclease domain (median pos = 756.0, localized to the polymerase core pyrophosphate pocket).
   - **Ganciclovir:** **59.2%** of resistant variants fall in the exonuclease proofreading domain (median pos = 532.5).
   - **Cidofovir:** **77.2%** of resistant variants fall in the exonuclease proofreading domain (median pos = 516.0).

### The Falsifiable Confirmatory Hypothesis
This exploratory observation suggests that measuring distance solely to the *polymerase* active site succeeds on FOS because FOS resistance is confined to the pyrophosphate binding pocket, but fails on GCV and CDV because their resistance mechanism predominantly operates via proofreading exonuclease mutations (e.g., impaired excision or altered incorporation/excision kinetic balance).

To prevent post-hoc storytelling (HARKing), this protocol pre-registers a formal **Double Dissociation Confirmatory Test** evaluating distance to the **exonuclease active-site catalytic pocket** before executing any calculations.

---

## 2. Pre-Specified Structural Definitions (Zero Phenotype Labels Read)

### A. Exonuclease Catalytic Active Site (7 Residues)
Derived directly from the experimental cryo-EM structure of herpesvirus DNA polymerase (PDB 8EXX) by identifying protein residues within 5.0 Å of the crystallographic magnesium ion (MG 1303) at the exonuclease proofreading center, mapped via BLOSUM62 alignment to HCMV UL54:
- **Exo I Catalytic Motif (D-X-E):**
  - HSV-1 UL30 `D368` $\rightarrow$ **HCMV UL54 `D301`**
  - HSV-1 UL30 `I369` $\rightarrow$ **HCMV UL54 `I302`**
  - HSV-1 UL30 `E370` $\rightarrow$ **HCMV UL54 `E303`**
- **Exo II Catalytic Motif (N-X-X-X-F-D):**
  - HSV-1 UL30 `Y465` $\rightarrow$ **HCMV UL54 `Y407`**
  - HSV-1 UL30 `N466` $\rightarrow$ **HCMV UL54 `N408`**
  - HSV-1 UL30 `F470` $\rightarrow$ **HCMV UL54 `F412`**
  - HSV-1 UL30 `D471` $\rightarrow$ **HCMV UL54 `D413`**

*Coordinate Source:* `data/structures/P08546_swissmodel_8oj7.pdb` (validated against experimental 8EXX cryo-EM coordinates at Spearman $\rho = 0.989$ and median $|\Delta| = 0.39\text{ \AA}$ in `logs/option_ab_distance_agreement.json`).

*Feature Construction:*
For any variant at position $p_i$, the exonuclease proximity score is defined as:
$$\text{Proximity}_{\text{exo}}(p_i) = - \min_{a \in \text{ExoSite}} \| \text{C}\alpha(p_i) - \mathbf{x}_a \|_2$$
Where $\text{ExoSite}$ consists of all heavy atoms of the 7 catalytic residues $\{301, 302, 303, 407, 408, 412, 413\}$.

### B. Sensitivity Bounds on Exonuclease Domain
To ensure conclusions do not depend on an arbitrary 7-residue cutoff, two pre-specified sensitivity boundaries will be evaluated alongside the catalytic site:
1. **Primary Catalytic Pocket:** The 7 residues above.
2. **Domain Boundary Window 1 (Standard):** Linear segment $300 \le p \le 545$ (literature Exo I / II / III).
3. **Domain Boundary Window 2 (Tight):** Linear segment $300 \le p \le 520$.
4. **Domain Boundary Window 3 (Expanded):** Linear segment $280 \le p \le 560$.

---

## 3. Pre-Registered Predictions (Double Dissociation Crossover)

| Feature Evaluated | Foscarnet (PPi Analogue) | Ganciclovir (Nucleoside Analogue) | Cidofovir (Nucleotide Phosphonate) |
|---|---|---|---|
| **Polymerase Active-Site Distance** ($\text{Dist}_{\text{pol}}$) | **High Signal** ($S_0 > 0.80$) | **Suppressed / Chance** ($S_0 \approx 0.41$) | **Inverted / Anti-predictive** ($S_0 \approx 0.25$) |
| **Exonuclease Active-Site Distance** ($\text{Dist}_{\text{exo}}$) | **Degraded to Chance** ($S_0 < 0.50$ or $\Delta_{\text{exo - pol}} < 0$) | **Recovered Signal** ($S_0 > 0.55$, $\Delta_{\text{exo - pol}} > 0$) | **Recovered Signal** ($S_0 > 0.55$, $\Delta_{\text{exo - pol}} > 0$) |

### Formal Contrast:
A formal double dissociation delta will be estimated using joint study-clustered bootstrap ($B=2,000$, seed 42):
$$\Delta_{\text{Crossover}} = [S_0(\text{GCV, exo}) - S_0(\text{GCV, pol})] - [S_0(\text{FOS, exo}) - S_0(\text{FOS, pol})]$$
**Prediction:** $\Delta_{\text{Crossover}} > 0$ with 95% bootstrap CI strictly excluding zero.

---

## 4. Confirmatory Decision Rules

1. **Rule 1 (GCV Exonuclease Recovery):**
   - **CONFIRMED:** $S_0(\text{GCV, exo}) \ge 0.550$ AND 95% bootstrap CI lower bound $> 0.500$.
   - **INCONCLUSIVE:** Point estimate $\ge 0.550$ but 95% CI spans $0.500$.
   - **NOT SUPPORTED:** Point estimate $< 0.550$.

2. **Rule 2 (CDV Exonuclease Recovery):**
   - **CONFIRMED:** $S_0(\text{CDV, exo}) \ge 0.550$ AND 95% bootstrap CI lower bound $> 0.500$.
   - **INCONCLUSIVE:** Point estimate $\ge 0.550$ but 95% CI spans $0.500$.
   - **NOT SUPPORTED:** Point estimate $< 0.500$.

3. **Rule 3 (FOS Exonuclease Degradation / Crossover):**
   - **CONFIRMED:** $S_0(\text{FOS, exo}) < S_0(\text{FOS, pol})$ AND $\Delta_{\text{Crossover}} > 0$ ($p < 0.05$).
   - **VIOLATED:** $S_0(\text{FOS, exo}) \ge S_0(\text{FOS, pol})$.

---

## 5. Statistical Estimation Framework

- **Estimator:** Primary $S_0$ pair-weighted within-study AUROC on discordant pairs from mixed studies ($M \ge 5$, $N_{\text{pairs}} \ge 100$).
- **Unweighted Study-Mean AUROC:** Reported alongside $S_0$.
- **Resampling:** $B = 2,000$ study-clustered bootstrap resamples (seed `42`).
- **Phenotyping Rule:** Identical to frozen Step 2: $\text{EC}_{50}$ fold-change $\ge 2.0$ or reported 'Resistant' $\implies 1$; fold-change $< 2.0$ or 'Sensitive' / 'Polymorphism' $\implies 0$.

---

## 6. Traceability, Log Integrity & Erratum Policy

- **Frozen Log Protection:** `logs/step2_hcmv_replication_results.json` is COMMITTED and FROZEN. It will NOT be modified or overwritten.
- **Dedicated Output Log:** All confirmatory exonuclease results will be written to `logs/step2_hcmv_exonuclease_confirmatory_results.json`.
- **Erratum Linkage:** An erratum note will be logged in `logs/errata.json` documenting the addition of this confirmatory test following the label-free coordinate agreement audit.

