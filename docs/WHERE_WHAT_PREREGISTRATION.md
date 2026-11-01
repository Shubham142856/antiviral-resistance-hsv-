# Pre-Registration: The "Where × What" Tolerance-Window Hypothesis and Dual-Buffered Domain Evaluation Protocol

**Document Version:** 1.0 (Frozen Pre-Registration)  
**Date:** 2026-11-02  
**Target Study:** HSV-1 Antiviral Resistance Benchmark & Mechanism Evaluation  
**Status:** Pre-registered BEFORE Execution. The Git commit SHA and timestamp of this document serve as cryptographic proof of prior registration.  
**Repository:** [`https://github.com/Shubham142856/antiviral-resistance-hsv-`](https://github.com/Shubham142856/antiviral-resistance-hsv-)

---

## 1. Scientific Background & Rationale

Machine learning models for viral drug resistance are frequently evaluated using position-grouped cross-validation, which fails to neutralize active-site spatial interpolation shortcuts. When contiguous sequence blocks or structural domains are held out, tabular classifiers collapse to chance (AUROC = 0.5000), confirming that the holdout eliminates position memorization.

However, zero-shot evolutionary scores and chemical substitution matrices exhibit a striking divergence between drug targets:
- In **thymidine kinase (UL23)**, a non-essential viral prodrug activator, amino acid radicalness (BLOSUM62) predicts resistance with AUROC $\approx 0.6380$.
- In **DNA polymerase (UL30)**, an essential replicative enzyme, radicalness operates at chance (AUROC $\approx 0.4851$–$0.5081$).

Standard AUROC assumes a strictly monotonic relationship between score and resistance probability. We hypothesize that this apparent divergence is not an artifact or a simple absence of signal, but arises from fundamentally distinct biophysical constraint geometries:
1. **UL23 (Dispensable Activator):** Resistance is driven by enzymatic loss-of-function. Any mutation that destabilizes the protein or disrupts the active site confers resistance without affecting viral replication competence in cell culture. Resistance probability increases **monotonically** with mutation severity.
2. **UL30 (Essential Polymerase):** The viral polymerase cannot tolerate severe catalytic disruption; inactivating mutations are lethal to the virus. Conversely, mild substitutions fail to prevent drug binding. Therefore, at structurally and evolutionarily constrained positions, resistance is restricted to a **narrow tolerance window of intermediate severity**—subtle conformational remodeling that permits continued replication while discriminating against inhibitors.

This document formally pre-registers the evaluation protocol, feature extraction methods, statistical models, decision thresholds, and explicit falsification criteria before execution.

---

## 2. Structural Domain Boundaries

Domain definitions are anchored to published structural determinations of the HSV-1 DNA polymerase holoenzyme and thymidine kinase, corroborated by CATH and InterPro coordinate annotations:
- **UL30 Reference:** Shenping Liu et al. (*J. Biol. Chem.* 2006; 281:18193–18200; PDB `2GV9`), Sundaresh Shankar et al. (*Cell* 2024; 187:5572–5586.e15; PDB `8V1Q`/`8EXX`), and CATH/InterPro accessions (`G3DSA:3.30.342.10`, `G3DSA:3.30.420.10`, `G3DSA:3.90.1600.10`, `G3DSA:1.10.287.690`, `G3DSA:1.10.132.60`).
- **UL23 Reference:** Bennett et al. (*FEBS Lett.* 1999; 443:121–125; PDB `1OF1`) and InterPro (`IPR002484`).

### Table 1: Structural Domain Definitions

| Gene | Domain / Subdomain | Primary Residue Range | Length (aa) | Structural & Functional Role |
| :--- | :--- | :---: | :---: | :--- |
| **UL30** | **Pre-NH2 & N-Terminal (NTD)** | 1–362 | 362 | Extended herpesvirus N-terminal architecture (1–140 pre-NH2; 141–362 core NTD). |
| | **3′–5′ Exonuclease (Exo)** | 363–593 | 231 | Proofreading exonuclease active site; excised mispaired nucleotides. |
| | **Palm Domain (Segment 1)** | 594–766 | 173 | Catalytic core segment 1; coordinates active-site aspartates. |
| | **Fingers Domain** | 767–825 | 59 | Mobile subdomain; samples open/closed states; nucleotide/analog coordination. |
| | **Palm Domain (Segment 2)** | 826–956 | 131 | Catalytic core segment 2; polymerase active-site motifs. |
| | **Thumb & C-Terminal** | 957–1235 | 279 | Duplex DNA binding (957–1197) and processivity factor UL42 binding (1200–1235). |
| **UL23** | **N-Terminal Region** | 1–45 | 45 | Flexible non-catalytic N-terminal tail. |
| | **ATP-Binding Core** | 46–160 | 115 | Walker A motif / phosphate-binding loop (P-loop: residues 56–63). |
| | **Nucleoside Pocket** | 161–240 | 80 | Thymidine/aciclovir nucleoside coordination pocket (Arg163, Glu225). |
| | **Lid / Dimerization** | 241–320 | 80 | Dimer interface and flexible lid domain enclosing active site. |
| | **C-Terminal Region** | 321–376 | 56 | Structural scaffolding C-terminus. |

*Composite 5-Domain Grouping for UL30:* To align with the classic Family B right-hand architecture, Palm Segment 1 (594–766) and Palm Segment 2 (826–956) are treated as a single composite **Palm Domain** ($N_{\text{length}} = 304$ aa).

---

## 3. Dataset Distribution & Evaluation Rules

The primary missense cohort ($N=929$; `data/processed/primary_position_folds.tsv`, SHA-256: `7c4d0a3de2cd171b07230757f44e0dfbd61c6a0a8570dc8f263195b354c66909`) stratifies across these structural domains as follows:

### Table 2: Domain Class Counts

| Target Pair | Domain | Sensitive ($y=0$) | Resistant ($y=1$) | Total ($N$) | Percent Resistant | Standalone AUROC Status |
| :--- | :--- | :---: | :---: | :---: | :---: | :--- |
| **UL30-Aciclovir** | NTD (1–362) | 84 | 0 | 84 | 0.0% | Undefined ($N_{\text{pos}}=0$) |
| ($N=332$, Pos=62) | Exo (363–593) | 29 | 8 | 37 | 21.6% | Evaluable |
| | Palm (594–766, 826–956) | 79 | 32 | 111 | 28.8% | Evaluable |
| | Fingers (767–825) | 2 | 16 | 18 | 88.9% | Evaluable |
| | Thumb & C-term (957–1235) | 76 | 6 | 82 | 7.3% | Evaluable |
| **UL30-Foscarnet** | NTD (1–362) | 81 | 0 | 81 | 0.0% | Undefined ($N_{\text{pos}}=0$) |
| ($N=311$, Pos=33) | Exo (363–593) | 31 | 2 | 33 | 6.1% | Underpowered ($N_{\text{pos}} < 5$) |
| | Palm (594–766, 826–956) | 80 | 20 | 100 | 20.0% | Evaluable |
| | Fingers (767–825) | 6 | 8 | 14 | 57.1% | Evaluable |
| | Thumb & C-term (957–1235) | 80 | 3 | 83 | 3.6% | Underpowered ($N_{\text{pos}} < 5$) |
| **UL23-Aciclovir** | N-term (1–45) | 32 | 0 | 32 | 0.0% | Undefined ($N_{\text{pos}}=0$) |
| ($N=286$, Pos=120) | ATP Core (46–160) | 38 | 52 | 90 | 57.8% | Evaluable |
| | Nucleoside Pocket (161–240) | 15 | 45 | 60 | 75.0% | Evaluable |
| | Lid / Dimer (241–320) | 45 | 17 | 62 | 27.4% | Evaluable |
| | C-term (321–376) | 36 | 6 | 42 | 14.3% | Evaluable |

### Protocol Rules for Sparse & Single-Class Folds:
1. **Primary Metric:** The primary evaluation metric is **Pooled Out-of-Fold AUROC**, computed by concatenating all out-of-fold predictions across test domains and calculating AUROC on the full cohort ($N=332$ for UL30-ACV, $N=311$ for UL30-FOS, $N=286$ for UL23-ACV).
2. **Confidence Intervals:** 95% confidence intervals are computed via **Domain-Clustered Bootstrap** (1,000 resamples of structural domains with replacement).
3. **Fold-Mean Caveat:** Any held-out domain with $N_{\text{pos}} < 5$ (e.g. NTD with 0 positives; FOS Exo with 2; FOS Thumb with 3) cannot sustain a reliable fold-level AUROC. Fold-averaged AUROCs are reported strictly on evaluable domains ($N_{\text{pos}} \ge 5$), and the number of evaluable domains is explicitly documented.

---

## 4. Dual-Buffer Protocol & Split-Validity Rule

To guarantee that leave-one-domain-out holdout strictly eliminates both 1D sequence-adjacency and 3D spatial-proximity leakage, we implement a **dual exclusion buffer**. When evaluating held-out test domain $D$, training residues in surrounding regions are excluded from the training set.

### Fixed Buffer Search Grid:
- Sequence buffer $\Delta_{\text{seq}} \in \{5, 10, 15\}$ residues flanking domain boundaries.
- 3D spatial buffer $\Delta_{\text{3D}} \in \{8, 10, 12\}\text{ \AA}$ (minimum Euclidean distance between any training $C_\alpha$ atom and any test domain $C_\alpha$ atom in `8V1Q` or `1OF1`).

### Operational Validity Check:
- A "Position-Only Random Forest" (trained exclusively on the residue index `pos`) is evaluated on the buffered split.
- **Validity Threshold:** Position-Only RF AUROC must fall within $[0.45, 0.55]$ (consistent with pure chance).
- **Deterministic Selection Rule:** We pre-register that we choose the **minimal buffer** (smallest sequence exclusion + smallest spatial exclusion) that satisfies the validity check on the primary target (UL30-Aciclovir). That exact selected buffer is then applied fixed without re-tuning to UL30-Foscarnet and UL23-Aciclovir.

---

## 5. "Where × What" Feature Definitions

Every variant is decomposed along two orthogonal biophysical axes:

### A. "Where" (Positional Evolutionary & Structural Constraint)
To avoid sample-size noise and selection bias from averaging only observed clinical variants, positional constraint is calculated across **all 19 possible non-wildtype substitutions** scored by ESM-2 (650M):
$$\mu_{\text{PLM}}(p) = \frac{1}{19} \sum_{a \ne \text{wt}(p)} \Delta\text{LLR}(p, a)$$
where $\Delta\text{LLR}(p, a) = \log P(\text{wt} \mid \text{seq}) - \log P(a \mid \text{seq})$. Higher values represent residues where substitutions are broadly detrimental (high evolutionary constraint).

Additionally, structural distance $d_{\text{site}}(p)$ is computed as the minimum Euclidean distance (in Å) from $C_\alpha(p)$ to:
- Foscarnet ligand in `8EXX` and the catalytic carboxylates (Asp717, Asp888, Glu889) for UL30.
- (South)-methanocarba-thymidine ligand SCT in `1OF1` for UL23.

A single composite positional constraint score is defined as:
$$\text{Where}(p) = Z(\mu_{\text{PLM}}(p)) - Z(d_{\text{site}}(p))$$
where $Z(\cdot)$ denotes standardization to mean 0, variance 1. High `Where` indicates high evolutionary conservation and close proximity to the active site.

### B. "What" (Mutation Severity within Position)
Within-residue substitution severity is calculated relative to the position mean:
$$\text{Severity}_{\text{PLM}}(p, \text{mut}) = \Delta\text{LLR}(p, \text{mut}) - \mu_{\text{PLM}}(p)$$
and compared against chemical dissimilarity:
$$\text{Severity}_{\text{BLOSUM}}(p, \text{mut}) = -\text{BLOSUM62}(\text{wt}, \text{mut})$$
Both severity scores are standardized ($Z$-scored) across the dataset.

---

## 6. Statistical Models & Explicit Falsification Criteria

### Model Formulation:
To prevent overfitting on targets with limited events, we fit an $L_2$-penalized logistic regression model (Ridge, $C=1.0$):
$$\text{logit}(P(\text{Resistant})) = \beta_0 + \beta_1 \cdot \text{Where} + \beta_2 \cdot \text{Severity} + \beta_3 \cdot \text{Severity}^2 + \beta_4 \cdot (\text{Where} \times \text{Severity}) + \beta_5 \cdot (\text{Where} \times \text{Severity}^2)$$

Standard errors and 95% confidence intervals are computed via position-clustered bootstrap (1,000 resamples).

### Statistical Power Hierarchy:
- **Primary Hypothesis Target:** **UL30-Aciclovir** ($N=332$, 62 resistant events). With 6 parameters, this provides ~10.3 events per parameter, meeting standard guidelines for penalized logistic modeling.
- **Secondary Target:** **UL30-Foscarnet** ($N=311$, 33 resistant events). With ~5.5 events per parameter, FOS is acknowledged upfront as underpowered; it is evaluated as secondary confirmatory evidence.
- **Comparative Baseline Target:** **UL23-Aciclovir** ($N=286$, 120 resistant events; ~20 events per parameter).

### Pre-Registered Predictions:
1. **UL23 (Prodrug Activator):**
   - $\beta_2 > 0$: Resistance increases with severity.
   - $\beta_3 \approx 0$ and $\beta_5 \approx 0$: No inverted-U curvature; relationship is strictly monotonic.
2. **UL30 (Replicative Polymerase):**
   - At high constraint sites ($\text{Where} > 0$), the interaction curvature $\beta_5 < 0$: intermediate severity maximizes resistance, while extreme severity is non-viable.

### Explicit Falsification Rules:
- **UL23 Falsifier:** If the 95% CI for $\beta_3$ or $\beta_5$ strictly excludes 0 in the negative direction, the monotonic loss-of-function claim for UL23 is **falsified**.
- **UL30 Tolerance Window Falsifier:** If the 95% CI for the interaction curvature $\beta_5$ **includes 0 or is positive** on the primary target (UL30-Aciclovir), the tolerance window hypothesis is **falsified**. In that event, the paper will report the result as a refutation of the tolerance window hypothesis and retain the linear/monotonic baseline framing.

---

## 7. Cryptographic Freeze

Before executing the feature extraction and modeling scripts, this file is committed to Git and its SHA-256 hash recorded. All experimental results will be reported against these frozen criteria regardless of outcome.
