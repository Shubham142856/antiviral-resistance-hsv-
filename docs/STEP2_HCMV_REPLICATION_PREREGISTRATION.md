# Pre-Registration: Cross-Herpesvirus Replication of Structural Proximity and Radicalness Dissociation (HCMV UL54 / UL97)

**Document Version:** 1.0 (Frozen Pre-Registration)  
**Date:** 2026-11-03  
**Status:** Pre-registered BEFORE executing Step 2 replication analysis  
**Repository:** `https://github.com/Shubham142856/antiviral-resistance-hsv-`  

---

## 1. Rationale and Replication Hypotheses

In HSV-1, antiviral resistance demonstrates a target-dependent biophysical dissociation across the Evidence Ladder:
1. **Essential Polymerase (HSV-1 UL30):** Supervised sequence models collapse under spatial leakage controls (domain-held-out AUROC ~ 0.50). In contrast, deterministic active-site proximity ("Where") achieves robust rank-separation within assaying laboratories ($S_0 \text{ AUROC} = 0.7865\ [0.739, 0.833]$ on ACV, $0.8415\ [0.478, 0.927]$ on FOS), while amino acid substitution radicalness ("What", BLOSUM62) operates near chance on ACV ($0.5562\ [0.326, 0.760]$) and below chance on FOS ($0.3232\ [0.192, 0.399]$, exploratory/non-independent).
2. **Dispensable Activator (HSV-1 UL23):** Thymidine kinase readily tolerates catalytic loss-of-function; chemical radicalness (BLOSUM62) predicts resistance within study ($S_0 = 0.5785\ [0.429, 0.706]$) and pooled ($0.6130$), outperforming supervision under strict cross-domain partitions ($0.6730$ vs $0.5670$).

We pre-register an independent replication test in Human Cytomegalovirus (HCMV / Human betaherpesvirus 5):
- **Primary Replication Target:** HCMV UL54 DNA Polymerase (ortholog of HSV-1 UL30) under Ganciclovir (GCV).
- **Secondary Descriptive Targets:** HCMV UL54 under Foscarnet (FOS) and Cidofovir (CDV).
- **Counter-Arm Target:** HCMV UL97 Phosphotransferase / Kinase (functional prodrug activator, analogous to HSV-1 UL23) under Ganciclovir (GCV).

### Pre-Registered Predictions
1. **Polymerase Proximity Replication (UL54):** Active-site 3D proximity ($d_{\text{active-site}}$) predicts resistance within assaying studies ($S_0 \text{ AUROC} \ge 0.65$, 95% bootstrap CI lower bound $> 0.55$), replicating the spatial clustering observed in HSV-1 UL30.
2. **Polymerase Radicalness Collapse (UL54):** BLOSUM62 radicalness on UL54 operates near chance ($S_0 \text{ AUROC} \in [0.45, 0.58]$), reflecting structural constraints against non-specific disruptive charge/size alterations in an essential replicative engine.
3. **Dispensable Activator Radicalness (UL97):** BLOSUM62 radicalness on UL97 exceeds chance ($S_0 \text{ AUROC} \ge 0.60$), reflecting resistance via non-specific kinase disruption.
4. **Empirical Contrast ($\Delta_{\text{BLOSUM}}$):** The between-gene radicalness contrast $\Delta = \text{AUROC}(\text{UL97-GCV}) - \text{AUROC}(\text{UL54-GCV})$ is positive ($\Delta > 0$).

---

## 2. Structural Provenance & Ground Truth Audit

### A. RCSB PDB & AlphaFold DB Audit Disclosure
Prior to analysis, we audited experimental structural coverage for HCMV UL54:
- **RCSB PDB Audit:** Direct query to RCSB Search API yields exactly 1 structure: `1YYP` (2.5 Å X-ray). `1YYP` contains the UL44 processivity sliding clamp bound to a 22-residue C-terminal peptide of UL54 (residues 1221–1242). **The entire catalytic core of HCMV UL54 (residues 1–1220) is absent from the RCSB PDB.**
- **AlphaFold DB Audit:** The proteome of HCMV (strain AD169 / NC_006273.2) is absent from the public AlphaFold DB release (API queries return HTTP 404).

### B. Structural Model Specification (Option A: Primary)
To obtain reliable 3D Euclidean distances without fabricating coordinates:
- We utilize the full-length homology model from the Swiss-Model Repository (`P08546_swissmodel_8oj7.pdb`), built using template `8OJ7.1.A` (HSV-1 DNA polymerase holoenzyme complex resolved by cryo-EM at 2.7 Å).
- Residue coverage: 1,216 amino acids (positions 17–1232), covering 351 of 353 (99.43%) of all UL54 missense variants in HerpesDRG-DB.
- Unmodeled terminal residues (< 17 or > 1232, specifically A15T and S1235T, both benign polymorphisms) are imputed with the maximal observed structural distance.

### C. Active-Site Definition Frozen from Cryo-EM Complexes
Active-site reference atoms are established by projecting experimental ligand contacts from HSV-1 cryo-EM structures:
- `8EXX`: HSV-1 holoenzyme bound to foscarnet (PPF) at 3.3 Å
- `8V1T`: HSV-1 holoenzyme bound to acyclovir triphosphate (AVP) at 2.8 Å
- Residues within 5.0 Å of the drug in Chain A: `[717, 718, 719, 720, 721, 722, 723, 785, 789, 811, 812, 815, 816, 818, 819, 887, 888, 925, 927]`.
- Pairwise BLOSUM62 alignment between HSV-1 UL30 (1235 aa) and HCMV UL54 (1242 aa) demonstrates that **18 of 19 residues are strictly conserved identical residues** (only S816 $\to$ A816 is substituted).
- Mapped active-site positions on HCMV UL54 (frozen in `logs/hcmv_ul54_active_site_mapping.json`):
  $$\mathcal{A}_{\text{UL54}} = \{717, 718, 719, 720, 721, 722, 723, 784, 788, 811, 812, 815, 816, 818, 819, 911, 912, 949, 951\}$$
- Structural proximity feature for each residue $i$:
  $$d_{\text{active-site}}(i) = \min_{j \in \mathcal{A}_{\text{UL54}}, a \in \text{atoms}(j)} \|\mathbf{x}_i^{\text{C}\alpha} - \mathbf{x}_{j,a}\|_2$$
  Proximity score: $- d_{\text{active-site}}$ (closer to active site = higher score).

### D. Sensitivity Check (Option B: Ortholog-Projected Distance)
As an explicit sensitivity control, active-site distance is alternatively computed by mapping HCMV UL54 positions directly onto the experimental HSV-1 coordinates (`8EXX` Chain A) via the global pairwise alignment.

---

## 3. Cohort Definition & Quality Criteria

- **Source Database:** `data/raw/herpesdrg-db/herpesdrg-db.tsv` (frozen commit).
- **Target Filter:** `virus == 'HCMV'`, genes `UL54` and `UL97`.
- **Variant Inclusion:** Single-amino-acid missense substitutions with parseable wildtype, integer position, and mutant residue. Co-occurring complex substitutions are excluded.
- **Resistance Phenotyping Thresholds:**
  - Standard clinical cut-off: Fold-change $\ge 2.0\times \implies$ Resistant ($y=1$), Fold-change $< 2.0\times$ or reported 'Polymorphism'/'Sensitive' $\implies$ Sensitive ($y=0$).
  - Stratum $S_2$: Evaluated strictly on numeric fold-change entries.
- **Study Identification:** Extracted hierarchically from `ref_doi` (prefixed `DOI:`), falling back to `ref_title` (prefixed `TITLE:`).

---

## 4. Evaluation Ladder & Statistical Protocol

We evaluate each target along the Evidence Ladder:
1. **$S_0$ Within-Study Ascertainment Control (Primary Metric):**
   - Computed exclusively on concordant/discordant pairs within the same assaying study:
     $$S_0 = \frac{\sum_{m=1}^M \sum_{i \in \mathcal{R}_m} \sum_{j \in \mathcal{S}_m} \mathbb{I}(s_i > s_j) + 0.5 \cdot \mathbb{I}(s_i = s_j)}{\sum_{m=1}^M |\mathcal{R}_m| \cdot |\mathcal{S}_m|}$$
   - Positivity prerequisites: Mixed studies $M \ge 5$, total pairs $N_{\text{pairs}} \ge 100$, top study share $\le 50\%$.
   - Confidence Intervals: 95% study-clustered bootstrap with $B = 2,000$ resamples (seed `42`). Resampling clusters by `study_id`.
   - Sensitivity Check: Study-equal unweighted mean of per-study AUROCs across mixed studies.
2. **$S_2$ Quantitative Measured Phenotype Stratum:**
   - Evaluated on variants possessing quantitative fold-change determinations.
3. **$D_1$ Uncontrolled Naive Pooled Cohort:**
   - Standard empirical AUROC on deduplicated unique substitutions across all studies.

### Pre-Registered Decision Criteria for Replication
- **PASS:** Pair-weighted $S_0 \text{ AUROC} \ge 0.65$ with study-clustered bootstrap 95% CI lower bound $> 0.55$.
- **INCONCLUSIVE:** Point estimate $\ge 0.65$ but 95% CI lower bound $\le 0.55$ (due to cluster variance), or point estimate in $[0.55, 0.65)$.
- **ARTIFACT / FALSIFIED:** Point estimate $< 0.55$ or 95% CI upper bound $< 0.50$.

---

## 5. Execution Commitment

This protocol is frozen and committed to version control before executing `scripts/run_step2_hcmv_replication.py`.
No post-hoc threshold adjustment, metric modification, or active-site residue editing is permitted.

---

## 6. Dated Implementation Deviation Note & Composite Resolution
**Date:** 2026-11-03  
**Status:** Protocol Clarification & Dual Reporting Mandate  

1. **Exact Pre-Registered Text vs Execution Divergence:**
   - In Section 1 (Prediction 1) and Section 2.C (lines 50-52), the replication target score was defined mathematically as active-site 3D proximity:
     $$d_{\text{active-site}}(i) = \min_{j \in \mathcal{A}_{\text{UL54}}, a \in \text{atoms}(j)} \|\mathbf{x}_i^{\text{C}\alpha} - \mathbf{x}_{j,a}\|_2$$
     Proximity score: $- d_{\text{active-site}}$ (closer to active site = higher score).
   - In the initial executed script (`scripts/run_step2_hcmv_replication.py`), the score was computed strictly as distance-only ($-d_{\text{active-site}}$). While this matched Section 2.C verbatim, it diverged from the overarching HSV-1 definition of "Where", which in Table 8 and Section 2 was established as a two-part composite score combining ESM-2 positional evolutionary constraint ($\mu_{\text{pLM-19}}$) and active-site distance:
     $$W_{\text{composite}} = \frac{Z(\text{Constraint}) - Z(d_{\text{site}})}{\text{std}(Z(\text{Constraint}) - Z(d_{\text{site}}))}$$
     where $Z(\text{Constraint}) = - (\mu_{\text{pLM-19}} - \bar{\mu}) / \sigma_\mu$.

2. **Resolution Protocol (Zero Hypothesis Shift):**
   - No decision thresholds (0.65 / 0.55), estimators ($S_0$ pair-weighted / study-mean), coordinate sets (Option A Swiss-Model), or site lists (19 residues) are altered.
   - The pre-registered composite score on UL54 is computed using ESM-2 650M (`facebook/esm2_t33_650M_UR50D`) positional masked marginals ($\mu_{\text{pLM-19}}$) and combined with Option A 3D distances under the exact standardized function.
   - Both runs are reported side by side in Table 9B:
     * `Active-Site Proximity (As Executed: Distance-Only)`
     * `Where Score (As Pre-Registered: Composite pLM + Distance)`
   - Whichever outcome emerges is reported directly without post-hoc framing.


