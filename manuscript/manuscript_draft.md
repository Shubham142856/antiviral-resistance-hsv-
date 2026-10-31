# Spatial Shortcut Learning Inflates Viral Resistance Prediction: A Domain-Generalization Benchmark for Herpes Simplex Virus

**Authors:** [Author List to be confirmed per PI sign-off]  
**Correspondence:** [Corresponding Author]  
**Target Venue:** *Briefings in Bioinformatics* / *PLOS Computational Biology*  
**Repository & Data:** [`https://github.com/shubhamraikwar142/antiviral-resistance-hsv-1`](https://github.com/shubhamraikwar142/antiviral-resistance-hsv-1)  
**Dataset Integrity Hash:** SHA-256 `7C4D0A3DE2CD171B07230757F44E0DFBD61C6A0A8570DC8F263195B354C66909`  
**Audit Protocol:** `AUDIT_PROTOCOL.md` (Strict Version 4, Frozen)

---

## Abstract

Machine learning models for viral drug resistance are frequently evaluated using position-grouped cross-validation, where individual sequence positions are held out during training. Here, we demonstrate that position-grouped cross-validation fails to remove spatial interpolation shortcuts in Herpes Simplex Virus (HSV-1 and HSV-2) resistance benchmarks. A baseline random forest model trained exclusively on sequence position index achieves an apparent AUROC of **0.8715 $\pm$ 0.0171** under position-grouped cross-validation on HSV-1 DNA polymerase (UL30) treated with aciclovir, but collapses completely to **0.5001 $\pm$ 0.0009** (chance) under a contiguous-block structural holdout that isolates protein domains.

To establish an evaluation standard, we constructed a verified, deduplicated benchmark of 929 single-site missense substitutions across thymidine kinase (UL23; aciclovir, $N=286$) and DNA polymerase (UL30; aciclovir, $N=332$; foscarnet, $N=311$). Under this contiguous-block protocol, zero-shot protein language models (PLMs) maintain robust predictive signal where tabular machine learning models collapse to chance. ESM-2 (650M) achieves contiguous-block mean AUROCs of **0.7397 $\pm$ 0.1565** (UL23-aciclovir), **0.7875 $\pm$ 0.1772** (UL30-aciclovir), and **0.8809 $\pm$ 0.1250** (UL30-foscarnet), passing our pre-registered stopping rule (mean $\ge 0.65$, 95% $t$-CI lower bound $> 0.50$) on all three target pairs.

By benchmarking against multi-herpesvirus evolutionary baselines and BLOSUM62 amino acid radicalness, we uncover a statistically robust biological dissociation between the two resistance targets. **BLOSUM62 radicalness alone separates UL23 from UL30**: on thymidine kinase, radicalness achieves **0.6380** (95% CI [0.5524, 0.7237], every block above 0.55), while on DNA polymerase it collapses to chance (**0.5074 to 0.5494**, CIs spanning 0.50). This divergence reflects the distinct resistance mechanisms of the two enzymes: UL23 resistance arises primarily through disruptive substitutions that abolish prodrug activation, whereas UL30 resistance requires subtle conformational changes that preserve replication while discriminating against inhibitors. Under contiguous-block holdout, zero-shot PLMs also outperform alignment-based conservation baselines on UL30 by **$+0.0563$ to $+0.0876$**, but this residual advantage is **directional rather than statistically decisive** (UL30-ACV $p=0.174$; UL30-FOS $p=0.373$; both 95% CIs span zero at $n=4$ evaluable blocks).

In a prospective temporal holdout trained on pre-2016 variants ($N=622$) and tested on post-2016 variants ($N=307$), performance stratifies sharply: models achieve AUROCs of **0.8872 to 0.9286** at previously seen positions, but drop to **0.6996 to 0.8126** at completely novel positions, broadly consistent with our contiguous-block holdout estimates. Truncating and frameshift variants represent **44.4%** of cataloged clinical UL23 resistance; our benchmark addresses the remaining **55.6%** missense class where predictive modeling is biologically non-trivial.

---

## 1. Introduction

Antiviral drug resistance in Herpes Simplex Virus types 1 and 2 (HSV-1 and HSV-2) poses significant clinical challenges, particularly in immunocompromised individuals, hematopoietic stem cell transplant recipients, and patients with recurrent herpes stromal keratitis. First-line therapy relies on aciclovir (ACV), a guanosine nucleoside analogue that requires initial monophosphorylation by the viral thymidine kinase (TK, encoded by the *UL23* gene), followed by conversion by host cellular kinases to aciclovir triphosphate, which selectively inhibits and chain-terminates the viral DNA polymerase catalytic subunit (encoded by the *UL30* gene). Second-line therapy for ACV-resistant infections primarily employs foscarnet (FOS, phosphonoformic acid), a pyrophosphate analogue that directly blocks the pyrophosphate-binding pocket of UL30 without requiring phosphorylation by UL23.

While phenotypic susceptibility testing via plaque reduction assays (PRA) or dye uptake remains the diagnostic gold standard, phenotypic assays require 1–2 weeks, prompting extensive efforts toward sequence-based genotypic resistance prediction. In recent years, numerous machine learning classifiers and deep learning models have been proposed for variant effect and resistance prediction across viral pathogens. However, variant-level machine learning benchmarks frequently suffer from unrecognized data leakage.

### The Spatial Leakage Problem in Resistance Modeling
Cross-validation protocols that partition datasets randomly across rows fail catastrophically in mutation benchmarks because different amino acid substitutions frequently occur at the same sequence position. When the same position appears in both training and test sets, models easily memorize position identity rather than learning generalizable principles of mutation effect. To address this, "position-grouped" cross-validation—holding out individual positions entirely—has become the standard recommendation.

However, position-grouped holdout does not eliminate spatial shortcut learning. The PRIME framework previously demonstrated that position-stratified validation exposes a related spatial leakage vulnerability in SARS-CoV-2 variant prediction. In this work, we demonstrate this vulnerability in Herpes Simplex Virus antiviral resistance and show that contiguous-block holdout provides the necessary structural evaluation standard. In globular viral enzymes such as DNA polymerase, resistance-conferring mutations cluster within localized structural domains (such as the catalytic Palm and Finger domains). Tree-based models or overparameterized neural networks can interpolate resistance phenotypes from adjacent or spatially neighboring positions. When tested on held-out positions that reside within known active-site clusters, models exploit this spatial interpolation shortcut, achieving high apparent test AUROCs while failing to generalize to novel regions.

### Contributions of this Work
1. **Empirical Demonstration of Spatial Shortcut Learning:** We show that a "position-only" random forest, trained solely on sequence index, scores **0.8715** under position-grouped cross-validation, but collapses to **0.5001** under a domain-level contiguous-block structural holdout.
2. **Standardized Missense Resistance Benchmark:** We establish a curated, deduplicated primary cohort of 929 single-site missense substitutions across UL23 and UL30, verified with cryptographic provenance (SHA-256 hash `7C4D...`).
3. **Radicalness Dissociation between Resistance Targets (Statistically Significant):** BLOSUM62 amino acid radicalness separates the two targets: it achieves **0.6380** (95% CI [0.5524, 0.7237], every block $> 0.55$) on UL23, but collapses to chance (**0.5074 to 0.5494**) on UL30. This is the statistically robust component of the biological dissociation and reflects the distinct resistance mechanisms of prodrug-activating kinase versus replicative polymerase.
4. **Zero-Shot PLMs Generalize Across Structural Blocks:** Zero-shot protein language models (ESM-2 650M and ProtBERT) maintain high predictive performance under contiguous-block holdout (0.74 to 0.88 AUROC) where tabular models fail. On UL23, PLMs tie alignment baselines; on UL30, they retain a directional advantage over all alignment baselines tested, but this advantage is not individually statistically significant at $n=4$ evaluable blocks.
5. **Prospective Temporal Validation:** Evaluating on prospective post-2016 variants confirms that performance at novel sequence positions matches our contiguous-block estimates, establishing the real-world predictive validity of domain holdout.
6. **Quantification of Clinical Scope:** We document that truncating and frameshift variants represent 44.4% of cataloged UL23 resistance, establishing that our benchmark rigorously evaluates the missense class where predictive modeling is genuinely needed.

---

## 2. Materials and Methods

### 2.1 Dataset Curation and Deduplication
The benchmark dataset was constructed by integrating the curated HerpesDRG database with comprehensive phenotypic resistance compilations from clinical literature (including the CMI 2025 systemic review by Döhne et al., Sauerbrei et al. JAC 2016, and Burrel et al. Antiviral Res 2010). All entries were merged on unique `(Virus, Gene, Drug, Mutation)` tuples.

Phenotypic resistance labels were standardized against established clinical consensus criteria:
- **Aciclovir (PRA / Dye Uptake):** The primary operational criterion is fold-change $\ge 3.0\times$ in $\text{IC}_{50}$ relative to the wild-type reference strain, or qualitative curation as Resistant versus Sensitive/Polymorphism in HerpesDRG and CMI 2025. In the clinical literature, the foundational absolute PRA cut-off is $\text{IC}_{50} \ge 2.0\text{ }\mu\text{g/mL}$ (established by Safrin et al., *Antimicrob Agents Chemother* 1994;38:1246–1250 [PMID: 8092821], where it demonstrated a 95% predictive value for lesion failure to heal, and compiled in reviews such as Sauerbrei et al., *J Antimicrob Chemother* 2016;71:6–16).
- **Foscarnet (PRA / Dye Uptake):** The primary operational criterion is fold-change $\ge 3.0\times$ in $\text{IC}_{50}$ relative to the wild-type reference strain, or qualitative curation as Resistant versus Sensitive/Polymorphism in HerpesDRG and CMI 2025. The absolute cut-off most consistently documented in clinical trials is $\text{IC}_{50} \ge 100\text{ }\mu\text{g/mL}$, first established in the landmark 115-patient HIV cohort of Safrin et al. (1994), where $\text{IC}_{50} \ge 100\text{ }\mu\text{g/mL}$ predicted an 88% failure-to-heal rate versus 0% below it (recently reviewed in Wang, *J Clin Microbiol* 2024;62:e00728-23). Published foscarnet cut-offs are stated in mass concentration ($\mu\text{g/mL}$) without specifying chemical hydration state; we strictly avoid converting them to molar units. Because absolute concentration thresholds depend on laboratory assay conditions, our benchmark evaluates the normalized fold-change metric ($\ge 3.0\times$), the standard reproducible criterion across multicenter studies.

Conflicted entries (where identical variants had discordant reported phenotypes across laboratories) were logged and strictly excluded from the primary cohort. Co-occurring combinations and multi-mutation entries were segregated into a dedicated epistasis set (`co_occurrence_entries.tsv`, $N=125$).

### 2.2 Primary Evaluation Cohort and Scope Limitation
The primary cohort was restricted strictly to single-site missense substitutions in HSV-1 reference coordinates (UL23 reference NCBI `NC_001806.2`, 376 aa; UL30 reference NCBI `NC_001806.2`, 1,235 aa). The primary cohort comprises 929 rows:
- **UL23-Aciclovir:** 286 rows (120 Resistant, 166 Sensitive; resistance rate 42.0%; 210 unique positions).
- **UL30-Aciclovir:** 332 rows (62 Resistant, 270 Sensitive; resistance rate 18.7%; 297 unique positions).
- **UL30-Foscarnet:** 311 rows (33 Resistant, 278 Sensitive; resistance rate 10.6%; 280 unique positions).

The complete dataset is frozen and tracked by its SHA-256 digest: `7C4D0A3DE2CD171B07230757F44E0DFBD61C6A0A8570DC8F263195B354C66909`.

#### Scope Quantification (The Truncating Class)
In our curated unified database of 941 UL23 entries, 252 entries are non-missense variants (168 nonsense stop codons, 68 frameshifts, 16 indels). Among all cataloged resistant UL23 variants, truncating/frameshift variants account for **186 / 419 (44.4%)**. These variants were excluded by construction from the benchmark cohort because protein language model masked-marginal likelihood ratios ($\Delta\text{LLR}$) are undefined for non-single-substitution variants. Clinically, premature termination abolishes thymidine kinase catalytic function; resistance in this class is deterministically known without machine learning. Our benchmark evaluates missense variant prediction, where clinical interpretation is actively challenging.

### 2.3 Evaluation Schemes and Stopping Rule

#### 1. 5-Fold Position-Grouped Cross-Validation (`cv_fold_id`)
Positions were partitioned into 5 disjoint folds using deterministic modulo hashing on position index: $\text{fold} = (\text{pos} \bmod 5)$. No residue position appears in both training and test sets.

#### 2. 5-Fold Contiguous-Block Cross-Validation (`block_fold_id`)
To evaluate structural domain holdout, sequence positions were partitioned into 5 contiguous structural blocks:
- **UL23 (376 aa):** Block $b = \min(4, \lfloor(\text{pos} - 1) / 75\rfloor)$ (75 aa per block).
- **UL30 (1,235 aa):** Block $b = \min(4, \lfloor(\text{pos} - 1) / 247\rfloor)$ (247 aa per block).

In UL30, catalytic resistance mutations concentrate heavily within the Palm and Finger subdomains (Blocks 2 and 3, aa 495–988). Block 0 (aa 1–247) contains 0 resistant mutations (all 63 rows are sensitive polymorphisms), rendering AUROC undefined in that block. Block 0 was treated as single-class and omitted from block AUROC averages, yielding 4 evaluable blocks for UL30.

#### Primary Evaluation Metric & Confidence Intervals
The primary metric is the **contiguous-block mean AUROC** across evaluable blocks. Confidence intervals are computed using a **block-level Student's $t$-distribution with $n - 1$ degrees of freedom**:
$$\text{CI}_{95\%} = \bar{\text{AUROC}} \pm t_{0.975, n-1} \cdot \frac{s}{\sqrt{n}}$$
where $n=5$ for UL23 and $n=4$ for UL30. Pooled AUROCs across the entire dataset and position-cluster bootstrap CIs are reported strictly as secondary descriptive statistics.

#### Pre-Registered Stopping Rule
A model passes the benchmark if it achieves a contiguous-block mean AUROC $\ge 0.65$ with the lower bound of its 95% Student's $t$ confidence interval strictly exceeding 0.50, on at least 2 of the 3 primary pairs.

### 2.4 Multi-Herpesvirus Multiple Sequence Alignment
To test whether sequence models capture structural signal beyond phylogenetic conservation, we constructed reference-anchored multiple sequence alignments from complete viral genomes:
1. **HSV-1:** Human alphaherpesvirus 1 strain 17 (`NC_001806.2`)
2. **HSV-2:** Human alphaherpesvirus 2 strain HG52 (`NC_001798.2`)
3. **VZV:** Human alphaherpesvirus 3 strain Dumas (`NC_001348.1`; ORF36 for UL23, ORF28 for UL30)
4. **HCMV:** Human betaherpesvirus 5 strain Merlin (`NC_006273.2`; UL54 for UL30)
5. **HHV-6B:** Human betaherpesvirus 6B (`NC_000898.1`; U38 for UL30)

Alignments were computed using global dynamic programming with a BLOSUM62 matrix (gap open $-10.0$, gap extend $-0.5$). For each position $i$ in the HSV-1 reference, we computed four evolutionary and profile baselines:
1. **Family Conservation Fraction:** $C(i) = \frac{\sum_{k} \mathbb{I}(\text{res}_{k, i} == \text{wt})}{\text{non-gap orthologs}}$
2. **Shannon Entropy:** $H(i) = -\sum_{a \in \mathcal{A}} p(a) \log_2 p(a)$
3. **Multi-Herpesvirus Weighted Conservation:** $\text{Score}(i, \text{wt} \to \text{mut}) = \left(1.0 - \frac{H(i)}{\log_2(N)}\right) \times 4.0 - \text{BLOSUM62}(\text{wt}, \text{mut})$
   where $N$ is the total number of aligned ortholog genomes ($N=3$ for UL23: HSV-1, HSV-2, VZV; $N=5$ for UL30: HSV-1, HSV-2, VZV, HCMV, HHV-6B).
4. **Position-Specific Scoring Matrix (PSSM):** Dirichlet pseudocount-smoothed profile ($\alpha=1.0$) with Robinson & Robinson (1991) background amino acid frequencies $q(a)$:
   $$p_i(a) = \frac{c_i(a) + \alpha \cdot q(a)}{N_i + \alpha}, \quad \text{Score}(i, \text{wt} \to \text{mut}) = \log_2 \frac{p_i(\text{wt})}{p_i(\text{mut})}$$
Each of these alignment-derived scores is evaluated as an independent baseline under the same contiguous-block protocol applied to the language models, and we compare zero-shot PLMs against the strongest baseline on each target pair rather than against a single nominated comparator.

### 2.5 Zero-Shot Protein Language Models
We evaluated three masked language models without fine-tuning:
1. **ESM-2 650M** (`facebook/esm2_t33_650M_UR50D`, 33 layers, 650M parameters)
2. **ProtBERT** (`Rostlab/prot_bert`, 30 layers, 420M parameters)
3. **ESM-2 35M** (`facebook/esm2_t12_35M_UR50D`, 12 layers, 35M parameters)

Predictions were computed using masked-marginal log-likelihood ratios ($\Delta\text{LLR}$):
$$\Delta\text{LLR} = \log P(\text{mut} \mid \mathbf{x}_{\backslash i}) - \log P(\text{wt} \mid \mathbf{x}_{\backslash i})$$
The pre-registered direction specifies that deleterious mutations increase resistance likelihood; thus, the decision score is $-\Delta\text{LLR}$.

### 2.6 Supervised Predictor and Temporal Validation
Supervised models were trained using eight engineered features: zero-shot ESM-2 650M $-\Delta\text{LLR}$, ProtBERT $-\Delta\text{LLR}$, multi-herpesvirus weighted conservation, family conservation fraction, BLOSUM62 radicalness, $\Delta\text{Hydropathy}$ (Kyte-Doolittle), $\Delta\text{Molecular Weight}$, and $\Delta\text{Charge}$ (at pH 7.4). Models were trained using L2-regularized Logistic Regression ($C=0.1$) and Gradient Boosting (`HistGradientBoostingClassifier`).

Prospective temporal validation was conducted by training on variants published prior to 2016 ($N=622$) and evaluating on variants published from 2016 onward ($N=307$). Test performance was stratified into novel substitutions at previously seen positions versus variants at completely novel sequence positions.

---

## 3. Results

### 3.1 Position-Grouped Cross-Validation Masks Spatial Leakage
To measure whether position-grouped cross-validation prevents spatial shortcut learning, we trained a "Position-Only Random Forest"—a decision tree ensemble provided exclusively with the sequence index (`pos`) and zero biological features—and compared it against an All-Features Random Forest and zero-shot PLMs (Table 1).

**Table 1: The Dual-Split Leakage Finding.** Cross-validation performance under Position-Grouped CV versus Contiguous-Block CV across 929 single-site missense mutations.

| Target Pair | Model Architecture | Position-Grouped Mean $\pm$ SD | Contiguous-Block Mean $\pm$ SD | 95% Student's $t$ CI (Block) | Status under Block Holdout |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **UL30-Aciclovir** | Position-Only RF | **0.8715 $\pm$ 0.0171** | **0.5001 $\pm$ 0.0009** | — | **Collapses to chance** |
| | All-Features RF | 0.7328 $\pm$ 0.0365 | 0.4089 $\pm$ 0.0198 | — | Sub-chance |
| | **ESM-2 650M (Zero-shot)** | **0.8273 $\pm$ 0.0632** | **0.7875 $\pm$ 0.1772** | **[0.5055, 1.0000]** | **Robust generalization** |
| **UL30-Foscarnet** | Position-Only RF | **0.8684 $\pm$ 0.0282** | **0.5166 $\pm$ 0.0105** | — | **Collapses to chance** |
| | All-Features RF | 0.6756 $\pm$ 0.0446 | 0.6372 $\pm$ 0.0181 | — | Degraded |
| | **ESM-2 650M (Zero-shot)** | **0.8454 $\pm$ 0.1050** | **0.8809 $\pm$ 0.1250** | **[0.6821, 1.0000]** | **Robust generalization** |

Under position-grouped cross-validation, the position-only model achieves AUROCs of **0.8715** (UL30-ACV) and **0.8684** (UL30-FOS). Because catalytic resistance mutations reside in spatial proximity, the tree interpolates labels from neighboring active-site residues. Under contiguous-block holdout, this shortcut is neutralized, and the model collapses to **0.5001** and **0.5166**. This proves that position-grouped cross-validation is vulnerable to spatial confounding.

---

### 3.2 Zero-Shot Protein Language Models Generalize Beyond Spatial Shortcuts
In contrast to tabular models, zero-shot protein language models evaluate sequence context without memorizing regional position labels. Table 2 summarizes benchmark performance across all evaluated models under contiguous-block holdout.

**Table 2: Zero-Shot PLM Benchmark under Contiguous-Block Holdout.**

| Model | Target Pair | Contiguous-Block Mean $\pm$ SD | 95% Student's $t$ CI | Pooled AUROC [95% Boot CI] | Stopping Rule Status |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **ESM-2 650M** | **UL23-Aciclovir** | **0.7397 $\pm$ 0.1565** | **[0.5455, 0.9340]** | 0.7881 [0.7134, 0.8422] | **PASS** |
| *(Primary)* | **UL30-Aciclovir** | **0.7875 $\pm$ 0.1772** | **[0.5055, 1.0000]** | 0.8240 [0.7570, 0.8857] | **PASS** |
| | **UL30-Foscarnet** | **0.8809 $\pm$ 0.1250** | **[0.6821, 1.0000]** | 0.8409 [0.7433, 0.9162] | **PASS (3/3 Pairs)** |
| **ProtBERT 420M** | UL23-Aciclovir | 0.7530 $\pm$ 0.0713 | [0.6645, 0.8415] | 0.7243 [0.6533, 0.7858] | PASS |
| *(Secondary)* | UL30-Aciclovir | 0.6957 $\pm$ 0.2416 | [0.3113, 1.0000] | 0.8398 [0.7740, 0.8938] | Mean pass ($\ge 0.65$) |
| | UL30-Foscarnet | 0.7423 $\pm$ 0.1593 | [0.4888, 0.9958] | 0.8545 [0.7956, 0.9136] | Mean pass ($\ge 0.65$) |
| **ESM-2 35M** | UL23-Aciclovir | 0.6110 $\pm$ 0.2001 | [0.3625, 0.8595] | 0.6640 [0.5748, 0.7292] | Fail |
| *(Baseline)* | UL30-Aciclovir | 0.6795 $\pm$ 0.2441 | [0.2911, 1.0000] | 0.6681 [0.5814, 0.7479] | Fail |
| | UL30-Foscarnet | 0.7481 $\pm$ 0.1842 | [0.4550, 1.0000] | 0.7381 [0.6205, 0.8351] | Fail |
| **BLOSUM62** | UL23-Aciclovir | 0.6380 $\pm$ 0.0690 | [0.5524, 0.7237] | 0.6130 [0.5484, 0.6728] | Non-zero baseline |
| *(Radicalness)*| UL30-Aciclovir | 0.5494 $\pm$ 0.2496 | [0.1523, 0.9466] | 0.4942 [0.4115, 0.5705] | Near chance |
| | UL30-Foscarnet | 0.5074 $\pm$ 0.2997 | [0.0305, 0.9844] | 0.4531 [0.3515, 0.5560] | Pure chance |

ESM-2 650M is the only model that satisfies the pre-registered stopping rule on all three pairs ($\ge 0.65$ mean and lower bound of 95% $t$-CI $> 0.50$). ProtBERT achieves high mean AUROCs (0.696 to 0.753), but exhibits higher fold-to-fold variance across UL30 blocks. ESM-2 650M was designated the primary model post hoc based on this result.

---

### 3.3 Amino Acid Radicalness Separates the Two Resistance Targets

The statistically robust component of the biological dissociation between thymidine kinase and DNA polymerase is amino acid radicalness. BLOSUM62 radicalness — a single scalar per substitution requiring no alignment, no model training, and no position-level information — separates the two targets cleanly (Table 2, BLOSUM62 rows):

- **Thymidine Kinase (UL23-Aciclovir):** Contiguous-block mean AUROC **0.6380 $\pm$ 0.0690** (95% $t$-CI: [0.5524, 0.7237]). All five individual blocks exceed 0.55. The lower bound of the confidence interval clears chance (0.50). Because ACV is a prodrug requiring TK phosphorylation, resistance arises primarily through disruptive substitutions that abolish kinase activity at conserved catalytic residues — and simple amino acid radicalness captures this.
- **DNA Polymerase (UL30-Aciclovir):** Contiguous-block mean AUROC **0.5494 $\pm$ 0.2496** (95% $t$-CI: [0.1523, 0.9466]). **UL30-Foscarnet:** **0.5074 $\pm$ 0.2997** (95% $t$-CI: [0.0305, 0.9844]). Both are indistinguishable from chance. The viral polymerase cannot simply be inactivated — it must preserve catalytic replication while selectively discriminating against nucleoside or pyrophosphate inhibitors — and raw substitution radicalness carries no signal for this biologically constrained process.

This divergence is the paper's most interpretable finding and is independent of any PLM or alignment method.

#### Alignment Baselines and the Residual PLM Margin

Benchmarking PLMs against four alignment-derived baselines provides biological context for the radicalness dissociation (Table 3). On UL23, all four baselines cluster tightly at **0.7237 to 0.7373** and tie ESM-2 650M (**0.7397**; paired $t$-test vs. PSSM: $t=0.251, p=0.814$), confirming that conservation fully accounts for PLM performance on the kinase. On UL30, zero-shot PLMs outperform the strongest alignment baseline (positional Shannon entropy) by **$+0.0876$** (ACV) and **$+0.0563$** (FOS), but both paired per-block 95% Student's $t$ confidence intervals span zero (UL30-ACV $p=0.174$; UL30-FOS $p=0.373$, at $n=4$ evaluable blocks). This residual PLM margin is therefore **directional rather than statistically decisive**.

**Table 3: Evolutionary Baselines on UL23 and UL30.** Contiguous-block mean AUROC $\pm$ SD. Significance column reports the paired per-block Student's $t$-test of ESM-2 650M minus the strongest baseline on that row (underlined).

| Target Pair | Multi-Herpesvirus Weighted Conservation | Family Conservation Fraction | Negative Positional Entropy | PSSM Profile (Log-Odds) | ESM-2 650M | Paired $\Delta$ vs. Strongest ($p$-value) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **UL23-ACV** | <u>0.7373 $\pm$ 0.0920</u> | 0.7302 $\pm$ 0.0918 | 0.7251 $\pm$ 0.0892 | 0.7237 $\pm$ 0.1009 | **0.7397 $\pm$ 0.1565** | $+0.0024$ ($p=0.98$; **tied**) |
| **UL30-ACV** | 0.6428 $\pm$ 0.3174 | 0.6994 $\pm$ 0.2462 | <u>0.6999 $\pm$ 0.2601</u> | 0.6507 $\pm$ 0.1920 | **0.7875 $\pm$ 0.1772** | $+0.0876$ ($p=0.174$; **n.s.**) |
| **UL30-FOS** | 0.6604 $\pm$ 0.2664 | 0.8223 $\pm$ 0.0832 | <u>0.8246 $\pm$ 0.0817</u> | 0.7832 $\pm$ 0.0949 | **0.8809 $\pm$ 0.1250** | $+0.0563$ ($p=0.373$; **n.s.**) |

---

### 3.4 Supervised Predictor and Honest Paired Comparison
We trained supervised classifiers combining PLM representations, conservation metrics, and physicochemical deltas. Table 4 reports the paired per-block difference ($\Delta = \text{AUROC}_{\text{supervised}} - \text{AUROC}_{\text{ESM-2}}$) across held-out contiguous blocks.

**Table 4: Paired Comparison of Supervised Predictor versus Zero-Shot ESM-2 650M.**

| Target Pair | Supervised Contig-Block Mean $\pm$ SD | Zero-Shot ESM-2 Contig-Block Mean $\pm$ SD | Mean Paired Delta ($\Delta$) $\pm$ SD | 95% Student's $t$ CI on $\Delta$ | Paired $t$-test ($p$-value) | Strict C3 Rule Passed? |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **UL23-Aciclovir** | **0.8188 $\pm$ 0.1121** | 0.7397 $\pm$ 0.1565 | **$+0.0791 \pm 0.0886$** | [-0.0310, +0.1891] | $t=1.994$ ($p=0.117$) | False (CI spans 0) |
| **UL30-Aciclovir** | 0.7823 $\pm$ 0.1944 | 0.7875 $\pm$ 0.1772 | **$-0.0052 \pm 0.0235$** | [-0.0426, +0.0321] | $t=-0.446$ ($p=0.686$) | False (Exact tie) |
| **UL30-Foscarnet** | 0.8845 $\pm$ 0.1003 | 0.8809 $\pm$ 0.1250 | **$+0.0036 \pm 0.0426$** | [-0.0642, +0.0713] | $t=0.167$ ($p=0.878$) | False (Exact tie) |

Pre-registered Condition C3 specified that the supervised predictor must outperform zero-shot scores alone with the confidence interval excluding zero. **Condition C3 was not met.** 
- On UL30, the paired difference is negligible ($-0.0052$ and $+0.0036$), indicating that zero-shot PLM representations already saturate model capacity.
- On UL23, supervision produces a substantial $+0.0791$ mean gain across 4 of 5 blocks (lifting AUROC to 0.8188), reflecting the complementary value of conservation and radicalness features in TK.

#### Probabilistic Calibration
Out-of-fold probabilistic predictions demonstrated good calibration across all three pairs. Brier scores, bootstrap intervals, and null baselines are all computed over the **evaluable-block cohort** — the rows contributing to block AUROC averaging, excluding single-class Block 0 ($N=286$ for UL23-ACV, $269$ for UL30-ACV, $250$ for UL30-FOS). The null baseline is the Brier score of a constant prediction at that cohort's base rate, $\bar{y}(1 - \bar{y})$:
- **UL23-Aciclovir:** Brier score = **0.1777** (95% bootstrap CI: [0.1543, 0.2006]; null baseline: 0.2435; Brier skill score: **$+0.271$**; indicative ECE: 0.0722).
- **UL30-Aciclovir:** Brier score = **0.1136** (95% bootstrap CI: [0.0878, 0.1406]; null baseline: 0.1774; Brier skill score: **$+0.359$**; indicative ECE: 0.0676).
- **UL30-Foscarnet:** Brier score = **0.0841** (95% bootstrap CI: [0.0602, 0.1131]; null baseline: 0.1146; Brier skill score: **$+0.266$**; indicative ECE: 0.0648).

Because the evaluable cohort excludes an all-sensitive block, its base rate is higher than that of the full cohort, and the corresponding null baseline is therefore larger. For transparency we note that referencing the same Brier scores against full-cohort nulls (0.2435, 0.1519, 0.0949) would yield skill scores of $+0.271$, $+0.252$ and $+0.113$; the evaluable-cohort pairing reported above is the internally consistent one, since the Brier scores themselves are computed on that cohort. Both sets are logged in `logs/phase1_deliverables_results.json` under `deliverable_1_2_brier_bootstrap_audit`. ECE is reported as an indicative secondary metric only, being subject to binning instability at these sample sizes.

---

### 3.5 Prospective Temporal Validation (Pre-2016 to Post-2016)
To evaluate real-world utility, models trained on pre-2016 variants were evaluated on post-2016 variants ($N=307$), stratified by whether the variant occurred at a previously characterized sequence position (Table 5).

**Table 5: Prospective Temporal Holdout Stratified by Position Novelty.**

| Target Pair | Cohort Split | Test Variants (Resistant) | Supervised AUROC | Zero-Shot ESM-2 AUROC | Supervised AUPRC (Base Rate) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **UL23-Aciclovir** | All Post-2016 | 131 (60) | **0.7641** | 0.7250 | 0.7918 (0.4580) |
| | Seen Positions | 33 (19) | **0.8872** | 0.9023 | — |
| | **Novel Positions** | 98 (41) | **0.6996** | 0.6500 | — |
| **UL30-Aciclovir** | All Post-2016 | 92 (17) | **0.8290** | 0.7647 | 0.5386 (0.1848) |
| | Seen Positions* | 15 (2) | **0.9231** | 0.9615 | — |
| | **Novel Positions** | 77 (15) | **0.8075** | 0.7462 | — |
| **UL30-Foscarnet** | All Post-2016 | 84 (9) | **0.8533** | 0.8667 | 0.3917 (0.1071) |
| | Seen Positions* | 16 (2) | **0.9286** | 0.9286 | — |
| | **Novel Positions** | 68 (7) | **0.8126** | 0.8454 | — |

*\*AUROC at seen positions for UL30 is evaluated on only $N=2$ resistant variants (out of 15 and 16 variants, respectively) and is reported strictly for descriptive completeness; it is not statistically stable.*

The stratified temporal evaluation directly confirms the paper's central finding:
1. When prospective variants occur at **previously seen positions**, models achieve apparent AUROCs of **0.8872 to 0.9286**, reflecting position-level memorization.
2. When evaluated on **completely novel sequence positions**, performance lands at **0.6996** (UL23-ACV), **0.8075** (UL30-ACV), and **0.8126** (UL30-FOS).
3. These prospective novel-position AUROCs **are broadly consistent with our contiguous-block cross-validation estimates** (0.8188, 0.7823, 0.8845), supporting contiguous-block holdout as a reasonable surrogate for real-world prospective generalization.

---

### 3.6 Drug Conditioning and Negative Control
A multi-task model trained jointly across ACV and FOS achieved a contiguous-block mean AUROC of **0.8695 $\pm$ 0.0855** (95% $t$-CI: [0.7634, 0.9756]).

To test target specificity, we performed a negative control evaluation on UL23-Foscarnet ($N=189$ in the unified database). Because FOS binds directly to DNA polymerase and does not require activation by thymidine kinase, mutations in UL23 do not confer true biological FOS resistance (180/189 entries, 95.2%, are sensitive). The drug-conditioned model correctly suppressed false transfer, predicting $P(\text{Resistant} \mid \text{UL23}, \text{FOS}) \approx 0$. Because the negative control cohort is 95.2% sensitive, this confirms operational sanity rather than demonstrating a novel pathway routing mechanism.

---

## 4. Discussion

### Methodological Implications for Viral Variant Benchmarking
Our findings demonstrate that position-grouped cross-validation, while superior to random row splits, is insufficient to prevent shortcut learning in viral resistance benchmarks. In proteins with compact active sites, holding out individual positions allows machine learning models to exploit spatial correlation across adjacent residues. Our contiguous-block holdout protocol establishes a rigorous standard for benchmarking resistance models against true regional generalization.

### Structural Biology and Mechanistic Interpretations
Recent cryo-EM structural determinations of the HSV-1 pre-insertion polymerase complex (Shankar et al., *Cell* 2024; PMID 39197451) revealed the atomic basis of substrate discrimination and drug coordination within the Palm and Finger domains. Our finding that zero-shot PLMs achieve high predictive power on UL30 without supervision is consistent with, but does not directly test, the structural remodeling mechanism proposed by Shankar et al. (Cell 2024). Unlike thymidine kinase, where resistance is frequently synonymous with complete loss of catalytic function, polymerase resistance requires delicate conformational remodeling that allows continued viral genome replication while discriminating against nucleoside or pyrophosphate inhibitors. While PLM representations implicitly encode these biophysical constraints, direct physical-chemical simulation and 3D structural modeling remain important directions for future work.

### The Radicalness Dissociation and the Narrowing PLM Margin
The statistically robust finding in our biological analysis is the radicalness divergence: BLOSUM62 amino acid radicalness — a zero-parameter, alignment-free feature — predicts resistance on UL23 (0.6380, 95% CI [0.5524, 0.7237]) but operates at chance on UL30 (0.5074–0.5494). This dissociation reflects the distinct resistance mechanisms of the two enzymes and holds under contiguous-block holdout with all individual UL23 blocks above 0.55.

An important methodological pattern complements this finding: the behavior of the residual PLM advantage on UL30 as the alignment baseline is strengthened. Against amino acid radicalness (chance on polymerase), language models demonstrate a large advantage. Against multi-herpesvirus weighted conservation, the advantage is +0.145 (ACV) and +0.221 (FOS). Against the strongest alignment baseline we evaluated — negative positional Shannon entropy — the advantage contracts to **+0.0876** on ACV (0.7875 vs 0.6999) and **+0.0563** on FOS (0.8809 vs 0.8246).

We state this trend plainly because it bears directly on how the result should be read: **with each stronger alignment baseline, the UL30 margin shrinks, and it shrinks most on the target where the margin was largest.** The FOS advantage falls by roughly three-quarters, from +0.221 to +0.056, once the comparator is positional entropy rather than the weighted conservation score. Notably, the weighted conservation score—which combines entropy with a BLOSUM62 radicalness term—is *weaker* on UL30 than entropy alone, because the radicalness component is uninformative on polymerase and dilutes the entropy signal. Nominating it as the single conservation comparator would therefore overstate the PLM advantage, and we avoid doing so.

Two honest conclusions follow. First, the qualitative dissociation between UL23 and UL30 is robust across multiple angles: on UL23 every alignment baseline ties ESM-2, while on UL30 radicalness is at chance and PLMs remain directionally ahead of every alignment baseline tested. Second, the *magnitude* of the UL30 PLM advantage is considerably smaller than a comparison against weighted conservation alone implies. When tested formally via a paired per-block Student's $t$-test against positional entropy, the mean paired differences are **$\Delta = +0.0877 \pm 0.0988$** on UL30-ACV ($t(3)=1.775$, $p=0.174$, 95% CI: $[-0.0696, +0.2449]$) and **$\Delta = +0.0563 \pm 0.1079$** on UL30-FOS ($t(3)=1.044$, $p=0.373$, 95% CI: $[-0.1154, +0.2281]$). Because both 95% confidence intervals span zero, the PLM advantage over pure positional conservation on polymerase is **directional and qualitative rather than statistically decisive**. The pattern is consistent across three independent lines of evidence (conservation gap, radicalness collapse, and lack of incremental supervision gain), but no single comparison is individually significant at this block sample size ($n=4$).

To directly test whether a profile-based model captures this residual margin, we evaluated a multi-herpesvirus position-specific scoring matrix (PSSM) baseline with Dirichlet pseudocount smoothing ($\alpha=1.0$, Robinson & Robinson background frequencies). The PSSM achieves contiguous-block mean AUROCs of **0.7237 $\pm$ 0.1009** (UL23-ACV), **0.6507 $\pm$ 0.1920** (UL30-ACV), and **0.7832 $\pm$ 0.0949** (UL30-FOS). On UL23, the PSSM ties zero-shot ESM-2 ($0.7237$ vs. $0.7397$, paired difference $+0.0160 \pm 0.1424$, $p=0.814$). On UL30, ESM-2 maintains a positive paired difference over the PSSM ($+0.1369 \pm 0.1177$ on ACV; $+0.0977 \pm 0.1154$ on FOS), but with $n=4$ evaluable blocks, both 95% Student's $t$ confidence intervals span zero ($p=0.103$ and $p=0.189$, respectively). Crucially, the PSSM performs below positional Shannon entropy on both polymerase targets (0.6507 vs. 0.6999 on ACV; 0.7832 vs. 0.8246 on FOS), leaving Shannon entropy as the strongest alignment comparator. Thus, evaluating a profile model confirms rather than contradicts our conclusion: zero-shot PLMs retain a directional advantage on polymerase, but the margin does not achieve statistical significance under domain holdout. Deeper alignment methods with explicit epistatic contact models (such as GEMME or EVE) and viral-specific fine-tuning (e.g., Sawhney et al., *PeerJ* 2025, PMID 40985030, which established that fine-tuning mitigates the underrepresentation of viral proteins in standard PLM training corpora) remain open directions for future exploration.

### Clinical Scope: The Truncating Class and the Missense Diagnostic Bottleneck
A critical scope consideration for translational genomics is that our benchmark evaluates single-site missense mutations. In our unified clinical database, truncating and frameshift mutations account for 44.4% (186/419) of cataloged resistance-conferring variants in thymidine kinase (UL23). Because premature truncation abolishes catalytic activity and nucleoside phosphorylation, clinical interpretation of nonsense and frameshift variants is deterministic—a clinician presented with an early stop codon does not require machine learning prediction to anticipate aciclovir resistance. The genuine clinical and computational bottleneck lies in interpreting novel missense variants, which represent 55.6% of UL23 and virtually all UL30 resistance, where substitutions can either represent benign polymorphisms or subtle resistance-conferring mutations. By focusing strictly on single-site missense substitutions, this benchmark isolates the specific class where algorithmic guidance is clinically essential.

### Study Limitations
1. **Scope Limitation:** Our benchmark evaluates single-site missense substitutions. While missense mutations represent 55.6% of UL23 and nearly 100% of UL30 resistance, roughly 44.4% of clinical UL23 resistance arises via frameshift or nonsense truncations.
2. **Sample Size in Block Splits:** Dividing sequences into 5 spatial blocks limits the effective sample size for Student's $t$ confidence intervals ($n=4$ or $5$), leading to relatively wide interval estimates despite strong point estimates.
3. **Temporal Covariate Shift:** Our domain classifier confirmed significant covariate shift between pre-2016 and post-2016 literature (AUROC 0.6794), reflecting evolving sequencing methodologies and clinical sampling practices.

---

## 5. Declarations and Data Availability

- **Data Availability:** The primary evaluation cohort (`primary_position_folds.tsv`, SHA-256 `7C4D0A3DE2CD171B07230757F44E0DFBD61C6A0A8570DC8F263195B354C66909`) and scored model outputs (`primary_position_folds_scored.tsv`) are permanently archived in the project repository.
- **Code Availability:** All model training, evaluation scripts, and cross-validation runners are available at: `https://github.com/shubhamraikwar142/antiviral-resistance-hsv-1`.
- **Pre-Submission Human Audit Gate:** Submission remains conditional on completion of the 50-row dual-annotator kappa curation (`logs/kappa_sample_50rows.tsv`), threshold citation review (`docs/threshold_citation_checklist.md`), and formal PI sign-off (`docs/pi_signoff_form.md`).

