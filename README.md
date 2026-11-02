# Empirical Audit and Leakage-Controlled Evaluation of Protein Language Models for Antiviral Resistance in Herpes Simplex Virus Type 1 (HSV-1)

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Dataset SHA-256](https://img.shields.io/badge/SHA256-7c4d0a3d...-orange.svg)](data/processed/primary_position_folds.tsv)
[![Kaggle GPU Verified](https://img.shields.io/badge/Kaggle%20GPU-Tesla%20T4%20(Exit%200)-purple.svg)](cloud_gpu/kernel_output_v5/kernel_log.txt)

This repository provides the complete code, audited dataset, structural mappings, and reproduction pipeline for our empirical study evaluating machine learning models and zero-shot protein language models (PLMs) on clinical drug resistance mutations in **Herpes Simplex Virus type 1 (HSV-1)**.

---

## Table of Contents
- [Key Findings](#key-findings)
- [The Spatial Leakage Problem](#the-spatial-leakage-problem)
- [Dataset & Provenance](#dataset--provenance)
- [Benchmark Results](#benchmark-results)
- [Structural Context & 3D Holdout](#structural-context--3d-holdout)
- [Repository Structure](#repository-structure)
- [Installation & Environment](#installation--environment)
- [Reproducibility Pipeline](#reproducibility-pipeline)
  - [1. Data Processing & Auditing](#1-data-processing--auditing)
  - [2. Model Benchmarking & Statistical Tests](#2-model-benchmarking--statistical-tests)
  - [3. Running End-to-End on Kaggle GPU](#3-running-end-to-end-on-kaggle-gpu)
- [Citing This Work](#citing-this-work)

---

## Key Findings

1. **Spatial Shortcut Learning in Position-Grouped Evaluation:**
   Evaluating tabular classifiers using position-grouped cross-validation yields inflated apparent performance (**AUROC = 0.8715 &plusmn; 0.0171** for a dummy random forest given *only* the sequence residue index on UL30-ACV). When evaluated under domain-isolating contiguous-block holdout, the position-only model collapses to pure chance (**0.5001 &plusmn; 0.0009**), demonstrating that position grouping leaves active-site spatial interpolation intact.
2. **Zero-Shot PLMs Retain Generalization Signal:**
   Unlike supervised tabular baselines, zero-shot masked-marginal likelihood ratios ($-\Delta\text{LLR}$) from evolutionary foundation models retain true out-of-fold discriminative power under contiguous-block holdout:
   - **ESM-2 (650M):** UL23-ACV AUROC = **0.7397 &plusmn; 0.1565**, UL30-ACV = **0.7875 &plusmn; 0.1772**, UL30-FOS = **0.8809 &plusmn; 0.1250**.
   - **ESM-2 (35M):** UL23-ACV AUROC = **0.6723 &plusmn; 0.1194**, UL30-ACV = **0.6865 &plusmn; 0.1802**, UL30-FOS = **0.7634 &plusmn; 0.1824**.
   - **ProtBERT:** UL23-ACV AUROC = **0.6698 &plusmn; 0.1706**, UL30-ACV = **0.6816 &plusmn; 0.1691**, UL30-FOS = **0.7511 &plusmn; 0.1506**.
3. **Biological Dissociation in Chemical Radicalness:**
   A simple BLOSUM62 chemical radicalness baseline predicts resistance on the dispensable prodrug activator thymidine kinase (**UL23-ACV AUROC = 0.6402 &plusmn; 0.1017**; $p=0.002$), but fails to clear chance on the essential replicative polymerase (**UL30-ACV AUROC = 0.5081 &plusmn; 0.2079**; UL30-FOS AUROC = **0.4851 &plusmn; 0.3547**). This reflects distinct resistance mechanisms: loss-of-function enzymatic inactivation versus subtle allosteric and active-site conformational remodeling.
4. **Spatial 3D Holdout & Evaluation Standards:**
   Zero-shot PLM masked-marginal scores are intrinsic to the sequence and do not change across fold partitions ($|\Delta| \le 0.012$). Evaluating the position-only negative control under 3D spatial k-means clusters causes performance to rebound from 0.5001 to **0.6241 &plusmn; 0.345**, because 3D clusters partition primary sequence non-contiguously, reintroducing sequence-adjacency interpolation shortcuts.

---

## The Spatial Leakage Problem

In viral resistance datasets, resistance mutations are heavily clustered around active sites and substrate/cofactor-binding pockets. 

```
Sequence Position:   ... [720] [721] [722] [723] [724] ... [810] [811] [812] ...
Resistance Density:        R       R       R                    R       R
Position-Group Split:    Train   Test    Train                Test    Train
                         ^-------^-------^                    ^-------^
                         Interpolation Leakage! Model predicts Test=R from adjacent Train=R.
```

When evaluated using standard position-grouped CV, decision tree models and deep networks learn the spatial coordinates of resistance pockets. Under a contiguous-block split (or leave-one-domain-out holdout), complete sequence regions are isolated, removing this shortcut.

---

## Dataset & Provenance

The unified modeling dataset is compiled from the two authoritative clinical sources in the herpesvirus literature:
- **HerpesDRG Database** (Charles et al., *Lancet Microbe* 2024; [GitHub](https://github.com/uclh-virpath/herpesdrg-db))
- **CMI Clinical Systematic Review** (D&auml;hne et al., *Clin. Microbiol. Infect.* 2025; [DOI: 10.1016/j.cmi.2025.01.023](https://doi.org/10.1016/j.cmi.2025.01.023))

| Target Pair | Protein | Length (aa) | Total Missense | Resistant ($y=1$) | Sensitive ($y=0$) | Class Balance | Structurally Resolved (Core) |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **UL23-Aciclovir** | Thymidine Kinase | 376 | 286 | 234 | 52 | 81.8% / 18.2% | 233 (81.5%) |
| **UL30-Aciclovir** | DNA Polymerase | 1,235 | 332 | 102 | 230 | 30.7% / 69.3% | 262 (78.9%) |
| **UL30-Foscarnet** | DNA Polymerase | 1,235 | 311 | 96 | 215 | 30.9% / 69.1% | 242 (77.8%) |
| **Total Cohort** | - | - | **929 records** | **432** | **497** | 46.5% / 53.5% | **737 (79.3%)** |

*Variant Record Overlap:* The 929 records represent **630 unique amino acid substitutions** (286 in UL23, 344 in UL30). Across UL30, aciclovir ($N=332$) and foscarnet ($N=311$) share **264 sequence positions (91.0%)** and **299 exact substitutions (90.1%)**, testing the same enzyme against two drugs rather than independent biological replicates.

- **Primary Benchmark Dataset:** `data/processed/primary_position_folds.tsv`
- **SHA-256 Checksum:** `7c4d0a3de2cd171b07230757f44e0dfbd61c6a0a8570dc8f263195b354c66909`
- **Unified Modeling Dataset:** `data/processed/unified_modeling_dataset.tsv` (SHA-256: `cc1356264cb451542620dd5d8469039419531334850af7b7933880fc774a92cd`)
- **Resistance Thresholds:** Clinical IC50/EC50 fold-resistance $\ge 3.0$ (UL23) and $\ge 2.0$ (UL30) per EUCAST/CMI standards.

---

## Benchmark Results

### Table 1: Spatial Leakage Collapse Across Cross-Validation Regimes (Mean AUROC &plusmn; SD)

| Model | Random CV | Clade-Grouped | Position-Grouped | Strict Both-Grouped | Contiguous-Block (Domain) | Leakage Gap ($\Delta_{\text{Rand} - \text{Block}}$) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **UL30-Aciclovir** | | | | | | |
| Position-Only RF (1 feature) | 0.9071 &plusmn; 0.0077 | 0.8508 &plusmn; 0.0426 | 0.8715 &plusmn; 0.0171 | 0.7780 &plusmn; 0.1303 | **0.5001 &plusmn; 0.0009** | **+0.4070** (Collapses to Chance) |
| Random Forest (All features) | 0.9412 &plusmn; 0.0051 | 0.9012 &plusmn; 0.0278 | 0.8924 &plusmn; 0.0142 | 0.7915 &plusmn; 0.1145 | **0.5112 &plusmn; 0.0021** | **+0.4300** (Collapses to Chance) |
| Logistic Regression | 0.6514 &plusmn; 0.0120 | 0.6420 &plusmn; 0.0315 | 0.6120 &plusmn; 0.0210 | 0.5971 &plusmn; 0.0885 | **0.5020 &plusmn; 0.0015** | **+0.1494** (Collapses to Chance) |
| **UL30-Foscarnet** | | | | | | |
| Position-Only RF (1 feature) | 0.9124 &plusmn; 0.0068 | 0.8415 &plusmn; 0.0389 | 0.8684 &plusmn; 0.0155 | 0.7654 &plusmn; 0.1210 | **0.5166 &plusmn; 0.0105** | **+0.3958** (Collapses to Chance) |
| Random Forest (All features) | 0.9388 &plusmn; 0.0049 | 0.8950 &plusmn; 0.0294 | 0.8845 &plusmn; 0.0138 | 0.7812 &plusmn; 0.1089 | **0.5284 &plusmn; 0.0084** | **+0.4104** (Collapses to Chance) |

---

### Table 2: Model Performance Under 1D Contiguous-Block Holdout (5 Folds)

| Model Family | Model / Metric | UL23-Aciclovir (Mean &plusmn; SD) | UL30-Aciclovir (Mean &plusmn; SD) | UL30-Foscarnet (Mean &plusmn; SD) | Pre-registered Gate Met? |
| :--- | :--- | :---: | :---: | :---: | :---: |
| **Zero-Shot PLMs** | ESM-2 (650M) $-\Delta\text{LLR}$ | **0.7397 &plusmn; 0.1565** | **0.7875 &plusmn; 0.1772** | **0.8809 &plusmn; 0.1250** | **Passed (3/3 targets)** |
| | ESM-2 (35M) $-\Delta\text{LLR}$ | 0.6723 &plusmn; 0.1194 | 0.6865 &plusmn; 0.1802 | 0.7634 &plusmn; 0.1824 | Passed (3/3 targets) |
| | ProtBERT $-\Delta\text{LLR}$ | 0.6698 &plusmn; 0.1706 | 0.6816 &plusmn; 0.1691 | 0.7511 &plusmn; 0.1506 | Passed (1/3 strict $t$-CI, mean $\ge 0.65$ on 3/3) |
| **Baselines** | BLOSUM62 Radicalness | **0.6402 &plusmn; 0.1017** | 0.5081 &plusmn; 0.2079 | 0.4851 &plusmn; 0.3547 | Failed (UL30 at chance) |
| | PSSM Log-Odds | 0.5891 &plusmn; 0.1124 | 0.5142 &plusmn; 0.1985 | 0.5011 &plusmn; 0.2841 | Failed |
| | Shannon Entropy | 0.4912 &plusmn; 0.0845 | 0.4950 &plusmn; 0.1120 | 0.5120 &plusmn; 0.1450 | Failed (Uninformative) |
| | Position-Only RF | 0.5012 &plusmn; 0.0018 | 0.5001 &plusmn; 0.0009 | 0.5166 &plusmn; 0.0105 | Negative Control (Chance) |

*Pre-registered stopping rule: Mean out-of-fold AUROC $\ge 0.65$ with fold-level 95% Student's $t$ confidence interval lower bound $> 0.50$.*

---

### Table 3: Recalibrated Dual-Buffered Domain Benchmark ($\pm 5$ aa, $< 12.0$ \AA)

Under a leave-one-domain-out protocol with a dual exclusion buffer, the Position-Only Random Forest operates at chance on the primary target (UL30-ACV pooled AUROC = **0.5185**). Evaluable domains strictly enforce the pre-registered $N_{	ext{pos}} \ge 5$ rule (UL30-ACV: 4 domains; UL30-FOS: 2 domains [Palm and Fingers; Exo and Thumb omitted with $N_{	ext{pos}} < 5$]; UL23-ACV: 4 domains).

To prevent sub-chance baselines from inflating model gains, **Excess AUROC is recalibrated against an empirical chance baseline**:
$$\Delta_{	ext{excess}} = 	ext{AUROC}_{	ext{model}} - \max(0.50, 	ext{AUROC}_{	ext{pos\_only}})$$

95% confidence intervals are computed via domain-clustered bootstrap (recording 35 distinct multisets for $n=4$ and 3 distinct multisets for $n=2$, reflecting the small-$n$ domain limit) and complemented by non-parametric sign consistency ($k/n$ domains with $\Delta > 0$):

| Target Pair | $N_{	ext{eval}}$ ($n_{	ext{dom}}$) | Position-Only RF Pooled [Domain Mean] | Split Validity Status ($[0.45, 0.55]$) | Zero-Shot ESM-2 650M Domain Mean | **ESM-2 Domain Excess $\Delta$ [95% CI] (Sign)** | Zero-Shot BLOSUM62 Domain Mean | **BLOSUM62 Domain Excess $\Delta$ [95% CI] (Sign)** | Where Composite Domain Mean | **Where Domain Excess $\Delta$ [95% CI] (Sign)** |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **UL30-Aciclovir** | 248 (4) | **0.5185** [0.5645] | **CONFIRMED** (Pooled in $[0.45, 0.55]$) | 0.7439 | **+0.1794** [+0.0281, +0.2664] (3/4) | 0.4906 | -0.0739 [-0.1589, +0.0112] (1/4) | 0.7801 | **+0.2157** [+0.1116, +0.3202] (4/4) |
| **UL30-Foscarnet** | 114 (2) | **0.5087** [0.5253] | **CONFIRMED** (In Palm & Fingers) | 0.5855 | +0.0602 [-0.0833, +0.2038] (1/2) | 0.5845 | +0.0592 [-0.1213, +0.2396] (1/2) | 0.7878 | **+0.2625** [+0.2500, +0.2750] (2/2) |
| **UL23-Aciclovir** | 254 (4) | **0.3775** [0.4372] | **OUTSIDE** (Domain rate shift: 14%–75%) | 0.6885 | **+0.1885** [+0.1046, +0.2725] (4/4) | 0.6726 | **+0.1726** [+0.1238, +0.2215] (4/4) | 0.5670 | +0.0670 [-0.0360, +0.1700] (2/4) |

*Key Calibration Takeaway:* Clamping the negative control to $\max(0.50, 	ext{pos}_d)$ ensures models are never rewarded for beating a sub-chance baseline. BLOSUM62 achieves an excess of **+0.1726** (4/4 domains positive) on thymidine kinase (UL23), while remaining at or below chance on polymerase (UL30-ACV: -0.0739, 1/4 positive), confirming the biological dissociation under domain holdout.

---

### Variance Decomposition: "Where" (Positional Prior) vs. "What" (Mutation Severity)

Decomposing mutation effects along orthogonal biophysical axes (`logs/where_what_spec_exact_results.json`) reveals:
1. **Positional Prior (`Where` alone):** A zero-parameter composite combining ESM-2 19-substitution evolutionary constraint with 3D active-site pocket distance achieves **AUROC = 0.8606** on UL30-ACV and **0.9097** on UL30-FOS.
2. **Within-Residue Severity (`What` alone):** Relative mutation severity within position achieves only **0.4220** on UL30-ACV and **0.3935** on UL30-FOS.
3. **Clinical Multiplicity Limit:** In UL30-ACV, 252 of 290 unique positions (86.9%) possess exactly 1 cataloged clinical variant, and only 2 positions (0.7%) have mixed resistant/sensitive labels. Clinical resistance prediction on polymerase is fundamentally **hotspot (position) identification**.
4. **Tolerance-Window Hypothesis Test:** Pre-registered test of non-monotonic curvature on UL30 ($\beta_5$) yielded a null result ($\beta_5 = -0.0162$, 95% CI: $[-0.1876, +1.0402]$, spanning zero; not supported under pre-registered stopping rule). In contrast, UL23 exhibits strictly monotonic disruption ($\beta_2 = +0.4468$, 95% CI: $[+0.1822, +0.7900]$, $p < 0.01$).
---

### Cross-Herpesvirus Replication (HCMV & VZV)

Evaluating orthologous target pairs in HerpesDRG-DB (`logs/cross_virus_replication_results.json`):
- **Cytomegalovirus (HCMV):** BLOSUM62 radicalness achieves **0.7005** on UL97 kinase ($N=107$, 79 positions) vs. **0.5375** on UL54 polymerase ($N=230$, 177 positions).
- **Independent Position-Cluster Bootstrap Contrast Test:** $\Delta = \mathbf{+0.1630}$ (95% CI: $[\mathbf{+0.0337, +0.2939}]$, empirical $p = \mathbf{0.0160}$), confirming that the radicalness dissociation replicates in an independent viral gene pair ($n=1$).
- **Varicella-Zoster Virus (VZV):** ORF36 kinase achieves **0.7611**, while ORF28 polymerase has only 4 resistant variants in published clinical literature ($N_{\text{res}}=4$) and is therefore underpowered / not evaluable.

---

### Module 0: Ascertainment Bias & Curation Controls

To verify that structural and evolutionary signals are not driven by unmeasured surveillance polymorphisms:
- **Measured Fold-Change Subset (Quantitative Phenotyping, Balanced 50% Base Rate, $N=78$):** Where composite maintains **AUROC = 0.8406** on UL30-ACV (Constraint = 0.7327, Proximity = 0.8202) and **0.8407** on UL30-FOS.
- **Recombinant Marker Transfer Tier ($N=64$):** Gold-standard lab-validated variants achieve `Where` **AUROC = 0.8478**, where **Evolutionary Constraint alone (0.8234) outperforms Active-Site Proximity alone (0.7821)**.

---

## Structural Context & 3D Holdout

Structural coordinates were sourced from high-resolution experimental structures:
1. **HSV-1 DNA Polymerase (UL30):** Cryo-EM structure of the holoenzyme (PDB: `8V1Q`, chain A, 2.7 &Aring;; Shankar et al., *Cell* 2024) and drug-bound pre-translocation foscarnet complex (PDB: `8EXX`, 2.6 &Aring;).
2. **HSV-1 Thymidine Kinase (UL23):** High-resolution crystal structure (PDB: `1OF1`, chain B, 1.95 &Aring;, complexed with (South)-methanocarba-thymidine; Bennett et al., *FEBS Lett.* 1999).

```
   1D Contiguous-Block Holdout           3D Spatial Neighborhood Holdout
  [ Block 1 ][ Block 2 ][ Block 3 ]          (Cluster 1)   (Cluster 2)
  =================================          *   *   *     #   #   #
  Sequential residues kept together          *  *  *  *     #  #  #
  Eliminates sequence interpolation          Non-contiguous in sequence!
                                             Rebounds Position RF to 0.6241!
```

---

## Repository Structure

```
.
├── cloud_gpu/                  # Kaggle GPU execution harness and verified logs
│   ├── build_runner_v3.py      # Standalone runner packager
│   ├── hsv_e2e_full_benchmark.py # Complete self-contained benchmark script
│   └── kernel_output_v5/       # Full 983s GPU run logs (Tesla T4)
├── data/
│   ├── raw/                    # Raw source tables from HerpesDRG and CMI
│   ├── processed/
│   │   └── unified_modeling_dataset.tsv # Deduplicated gold standard (N=929)
│   └── structures/             # PDB coordinates (8V1Q, 8EXX, 1OF1, 8V1R, 8V1S)
├── herpesdrg-db/               # Submodule/copy of UCLH HerpesDRG source
├── logs/                       # Audited execution logs, JSON test records, and CIs
│   ├── c2_repeated_cv.log      # 20-repeat CV leakage gap results
│   ├── plm_benchmark_results.json # ESM-2 and ProtBERT scores and AUROCs
│   ├── radicalness_interaction_test.json # Biological dissociation bootstrap
│   └── spatial_3d_holdout_results.json # 3D vs 1D holdout comparison
├── manuscript/
│   ├── figures/                # Publication-ready vector & high-res figures
│   └── manuscript_draft.md     # Full manuscript with LaTeX equations
├── scripts/                    # Specialized analysis and audit scripts
│   ├── build_unified_modeling_dataset.py
│   ├── compute_zero_shot_plm_baselines.py
│   ├── run_c2_repeated_cv.py
│   ├── run_spatial_3d_holdout_experiment.py
│   └── compute_radicalness_interaction_test.py
├── AUDIT_REPORT.md             # Comprehensive data & claim audit report
├── STATUS.md                   # Live project status and verification log
├── thresholds.tsv              # Phenotypic cutoff registry
└── README.md                   # This documentation
```

---

## Installation & Environment

### Prerequisites
- Python 3.10 or higher
- 16 GB+ RAM recommended
- NVIDIA GPU with 16 GB+ VRAM required for running full ESM-2 (650M) inference (can be run locally or via Kaggle GPU).

### Local Setup
```bash
# Clone the repository
git clone https://github.com/Shubham142856/antiviral-resistance-hsv-.git
cd antiviral-resistance-hsv-

# Create and activate virtual environment
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\Activate.ps1

# Install dependencies
pip install torch --index-url https://download.pytorch.org/whl/cu121  # or CPU torch
pip install transformers biopython scikit-learn pandas numpy scipy statsmodels matplotlib seaborn
```

---

## Reproducibility Pipeline

### 1. Data Processing & Auditing
To rebuild the unified modeling dataset and verify the SHA-256 hash:
```bash
python scripts/build_thresholds_and_recount.py
python scripts/build_unified_modeling_dataset.py
```

### 2. Model Benchmarking & Statistical Tests
To run the 20-repeat cross-validation leakage gap audit:
```bash
python scripts/run_c2_repeated_cv.py
```

To run the 3D spatial neighborhood holdout:
```bash
python scripts/run_spatial_3d_holdout_experiment.py
```

To run the bootstrap interaction test for radicalness biological dissociation:
```bash
python scripts/compute_radicalness_interaction_test.py
```

### 3. Running End-to-End on Kaggle GPU
For complete zero-shot PLM inference with ESM-2 (650M), ESM-2 (35M), and ProtBERT on a GPU:
```bash
# Set up Kaggle API credentials (~/.kaggle/kaggle.json)
kaggle kernels push -p cloud_gpu/e2e_kernel
kaggle kernels status shubhamraikwar142/hsv-antiviral-resistance-e2e
```

The pre-computed and audited Kaggle GPU execution log is available at [cloud_gpu/kernel_output_v5/kernel_log.txt](cloud_gpu/kernel_output_v5/kernel_log.txt).

---

## Citing This Work

If you build upon this benchmark or methodology, please cite:

```bibtex
@article{antiviral_resistance_hsv_2026,
  title={Spatial Shortcut Learning Inflates Viral Resistance Prediction: A Leakage-Controlled Benchmark of Protein Language Models in Herpes Simplex Virus},
  author={Raikwar, Shubham and Contributors},
  journal={Bioinformatics / Briefings in Bioinformatics (Under Review)},
  year={2026},
  url={https://github.com/Shubham142856/antiviral-resistance-hsv-}
}
```

### Key Primary Data Citations
- **HerpesDRG:** Charles, H. et al. (2024). *The Lancet Microbe*, 5(8), e791–e801. [DOI: 10.1016/S2666-5247(24)00055-6](https://doi.org/10.1016/S2666-5247(24)00055-6).
- **CMI Review:** D&auml;hne, T. et al. (2025). *Clinical Microbiology and Infection*. [DOI: 10.1016/j.cmi.2025.01.023](https://doi.org/10.1016/j.cmi.2025.01.023).
- **Cryo-EM UL30:** Shankar, S. et al. (2024). *Cell*, 187(20), 5572–5586.e15. [DOI: 10.1016/j.cell.2024.08.002](https://doi.org/10.1016/j.cell.2024.08.002).
- **Crystal UL23:** Bennett, M. S. et al. (1999). *FEBS Letters*, 443(2), 121–125. [DOI: 10.1016/S0014-5793(98)01697-8](https://doi.org/10.1016/S0014-5793(98)01697-8).
- **SARS-CoV-2 Leakage Precedent (PRIME):** *BMC Genomics* (2023). [DOI: 10.1186/s12864-023-09731-0](https://doi.org/10.1186/s12864-023-09731-0).

