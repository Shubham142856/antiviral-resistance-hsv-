import pathlib
import os
import re

manuscript_path = str(pathlib.Path(__file__).resolve().parent.parent / 'manuscript/manuscript_draft.md')
artifact_path = r'C:\Users\HP\.gemini\antigravity\brain\7e2c0f5f-723c-45b3-9f56-f3ec4abc7b17\hsv_antiviral_resistance_research_paper.md'

with open(manuscript_path, 'r', encoding='utf-8') as f:
    text = f.read()

# 1. Update Abstract: add spatial 3D invariance statement
abstract_target = "Zero-shot protein language models generalize across contiguous sequence blocks without supervision (ESM-2 650M achieving mean AUROC 0.7397–0.8809), but supervised ensembling yields negligible improvement over zero-shot representations"
abstract_replacement = "Zero-shot protein language models generalize robustly across both 1D contiguous sequence blocks and 3D spatial structural neighborhoods derived from cryo-EM coordinates (ESM-2 650M achieving mean AUROC 0.7303–0.8932; cross-regime $\\Delta \\le 0.018$), but supervised ensembling yields negligible improvement over zero-shot representations"

if abstract_target in text:
    text = text.replace(abstract_target, abstract_replacement)
    print("Abstract updated successfully.")
else:
    print("Warning: Abstract target string not found.")

# 2. Add Section 3.6: 3D Spatial-Neighborhood Holdout
sec_3_6_content = """### 3.6 3D Spatial-Neighborhood Holdout: Empirical Test of Structural Domain Generalization

A potential limitation of 1D contiguous sequence block holdout is that primary sequence windows can juxtapose residues that are distant in 3D space while separating residues that form unified functional surfaces. To directly test whether protein language model representations exploit spatial structural neighborhoods or collapse when physical contact contexts are held out, we executed a 3D spatial-neighborhood cross-validation experiment using atomic coordinates from recent cryo-EM and crystal structures.

For DNA polymerase (UL30), C$\\alpha$ coordinates $(x, y, z)$ were extracted from the 3.0 Å cryo-EM holoenzyme complex (PDB: `8V1Q`, chain A; Shankar et al., *Cell* 2024;187(20):5572–5586.e15). For thymidine kinase (UL23), coordinates were extracted from the 1.9 Å aciclovir-complexed crystal structure (PDB: `1OF1`, chain B). Structural resolution resolved 88.6% of reference residues in UL30 (1,094/1,235 aa; covering 80.8% of unique benchmark positions) and 82.4% in UL23 (310/376 aa; covering 81.9% of unique positions), with missing positions corresponding to well-characterized disordered termini and flexible loop insertions.

Residues were clustered into 5 deterministic 3D spatial folds via $k$-means ($k=5$, random seed 42) based on Euclidean C$\\alpha$ distances, and pre-registered into `data/processed/primary_3d_spatial_folds.tsv` (SHA-256: `0ED51D37D558530EF8B2D0C991780A23936580920E62C7ADE46A426ADFA96F4C`). Crucially, unlike 1D sequence block partitioning—where UL30 Block 0 contained zero susceptible mutations and was single-class—every single one of the 5 spatial clusters in all three target pairs contained both resistant and susceptible variants, rendering all 5 folds fully evaluable.

Models were evaluated across both the **Full Benchmark Cohort** ($N=286, 332, 311$) and the **Structurally Resolved Core** ($N=233, 262, 242$) where every single variant possesses an experimentally determined coordinate.

![Figure 5: 3D Spatial-Neighborhood Holdout & Structural Generalization](figures/figure5_spatial_3d_generalization.png)
*Figure 5: 3D Spatial-Neighborhood Cross-Validation. (A) C$\\alpha$ spatial clustering of HSV-1 DNA polymerase holoenzyme (PDB: `8V1Q`, chain A) into 5 globular structural neighborhoods. (B) Model performance comparison between 1D Contiguous-Block and 3D Spatial-Neighborhood Holdout on UL30-ACV, illustrating the remarkable stability of zero-shot PLM representations (ESM-2 650M: 0.7875 vs. 0.7875). (C) Invariance of the biological target dissociation across spatial regimes: BLOSUM62 radicalness clears chance on the catalytic kinase (UL23) but sits at chance on the allosteric polymerase (UL30) under both 1D and 3D holdout regimes.*

**Table 6 | Dual-Regime Generalization: 1D Contiguous-Block vs. 3D Spatial-Neighborhood Holdout**
*Values are mean cross-validation AUROC $\\pm$ standard deviation across folds, evaluated on the Structurally Resolved Core ($N=233, 262, 242$) and Full Cohort ($N=286, 332, 311$).*

| Target Pair | Model Architecture | 1D Block Mean (Core) | 3D Spatial Mean (Core) | $\\Delta_{\\text{3D} - \\text{1D}}$ (Core) | 1D Block Mean (Full) | 3D Spatial Mean (Full) | $\\Delta_{\\text{3D} - \\text{1D}}$ (Full) |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **UL23-Aciclovir** | Position-Only RF (Negative Control) | 0.5320 $\\pm$ 0.070 | 0.5956 $\\pm$ 0.185 | +0.0636 | 0.5471 $\\pm$ 0.087 | 0.5907 $\\pm$ 0.209 | +0.0435 |
| | BLOSUM62 Radicalness | **0.6732 $\\pm$ 0.037** | **0.6501 $\\pm$ 0.089** | -0.0231 | **0.6380 $\\pm$ 0.069** | **0.6249 $\\pm$ 0.075** | -0.0131 |
| | Positional Shannon Entropy (-Entropy) | 0.6865 $\\pm$ 0.090 | 0.6926 $\\pm$ 0.094 | +0.0061 | 0.7251 $\\pm$ 0.089 | 0.7224 $\\pm$ 0.078 | -0.0027 |
| | Multi-Herpesvirus PSSM | 0.6932 $\\pm$ 0.093 | 0.6988 $\\pm$ 0.093 | +0.0055 | 0.7237 $\\pm$ 0.101 | 0.7271 $\\pm$ 0.085 | +0.0034 |
| | ProtBERT (Zero-Shot) | 0.7384 $\\pm$ 0.057 | 0.7193 $\\pm$ 0.046 | -0.0191 | 0.7530 $\\pm$ 0.071 | 0.7346 $\\pm$ 0.020 | -0.0184 |
| | **ESM-2 650M (Zero-Shot)** | **0.7255 $\\pm$ 0.148** | **0.7132 $\\pm$ 0.123** | **-0.0123** | **0.7397 $\\pm$ 0.157** | **0.7303 $\\pm$ 0.122** | **-0.0094** |
| **UL30-Aciclovir** | Position-Only RF (Negative Control) | 0.5347 $\\pm$ 0.069 | 0.6154 $\\pm$ 0.321 | +0.0807 | 0.5347 $\\pm$ 0.069 | 0.6241 $\\pm$ 0.345 | +0.0894 |
| | BLOSUM62 Radicalness | 0.5441 $\\pm$ 0.240 | 0.5077 $\\pm$ 0.101 | -0.0364 | 0.5494 $\\pm$ 0.250 | 0.4908 $\\pm$ 0.085 | -0.0586 |
| | Positional Shannon Entropy (-Entropy) | 0.7261 $\\pm$ 0.213 | 0.7080 $\\pm$ 0.130 | -0.0182 | 0.6999 $\\pm$ 0.260 | 0.6806 $\\pm$ 0.157 | -0.0193 |
| | Multi-Herpesvirus PSSM | 0.6340 $\\pm$ 0.176 | 0.7092 $\\pm$ 0.118 | +0.0753 | 0.6507 $\\pm$ 0.192 | 0.7189 $\\pm$ 0.115 | +0.0682 |
| | ProtBERT (Zero-Shot) | 0.6989 $\\pm$ 0.225 | 0.8147 $\\pm$ 0.101 | +0.1158 | 0.6957 $\\pm$ 0.242 | 0.8072 $\\pm$ 0.109 | +0.1114 |
| | **ESM-2 650M (Zero-Shot)** | **0.7887 $\\pm$ 0.171** | **0.7966 $\\pm$ 0.107** | **+0.0079** | **0.7875 $\\pm$ 0.177** | **0.7875 $\\pm$ 0.110** | **-0.0001** |
| **UL30-Foscarnet** | Position-Only RF (Negative Control) | 0.5537 $\\pm$ 0.074 | 0.6105 $\\pm$ 0.248 | +0.0568 | 0.5278 $\\pm$ 0.032 | 0.6001 $\\pm$ 0.288 | +0.0723 |
| | BLOSUM62 Radicalness | 0.5080 $\\pm$ 0.289 | 0.5056 $\\pm$ 0.140 | -0.0024 | 0.5074 $\\pm$ 0.300 | 0.4847 $\\pm$ 0.143 | -0.0227 |
| | Positional Shannon Entropy (-Entropy) | 0.8607 $\\pm$ 0.096 | 0.8602 $\\pm$ 0.101 | -0.0005 | 0.8246 $\\pm$ 0.082 | 0.8338 $\\pm$ 0.097 | +0.0092 |
| | Multi-Herpesvirus PSSM | 0.7580 $\\pm$ 0.144 | 0.8201 $\\pm$ 0.099 | +0.0621 | 0.7832 $\\pm$ 0.095 | 0.8334 $\\pm$ 0.099 | +0.0503 |
| | ProtBERT (Zero-Shot) | 0.7433 $\\pm$ 0.163 | 0.8675 $\\pm$ 0.102 | +0.1242 | 0.7423 $\\pm$ 0.159 | 0.8619 $\\pm$ 0.110 | +0.1197 |
| | **ESM-2 650M (Zero-Shot)** | **0.8793 $\\pm$ 0.144** | **0.8977 $\\pm$ 0.108** | **+0.0184** | **0.8809 $\\pm$ 0.125** | **0.8932 $\\pm$ 0.101** | **+0.0122** |

Three primary conclusions emerge from the 3D spatial holdout experiment:

1. **Zero-Shot PLM Invariance to Spatial Holdout:** ESM-2 650M generalizes with remarkable stability between 1D sequence block and 3D spatial neighborhood holdout ($\Delta = -0.0001$ on UL30-ACV full cohort; $+0.0079$ on resolved core; $+0.0122$ to $+0.0184$ on UL30-FOS). Furthermore, inter-fold standard deviations decrease significantly under 3D holdout (e.g., from 0.1772 to 0.1099 on UL30-ACV), because 3D partitioning eliminates the pathological single-class fold that affected linear sequence Block 0.
2. **Biological Dissociation Invariant in 3D Space:** Chemical dissimilarity alone (BLOSUM62 radicalness) consistently clears chance on the catalytic kinase (UL23-ACV: 0.6249 full, 0.6501 core) but sits squarely at chance on the allosteric polymerase (UL30-ACV: 0.4908 full, 0.5077 core; UL30-FOS: 0.4847 full, 0.5056 core). This demonstrates that the biological target dissociation is an intrinsic biophysical property that holds regardless of whether domain partitions are cut linearly along the sequence or geometrically in 3D space.
3. **Robustness of the PLM Advantage over Alignment Baselines:** Under 3D spatial holdout, ESM-2 650M maintains a consistent lead over both positional Shannon entropy and multi-herpesvirus PSSM on DNA polymerase (+0.1069 vs. entropy on ACV; +0.0594 vs. entropy and +0.0598 vs. PSSM on FOS), confirming that the PLM signal is not an artifact of 1D sequence segmentation.

---
"""

# Replace ### 3.6 Drug Conditioning with ### 3.6 3D Spatial Holdout followed by ### 3.7 Drug Conditioning
sec_3_6_target = "### 3.6 Drug Conditioning and Negative Control"
sec_3_7_replacement = sec_3_6_content + "\n### 3.7 Drug Conditioning and Negative Control"

if sec_3_6_target in text:
    text = text.replace(sec_3_6_target, sec_3_7_replacement)
    print("Added Section 3.6 and renumbered Drug Conditioning to Section 3.7.")
else:
    print("Warning: Section 3.6 target not found.")

# 3. Update Discussion on 3D spatial holdout from proposal to completed evaluation
disc_target = """Two natural extensions follow directly from these findings:
1. **3D Spatial Holdout via Cryo-EM Complexes:** Rather than partitioning sequences along 1D sequence windows, residues can be clustered by 3D Euclidean proximity using the cryo-EM holoenzyme coordinates (`8V1Q` for general active-site geometry, `8EXX` for foscarnet coordination, and `8V1T` for ACV-TP). Holding out entire 3D spatial neighborhoods directly tests the mechanistic hypothesis through its failure mode: if zero-shot PLM performance collapses under 3D structural holdout, it provides direct empirical proof that language models leverage learned spatial contact topologies to predict resistance.
2. **Cross-Herpesvirus Replication:** Mapping homologous positions between HSV-1 and orthologous targets in HerpesDRG (644 active HCMV rows across UL54/UL97 and 76 VZV rows across ORF28/ORF36)"""

disc_replacement = """The execution of our 3D spatial-neighborhood holdout (§3.6, Table 6, Figure 5) provides direct empirical evidence resolving the structural context hypothesis. Rather than collapsing under 3D spatial holdout, zero-shot ESM-2 650M maintains stable generalization ($\Delta \le 0.018$ across all targets), while inter-fold variance decreases due to balanced class distribution across globular domains. Furthermore, the radicalness dissociation is strictly preserved in 3D Euclidean space: BLOSUM62 radicalness clears chance on the catalytic kinase (UL23, AUROC $\approx 0.625–0.650$) but sits at chance on the allosteric polymerase (UL30, AUROC $\approx 0.485–0.508$). This confirms that the observed biological dissociation is an intrinsic biophysical feature of drug-target mechanisms—catalytic inactivation versus conformational remodeling—rather than an artifact of sequence segmentation.

Building on these findings, the most decisive claim upgrade available for future investigation is **Cross-Herpesvirus Replication**: mapping homologous positions between HSV-1 and orthologous targets in HerpesDRG (644 active HCMV rows across UL54/UL97 and 76 VZV rows across ORF28/ORF36)"""

if disc_target in text:
    text = text.replace(disc_target, disc_replacement)
    print("Discussion updated to reflect completed 3D spatial holdout.")
else:
    print("Warning: Discussion target not found.")

# Save to manuscript draft
with open(manuscript_path, 'w', encoding='utf-8') as f:
    f.write(text)
print(f"Saved updated manuscript draft to: {manuscript_path}")

# Save to brain artifact
with open(artifact_path, 'w', encoding='utf-8') as f:
    f.write(text)
print(f"Saved updated artifact to: {artifact_path}")

