import re

with open("manuscript/manuscript_draft.md", "r", encoding="utf-8") as f:
    text = f.read()

# 1. Abstract fix: remove uninformative zero-shot 3D invariant null claim & replace "strictly the more stringent"
abstract_old = """Evaluating models under a 3D spatial-neighborhood holdout using cryo-EM holoenzyme coordinates (PDB: `8V1Q`, chain A) and thymidine kinase crystal coordinates (PDB: `1OF1`, chain B; 1.95 Å, complexed with (South)-methanocarba-thymidine) yields a methodological null result: zero-shot PLM performance is invariant to 3D spatial neighborhood partitioning (ESM-2 650M cross-regime $|\\Delta| \\le 0.012$), indicating that PLMs do not depend on localized physical contact clusters. Crucially, evaluating the position-only negative control under 3D spatial holdout causes its AUROC to rebound from chance (0.5001) to 0.6241 on UL30-ACV because 3D geometric clusters necessarily partition sequential residues across folds, reintroducing sequence-adjacency interpolation shortcuts. This establishes that 1D contiguous-block holdout is strictly the more stringent evaluation standard, while the radicalness dissociation remains invariant in 3D Euclidean space."""

abstract_new = """Evaluating models under a 3D spatial-neighborhood holdout using cryo-EM holoenzyme coordinates (PDB: `8V1Q`, chain A) and thymidine kinase crystal coordinates (PDB: `1OF1`, chain B; 1.95 Å, complexed with (South)-methanocarba-thymidine) reveals an important methodological hazard: the position-only negative control rebounds from chance (0.5001) to 0.6241 on UL30-ACV because 3D geometric clusters necessarily partition primary sequence non-contiguously, reintroducing sequence-adjacency interpolation shortcuts across training and test folds. This demonstrates that 1D contiguous-block holdout provides a more effective control against spatial shortcut learning for sequence models than unconstrained 3D spatial clustering. Zero-shot PLMs evaluate sequence context directly without supervision, while the biological dissociation between kinase and polymerase remains invariant in 3D Euclidean space."""

assert abstract_old in text, "Abstract old text not found!"
text = text.replace(abstract_old, abstract_new)

# 2. Contribution 4: Correct ProtBERT claim
contrib4_old = """4. **Zero-Shot PLMs Generalize Across Structural Blocks:** Zero-shot protein language models (ESM-2 650M and ProtBERT) maintain high predictive performance under contiguous-block holdout (0.74 to 0.88 AUROC) where tabular models fail. On UL23, PLMs tie alignment baselines; on UL30, they retain a directional advantage over all alignment baselines tested, but this advantage is not individually statistically significant at $n=4$ evaluable blocks."""

contrib4_new = """4. **Evaluation of Zero-Shot PLMs Across Structural Blocks:** Evaluated across identical contiguous-block splits, ESM-2 (650M) satisfies our pre-registered criteria ($\\ge 0.65$ mean AUROC, 95% $t$-CI lower bound $> 0.50$) across all three target pairs (0.74 to 0.88 AUROC). ProtBERT achieves high mean AUROCs (0.696 to 0.753) and satisfies criteria on UL23-ACV, but exhibits higher block-level variance across UL30 blocks. While PLMs display a directional advantage over alignment baselines on UL30, this difference is not individually statistically significant under block-level testing ($n=4$ evaluable blocks)."""

assert contrib4_old in text, "Contrib 4 old text not found!"
text = text.replace(contrib4_old, contrib4_new)

# 3. Contribution 7: Remove PLMs invariant null claim & strictly more stringent
contrib7_old = """7. **Methodological Audit of 3D Spatial Holdout:** Evaluating 3D spatial neighborhood holdout using cryo-EM holoenzyme coordinates (`8V1Q`) reveals that 3D geometric clustering re-introduces sequence-adjacency leakage (position-only RF rebounds from chance to 0.6241), demonstrating that 1D contiguous-block holdout is strictly the more stringent evaluation standard. Furthermore, zero-shot PLM representations remain invariant between 1D and 3D holdouts ($|\\Delta| \\le 0.012$), constituting a null result regarding structural contact dependence."""

contrib7_new = """7. **Vulnerability of 3D Spatial Clustering to Sequence-Adjacency Leakage:** Evaluating a position-only negative control under 3D spatial neighborhood clustering (using cryo-EM holoenzyme coordinates PDB: `8V1Q`) reveals an unexpected structural vulnerability: 3D geometric clusters partition primary sequence non-contiguously, allowing sequence-indexed models to rebound from chance (0.5001) to 0.6241 via sequence-adjacency interpolation across folds. Consequently, 1D contiguous-block holdout provides a more effective control against spatial interpolation shortcuts for sequence models than unconstrained 3D spatial clustering."""

assert contrib7_old in text, "Contrib 7 old text not found!"
text = text.replace(contrib7_old, contrib7_new)

# 4. Table 1: Clarify Table 1 vs Table 6 0.5001 vs 0.5347 and replace "proves"
table1_note_old = """Under position-grouped cross-validation, the position-only model achieves AUROCs of **0.8715** (UL30-ACV) and **0.8684** (UL30-FOS). Because catalytic resistance mutations reside in spatial proximity, the tree interpolates labels from neighboring active-site residues. Under contiguous-block holdout, this shortcut is neutralized, and the model collapses to **0.5001** and **0.5166**. This proves that position-grouped cross-validation is vulnerable to spatial confounding."""

table1_note_new = """Under position-grouped cross-validation, the position-only model achieves AUROCs of **0.8715** (UL30-ACV) and **0.8684** (UL30-FOS). Because catalytic resistance mutations reside in spatial proximity, the tree interpolates labels from neighboring active-site residues. Under contiguous-block holdout, this shortcut is neutralized, and the model collapses to **0.5001** and **0.5166** (20-repeat CV mean $\\pm$ SD). In the single primary 5-block partition, position-only RF produces chance predictions (AUROC = 0.5000) on 3 of 4 evaluable folds, with a fold-mean of 0.5347 $\\pm$ 0.0694 and a pooled out-of-fold AUROC of 0.4591 [95% CI: 0.3741, 0.5468], confirming collapse across both repeated and single-split regimes. This demonstrates that position-grouped cross-validation is vulnerable to spatial confounding."""

assert table1_note_old in text, "Table 1 note old text not found!"
text = text.replace(table1_note_old, table1_note_new)

# 5. Table 2: Remove (Primary), (Secondary), (Baseline) and fix post-hoc primary designation
text = text.replace("| *(Primary)* | **UL30-Aciclovir** |", "| | **UL30-Aciclovir** |")
text = text.replace("| *(Secondary)* | UL30-Aciclovir |", "| | UL30-Aciclovir |")
text = text.replace("| *(Baseline)* | UL30-Aciclovir |", "| | UL30-Aciclovir |")

table2_follow_old = """ESM-2 650M is the only model that satisfies the pre-registered stopping rule on all three pairs ($\\ge 0.65$ mean and lower bound of 95% $t$-CI $> 0.50$). ProtBERT achieves high mean AUROCs (0.696 to 0.753), but exhibits higher fold-to-fold variance across UL30 blocks. ESM-2 650M was designated the primary model post hoc based on this result."""

table2_follow_new = """All three zero-shot protein language models were evaluated under identical contiguous-block splits. ESM-2 650M satisfies the pre-registered stopping rule across all three targets ($\\ge 0.65$ mean and 95% $t$-CI lower bound $> 0.50$). ProtBERT achieves high mean AUROCs (0.6957 to 0.7530), satisfying criteria on UL23-ACV while exhibiting wider confidence intervals on UL30 due to fold-level variance across the 4 evaluable blocks. ESM-2 35M achieves moderate mean AUROCs (0.6110 to 0.7481) but fails the conservative $t$-CI threshold across targets."""

assert table2_follow_old in text, "Table 2 follow old text not found!"
text = text.replace(table2_follow_old, table2_follow_new)

# 6. Section 3.6: Update Table 6 note and Conclusions 1 & 2
table6_note_old = """*Values are mean cross-validation AUROC $\\pm$ standard deviation across folds, evaluated on the Structurally Resolved Core ($N=233, 262, 242$) and Full Cohort ($N=286, 332, 311$). Note on 1D Block baseline: Table 1 reports canonical 20-repeat CV estimates (0.5001 $\\pm$ 0.0009 and 0.5166 $\\pm$ 0.0105 for UL30-ACV and UL30-FOS). A single seed-42 partition yields 0.5347 and 0.5278 due to single-split variance; under 3D spatial holdout, the position-only negative control rebounds systematically to 0.6001–0.6241 due to sequence-adjacency leakage.*"""

table6_note_new = """*Values are mean cross-validation AUROC $\\pm$ standard deviation across folds, evaluated on the Structurally Resolved Core ($N=233, 262, 242$) and Full Cohort ($N=286, 332, 311$). Note on 1D Block baseline: Under the single 5-block partition, position-only RF predicts at pure chance (0.5000) on 3 of 4 evaluable folds with a fold mean of 0.5347 $\\pm$ 0.0694 (pooled out-of-fold AUROC = 0.4591 [0.3741, 0.5468]; Table 1 reports canonical 20-repeat CV average: 0.5001 $\\pm$ 0.0009). Under 3D spatial holdout, the position-only negative control rebounds systematically to 0.6001–0.6241 due to sequence-adjacency leakage.*"""

assert table6_note_old in text, "Table 6 note old text not found!"
text = text.replace(table6_note_old, table6_note_new)

sec36_c1_c2_old = """1. **3D Spatial Clustering Re-introduces Sequence-Adjacency Leakage:** Evaluating the position-only negative control under 3D spatial clustering exposes an unexpected structural vulnerability: the position-only random forest rebounds from near chance under 1D contiguous block holdout (0.5001 on UL30-ACV) to **0.6241 $\\pm$ 0.345** under 3D spatial holdout (and from 0.5166 to **0.6001 $\\pm$ 0.288** on UL30-FOS). Because 3D geometric clusters partition residues based on spatial contact rather than sequence continuity, sequential residues that are adjacent in 1D sequence frequently belong to different 3D spatial folds. A tree model indexed on sequence position can therefore interpolate between sequence-adjacent residues across training and test folds. This demonstrates empirically that **1D contiguous-block holdout is strictly the more stringent evaluation standard**, as it eliminates sequence-adjacency interpolation shortcuts that 3D geometric clustering inadvertently re-creates.
2. **Methodological Null Result Regarding PLM Structural Neighborhood Reliance:** Pre-registered evaluation specified that if zero-shot PLMs depend primarily on localized physical 3D contact environments, holding out 3D spatial neighborhoods would induce a significant performance drop relative to 1D sequence blocks. Instead, ESM-2 650M exhibits near-zero cross-regime divergence ($|\\Delta| \\le 0.012$ across all three target pairs; $0.7875 \\to 0.7875$ on UL30-ACV; $0.8809 \\to 0.8932$ on UL30-FOS). This confirms a **methodological null result**: 3D spatial grouping provides no incremental evaluation penalty over 1D contiguous-block holdout, indicating that zero-shot PLM predictive power is driven by broad evolutionary context along the primary sequence rather than reliance on localized 3D contact clusters."""

sec36_c1_c2_new = """1. **3D Spatial Clustering Re-introduces Sequence-Adjacency Leakage:** Evaluating the position-only negative control under 3D spatial clustering exposes an unexpected structural vulnerability: the position-only random forest rebounds from near chance under 1D contiguous block holdout (0.5001 on UL30-ACV) to **0.6241 $\\pm$ 0.345** under 3D spatial holdout (and from 0.5166 to **0.6001 $\\pm$ 0.288** on UL30-FOS). Because 3D geometric clusters partition residues based on spatial contact rather than sequence continuity, sequential residues that are adjacent in 1D sequence frequently belong to different 3D spatial folds. A tree model indexed on sequence position can therefore interpolate between sequence-adjacent residues across training and test folds. This demonstrates empirically that **1D contiguous-block holdout provides a more effective control against spatial interpolation shortcuts** than unconstrained 3D spatial clustering, as it eliminates sequence-adjacency leakage that 3D geometric clustering inadvertently re-creates.
2. **Zero-Shot Score Invariance Across Partition Schemes:** Because zero-shot PLM masked-marginal scores ($-\\Delta\\text{LLR}$) are computed deterministically from sequence without supervision or fold-dependent training, variant predictions are identical regardless of partition geometry. Evaluating zero-shot scores across 1D contiguous blocks versus 3D spatial clusters yields nearly identical fold-averaged AUROCs ($|\\Delta| \\le 0.012$ for ESM-2 650M), confirming that changing the test partition geometry does not distort out-of-fold AUROC aggregations."""

assert sec36_c1_c2_old in text, "Sec 3.6 C1/C2 old text not found!"
text = text.replace(sec36_c1_c2_old, sec36_c1_c2_new)

# 7. Discussion update (line 355)
disc_old = """Our evaluation of 3D spatial-neighborhood holdout (§3.6, Table 6, Figure 5) yields two vital methodological findings. First, rather than revealing an over-reliance on local physical contact environments, zero-shot ESM-2 650M exhibits near-zero cross-regime divergence between 1D sequence blocks and 3D Euclidean spatial clusters ($|\\Delta| \\le 0.012$ across all three targets). By our pre-registered evaluation criterion, this constitutes a **methodological null result**: 3D spatial grouping provides no incremental evaluation penalty beyond 1D contiguous blocks, indicating that zero-shot PLM representations do not depend on localized 3D contact geometries to make accurate predictions. Second, evaluating our negative control reveals that 3D geometric clustering re-introduces sequence-adjacency leakage: because 3D clusters necessarily divide contiguous primary sequence stretches, a tree model indexed on sequence position rebounds from pure chance under 1D block holdout ($0.5001$) to **0.6241** under 3D spatial holdout on UL30-ACV. This establishes that **1D contiguous-block holdout is strictly the more stringent evaluation standard**, as it eliminates sequence-adjacency interpolation shortcuts that 3D spatial clustering inadvertently re-creates."""

disc_new = """Our evaluation of 3D spatial-neighborhood holdout (§3.6, Table 6, Figure 5) yields two vital methodological findings. First, because zero-shot PLM scores are computed directly from sequence context without fold training, evaluating them across 1D sequence blocks versus 3D Euclidean spatial clusters yields essentially invariant aggregations ($|\\Delta| \\le 0.012$), reflecting the static nature of zero-shot representations. Second, evaluating the position-only negative control reveals that 3D geometric clustering re-introduces sequence-adjacency leakage: because 3D clusters divide primary sequence into non-contiguous segments, a tree model indexed on sequence position rebounds from chance under 1D block holdout (0.5001) to **0.6241** under 3D spatial holdout on UL30-ACV. This establishes that **1D contiguous-block holdout provides a more effective control against sequence-adjacency interpolation shortcuts** than 3D spatial clustering."""

assert disc_old in text, "Discussion old text not found!"
text = text.replace(disc_old, disc_new)

with open("manuscript/manuscript_draft.md", "w", encoding="utf-8") as f:
    f.write(text)

print("Successfully patched manuscript/manuscript_draft.md with all Phase 1 corrections!")

