import pathlib

def main():
    p = pathlib.Path(__file__).resolve().parent.parent / "manuscript/manuscript_draft.md"
    txt = p.read_text(encoding="utf-8")

    # 1. Scope fix in abstract
    old_scope = "in Herpes Simplex Virus (HSV-1 and HSV-2) resistance benchmarks"
    new_scope = "in Herpes Simplex Virus type 1 (HSV-1) resistance benchmarks"
    assert old_scope in txt, f"Missing: {old_scope}"
    txt = txt.replace(old_scope, new_scope)

    # 2. 1OF1 ligand & 3D holdout null / stringency fix in abstract
    old_abs_3d = (
        "Evaluating models under a 3D spatial-neighborhood holdout using cryo-EM holoenzyme coordinates "
        "(PDB: `8V1Q`) and thymidine kinase crystal coordinates (`1OF1`) demonstrates that zero-shot PLM representations "
        "are invariant to 3D structural domain shift (ESM-2 650M achieving mean AUROCs of 0.7303, 0.7875, and 0.8932; "
        "cross-regime $|\\Delta| \\le 0.012$), while the radicalness dissociation remains strictly invariant in 3D Euclidean space."
    )
    new_abs_3d = (
        "Evaluating models under a 3D spatial-neighborhood holdout using cryo-EM holoenzyme coordinates "
        "(PDB: `8V1Q`, chain A) and thymidine kinase crystal coordinates (PDB: `1OF1`, chain B; 1.95 Å, complexed with "
        "(South)-methanocarba-thymidine) yields a methodological null result: zero-shot PLM performance is invariant to "
        "3D spatial neighborhood partitioning (ESM-2 650M cross-regime $|\\Delta| \\le 0.012$), indicating that PLMs do "
        "not depend on localized physical contact clusters. Crucially, evaluating the position-only negative control under "
        "3D spatial holdout causes its AUROC to rebound from chance (0.5001) to 0.6241 on UL30-ACV because 3D geometric clusters "
        "necessarily partition sequential residues across folds, reintroducing sequence-adjacency interpolation shortcuts. "
        "This establishes that 1D contiguous-block holdout is strictly the more stringent evaluation standard, while the "
        "radicalness dissociation remains invariant in 3D Euclidean space."
    )
    assert old_abs_3d in txt, "Missing old_abs_3d"
    txt = txt.replace(old_abs_3d, new_abs_3d)

    # 3. Abstract radicalness divergence softening
    old_abs_rad = (
        "By benchmarking against multi-herpesvirus evolutionary baselines and BLOSUM62 amino acid radicalness, "
        "we uncover a statistically robust biological dissociation between the two resistance targets. "
        "**BLOSUM62 radicalness alone separates UL23 from UL30**: on thymidine kinase, radicalness achieves **0.6380** "
        "(95% CI [0.5524, 0.7237], every block above 0.55), while on DNA polymerase it collapses to chance "
        "(**0.5074 to 0.5494**, CIs spanning 0.50). This divergence reflects the distinct resistance mechanisms of the "
        "two enzymes: UL23 resistance arises primarily through disruptive substitutions that abolish prodrug activation, "
        "whereas UL30 resistance requires subtle conformational changes that preserve replication while discriminating against inhibitors."
    )
    new_abs_rad = (
        "By benchmarking against multi-herpesvirus evolutionary baselines and BLOSUM62 amino acid radicalness, "
        "we observe a clear biological dissociation between the two resistance targets. "
        "**BLOSUM62 radicalness alone separates UL23 from UL30**: on thymidine kinase, radicalness achieves **0.6380** "
        "(95% CI [0.5524, 0.7237], clearing chance across all blocks), while on DNA polymerase it collapses to chance "
        "(**0.5074 to 0.5494**, CIs spanning 0.50). This divergence reflects the distinct resistance mechanisms of the "
        "two enzymes: UL23 resistance arises primarily through disruptive substitutions that abolish prodrug activation "
        "in a non-essential enzyme, whereas UL30 resistance requires subtle conformational changes that preserve essential "
        "replication while discriminating against inhibitors. However, due to limited degrees of freedom across structural domains "
        "($n=5$ vs. $n=4$ blocks), this divergence is descriptive rather than statistically decisive under block-level testing "
        "($p=0.431$ and $p=0.406$), though secondary position-cluster resampling strictly excludes zero."
    )
    assert old_abs_rad in txt, "Missing old_abs_rad"
    txt = txt.replace(old_abs_rad, new_abs_rad)

    # 4. Contribution 3 and 7
    old_contrib = (
        "3. **Radicalness Dissociation between Targets (Tested via Interaction Bootstrap):** "
        "BLOSUM62 amino acid radicalness predicts resistance on UL23 (**0.6380**, 95% CI [0.5524, 0.7237]) but collapses to chance on "
        "UL30 (**0.5074 to 0.5494**). A formal interaction test (2,000 position-cluster bootstrap resamples) confirms that "
        "radicalness predictive power is significantly higher on UL23 than UL30 ($\\Delta = +0.1187$, 95% CI: $[+0.0387, +0.2021]$, "
        "empirical $p=0.0020$ vs. UL30-ACV; $\\Delta = +0.1598$, 95% CI: $[+0.0655, +0.2597]$, $p < 0.0001$ vs. UL30-FOS; "
        "`logs/radicalness_interaction_test.json`), though block-level variance remains wide at $n=4$ blocks."
    )
    new_contrib = (
        "3. **Radicalness Dissociation between Targets:** BLOSUM62 amino acid radicalness predicts resistance on UL23 "
        "(**0.6380**, 95% CI [0.5524, 0.7237], clearing chance across all blocks) but collapses to chance on UL30 "
        "(**0.5074 to 0.5494**). Under primary block-level testing, this divergence is descriptive rather than statistically decisive "
        "($p=0.431$ and $p=0.406$ across $n=5$ vs. $n=4$ blocks), while secondary position-cluster bootstrap sensitivity analysis "
        "strictly excludes zero ($\\Delta = +0.1187$, 95% CI: $[+0.0387, +0.2021]$ vs. UL30-ACV; $\\Delta = +0.1598$, 95% CI: $[+0.0655, +0.2597]$ vs. UL30-FOS). "
        "This divergence reflects the distinct biophysical demands of catalytic disruption in non-essential TK versus subtle conformational remodeling in essential polymerase."
    )
    assert old_contrib in txt, "Missing old_contrib"
    txt = txt.replace(old_contrib, new_contrib)

    old_contrib_end = (
        "6. **Quantification of Clinical Scope:** We document that truncating and frameshift variants represent 44.4% of cataloged UL23 "
        "resistance, establishing that our benchmark rigorously evaluates the missense class where predictive modeling is genuinely needed."
    )
    new_contrib_end = (
        "6. **Quantification of Clinical Scope:** We document that truncating and frameshift variants represent 44.4% of cataloged UL23 "
        "resistance, establishing that our benchmark rigorously evaluates the missense class where predictive modeling is genuinely needed.\n"
        "7. **Methodological Audit of 3D Spatial Holdout:** Evaluating 3D spatial neighborhood holdout using cryo-EM holoenzyme coordinates (`8V1Q`) "
        "reveals that 3D geometric clustering re-introduces sequence-adjacency leakage (position-only RF rebounds from chance to 0.6241), "
        "demonstrating that 1D contiguous-block holdout is strictly the more stringent evaluation standard. Furthermore, zero-shot PLM representations "
        "remain invariant between 1D and 3D holdouts ($|\\Delta| \\le 0.012$), constituting a null result regarding structural contact dependence."
    )
    assert old_contrib_end in txt, "Missing old_contrib_end"
    txt = txt.replace(old_contrib_end, new_contrib_end)

    # 5. Methods 2.1: Döhne -> Dähne & HerpesDRG citation
    old_m21 = "including the CMI 2025 systemic review by Döhne et al., Sauerbrei et al. JAC 2016, and Burrel et al. Antiviral Res 2010"
    new_m21 = "including the CMI 2025 systemic review by Dähne et al. (PMID: 40349973), the HerpesDRG database (Charles et al., 2024; PMID: 39192205), Sauerbrei et al. JAC 2016, and Burrel et al. Antiviral Res 2010"
    assert old_m21 in txt, "Missing old_m21"
    txt = txt.replace(old_m21, new_m21)

    # 6. Results 3.3 lead-in: block level primary
    old_r33 = (
        "### 3.3 Amino Acid Radicalness Separates the Two Resistance Targets\n\n"
        "The statistically robust component of the biological dissociation between thymidine kinase and DNA polymerase is amino acid radicalness. "
        "BLOSUM62 radicalness — a single scalar per substitution requiring no alignment, no model training, and no position-level information — "
        "separates the two targets cleanly (Table 2, BLOSUM62 rows):"
    )
    new_r33 = (
        "### 3.3 Amino Acid Radicalness Separates the Two Resistance Targets\n\n"
        "The primary biological divergence observed in our benchmark between thymidine kinase and DNA polymerase is amino acid radicalness. "
        "BLOSUM62 radicalness — a single scalar per substitution requiring no alignment, no model training, and no position-level information — "
        "separates the two targets cleanly under domain holdout (Table 2, BLOSUM62 rows):"
    )
    assert old_r33 in txt, "Missing old_r33"
    txt = txt.replace(old_r33, new_r33)

    old_r33_test = (
        "To test this divergence directly and avoid the difference-of-significance fallacy (Gelman & Stern, 2006), "
        "we performed a position-cluster bootstrap interaction test (2,000 resamples over unique sequence positions; "
        "`logs/radicalness_interaction_test.json`). The predictive advantage of radicalness on UL23 over UL30 is positive "
        "and strictly excludes zero across the full cohort:\n"
        "- **UL23-ACV vs. UL30-ACV:** Observed point delta $\\Delta = +0.1187$, 95% position-cluster bootstrap CI: **[+0.0387, +0.2021]** (empirical $p=0.0020$).\n"
        "- **UL23-ACV vs. UL30-FOS:** Observed point delta $\\Delta = +0.1598$, 95% position-cluster bootstrap CI: **[+0.0655, +0.2597]** (empirical $p < 0.0001$).\n"
        "- **UL23-ACV vs. UL30-Combined:** Observed point delta $\\Delta = +0.1328$, 95% position-cluster bootstrap CI: **[+0.0616, +0.2150]** (empirical $p < 0.0001$).\n\n"
        "At the block-mean level, however, fold-level variance is wide due to limited degrees of freedom ($n=5$ vs. $n=4$; "
        "Welch's $t=0.885, p=0.431$ for ACV; $t=0.949, p=0.406$ for FOS). Thus, the dissociation is statistically robust under "
        "position-cluster resampling across the cohort, but descriptive when stratified across five structural blocks.\n\n"
        "This divergence is the paper's most interpretable finding and is independent of any PLM or alignment method."
    )
    new_r33_test = (
        "When tested formally at the primary block-mean level, wide fold-level variance across the 5 structural blocks "
        "leaves the divergence non-significant (Welch's two-sample $t$-test: $\\Delta = +0.0887, t(5.4)=0.885, p=0.431$ for ACV; "
        "$\\Delta = +0.1307, t(3.8)=0.949, p=0.406$ for FOS; both 95% CIs span zero). The block-level divergence between targets "
        "is therefore **descriptive and qualitative rather than statistically decisive**, reflecting the inherent degrees-of-freedom "
        "constraint ($n=5$ vs. $n=4$) of structural domain holdout.\n\n"
        "To test whether the divergence holds under pooled sample resampling and avoid the difference-of-significance fallacy "
        "(Gelman & Stern, 2006), we also performed a secondary position-cluster bootstrap interaction test (2,000 resamples over "
        "unique sequence positions; `logs/radicalness_interaction_test.json`). Across the pooled cohort, the advantage of radicalness "
        "on UL23 over UL30 strictly excludes zero:\n"
        "- **UL23-ACV vs. UL30-ACV:** Observed point delta $\\Delta = +0.1187$, 95% position-cluster bootstrap CI: **[+0.0387, +0.2021]** (empirical $p=0.0020$).\n"
        "- **UL23-ACV vs. UL30-FOS:** Observed point delta $\\Delta = +0.1598$, 95% position-cluster bootstrap CI: **[+0.0655, +0.2597]** (empirical $p < 0.0001$).\n"
        "- **UL23-ACV vs. UL30-Combined:** Observed point delta $\\Delta = +0.1328$, 95% position-cluster bootstrap CI: **[+0.0616, +0.2150]**, 99% CI: **[+0.0308, +0.2394]** (empirical $p < 0.0001$).\n\n"
        "However, because the central thesis of this benchmark is that position-level resampling ignores spatial confounding, "
        "we treat the block-level test as the primary standard and report the position-cluster result strictly as secondary sensitivity evidence."
    )
    assert old_r33_test in txt, "Missing old_r33_test"
    txt = txt.replace(old_r33_test, new_r33_test)

    # 7. Section 3.6: 1OF1 ligand & description
    old_1of1 = "For thymidine kinase (UL23), coordinates were extracted from the 1.9 Å aciclovir-complexed crystal structure (PDB: `1OF1`, chain B)."
    new_1of1 = "For thymidine kinase (UL23), coordinates were extracted from the 1.95 Å crystal structure of HSV-1 TK complexed with the thymidine analogue (South)-methanocarba-thymidine (ligand code: SCT; PDB: `1OF1`, chain B; Bennett et al., *FEBS Lett.* 1999)."
    assert old_1of1 in txt, "Missing old_1of1"
    txt = txt.replace(old_1of1, new_1of1)

    # 8. Table 6: Footnote reconciling Table 1 and Table 6 Position-Only RF
    old_t6_lead = (
        "**Table 6 | Dual-Regime Generalization: 1D Contiguous-Block vs. 3D Spatial-Neighborhood Holdout**\n"
        "*Values are mean cross-validation AUROC $\\pm$ standard deviation across folds, evaluated on the Structurally Resolved Core ($N=233, 262, 242$) and Full Cohort ($N=286, 332, 311$).*"
    )
    new_t6_lead = (
        "**Table 6 | Dual-Regime Generalization: 1D Contiguous-Block vs. 3D Spatial-Neighborhood Holdout**\n"
        "*Values are mean cross-validation AUROC $\\pm$ standard deviation across folds, evaluated on the Structurally Resolved Core ($N=233, 262, 242$) and Full Cohort ($N=286, 332, 311$). "
        "Note on 1D Block baseline: Table 1 reports canonical 20-repeat CV estimates (0.5001 $\\pm$ 0.0009 and 0.5166 $\\pm$ 0.0105 for UL30-ACV and UL30-FOS). "
        "A single seed-42 partition yields 0.5347 and 0.5278 due to single-split variance; under 3D spatial holdout, the position-only negative control rebounds "
        "systematically to 0.6001–0.6241 due to sequence-adjacency leakage.*"
    )
    assert old_t6_lead in txt, "Missing old_t6_lead"
    txt = txt.replace(old_t6_lead, new_t6_lead)

    # 9. Conclusions of 3.6
    old_36_concl = (
        "Three primary conclusions emerge from the 3D spatial holdout experiment:\n\n"
        "1. **Zero-Shot PLM Invariance to Spatial Holdout:** ESM-2 650M generalizes with remarkable stability between 1D sequence block and 3D spatial neighborhood holdout "
        "($\\Delta = -0.0001$ on UL30-ACV full cohort; $+0.0079$ on resolved core; $+0.0122$ to $+0.0184$ on UL30-FOS). Furthermore, inter-fold standard deviations "
        "decrease significantly under 3D holdout (e.g., from 0.1772 to 0.1099 on UL30-ACV), because 3D partitioning eliminates the pathological single-class fold "
        "that affected linear sequence Block 0.\n"
        "2. **Biological Dissociation Invariant in 3D Space:** Chemical dissimilarity alone (BLOSUM62 radicalness) consistently clears chance on the catalytic kinase "
        "(UL23-ACV: 0.6249 full, 0.6501 core) but sits squarely at chance on the allosteric polymerase (UL30-ACV: 0.4908 full, 0.5077 core; UL30-FOS: 0.4847 full, "
        "0.5056 core). This demonstrates that the biological target dissociation is an intrinsic biophysical property that holds regardless of whether domain partitions "
        "are cut linearly along the sequence or geometrically in 3D space.\n"
        "3. **Robustness of the PLM Advantage over Alignment Baselines:** Under 3D spatial holdout, ESM-2 650M maintains a consistent lead over both positional "
        "Shannon entropy and multi-herpesvirus PSSM on DNA polymerase (+0.1069 vs. entropy on ACV; +0.0594 vs. entropy and +0.0598 vs. PSSM on FOS), "
        "confirming that the PLM signal is not an artifact of 1D sequence segmentation."
    )
    new_36_concl = (
        "Three critical conclusions emerge from the 3D spatial holdout experiment:\n\n"
        "1. **3D Spatial Clustering Re-introduces Sequence-Adjacency Leakage:** Evaluating the position-only negative control under 3D spatial clustering exposes "
        "an unexpected structural vulnerability: the position-only random forest rebounds from near chance under 1D contiguous block holdout (0.5001 on UL30-ACV) "
        "to **0.6241 $\\pm$ 0.345** under 3D spatial holdout (and from 0.5166 to **0.6001 $\\pm$ 0.288** on UL30-FOS). Because 3D geometric clusters partition residues "
        "based on spatial contact rather than sequence continuity, sequential residues that are adjacent in 1D sequence frequently belong to different 3D spatial folds. "
        "A tree model indexed on sequence position can therefore interpolate between sequence-adjacent residues across training and test folds. This demonstrates empirically "
        "that **1D contiguous-block holdout is strictly the more stringent evaluation standard**, as it eliminates sequence-adjacency interpolation shortcuts that 3D geometric "
        "clustering inadvertently re-creates.\n"
        "2. **Methodological Null Result Regarding PLM Structural Neighborhood Reliance:** Pre-registered evaluation specified that if zero-shot PLMs depend primarily on "
        "localized physical 3D contact environments, holding out 3D spatial neighborhoods would induce a significant performance drop relative to 1D sequence blocks. "
        "Instead, ESM-2 650M exhibits near-zero cross-regime divergence ($|\\Delta| \\le 0.012$ across all three target pairs; $0.7875 \\to 0.7875$ on UL30-ACV; "
        "$0.8809 \\to 0.8932$ on UL30-FOS). This confirms a **methodological null result**: 3D spatial grouping provides no incremental evaluation penalty over "
        "1D contiguous-block holdout, indicating that zero-shot PLM predictive power is driven by broad evolutionary context along the primary sequence rather than "
        "reliance on localized 3D contact clusters.\n"
        "3. **Biological Dissociation Invariant in 3D Space:** Chemical dissimilarity alone (BLOSUM62 radicalness) consistently clears chance on the catalytic kinase "
        "(UL23-ACV: 0.6249 full, 0.6501 core) but sits squarely at chance on the allosteric polymerase (UL30-ACV: 0.4908 full, 0.5077 core; UL30-FOS: 0.4847 full, "
        "0.5056 core). This demonstrates that the biological target dissociation is an intrinsic biophysical feature—catalytic inactivation versus conformational "
        "remodeling—that remains invariant whether domains are segmented linearly along the sequence or geometrically in 3D space."
    )
    assert old_36_concl in txt, "Missing old_36_concl"
    txt = txt.replace(old_36_concl, new_36_concl)

    # 10. Discussion 3D holdout paragraph
    old_disc_3d = (
        "The execution of our 3D spatial-neighborhood holdout (§3.6, Table 6, Figure 5) provides direct empirical evidence resolving the structural context hypothesis. "
        "Rather than collapsing under 3D spatial holdout, zero-shot ESM-2 650M maintains stable generalization ($\\Delta \\le 0.018$ across all targets), while "
        "inter-fold variance decreases due to balanced class distribution across globular domains. Furthermore, the radicalness dissociation is strictly preserved "
        "in 3D Euclidean space: BLOSUM62 radicalness clears chance on the catalytic kinase (UL23, AUROC $\\approx 0.625–0.650$) but sits at chance on the "
        "allosteric polymerase (UL30, AUROC $\\approx 0.485–0.508$). This confirms that the observed biological dissociation is an intrinsic biophysical feature "
        "of drug-target mechanisms—catalytic inactivation versus conformational remodeling—rather than an artifact of sequence segmentation."
    )
    new_disc_3d = (
        "Our evaluation of 3D spatial-neighborhood holdout (§3.6, Table 6, Figure 5) yields two vital methodological findings. "
        "First, rather than revealing an over-reliance on local physical contact environments, zero-shot ESM-2 650M exhibits near-zero cross-regime divergence "
        "between 1D sequence blocks and 3D Euclidean spatial clusters ($|\\Delta| \\le 0.012$ across all three targets). By our pre-registered evaluation criterion, "
        "this constitutes a **methodological null result**: 3D spatial grouping provides no incremental evaluation penalty beyond 1D contiguous blocks, indicating that "
        "zero-shot PLM representations do not depend on localized 3D contact geometries to make accurate predictions. "
        "Second, evaluating our negative control reveals that 3D geometric clustering re-introduces sequence-adjacency leakage: because 3D clusters necessarily divide "
        "contiguous primary sequence stretches, a tree model indexed on sequence position rebounds from pure chance under 1D block holdout ($0.5001$) to **0.6241** "
        "under 3D spatial holdout on UL30-ACV. This establishes that **1D contiguous-block holdout is strictly the more stringent evaluation standard**, as it eliminates "
        "sequence-adjacency interpolation shortcuts that 3D spatial clustering inadvertently re-creates. Crucially, the biological dissociation remains invariant across "
        "both spatial regimes: BLOSUM62 radicalness clears chance on the catalytic kinase (UL23, AUROC $\\approx 0.625–0.650$) while remaining strictly at chance on the "
        "allosteric polymerase (UL30, AUROC $\\approx 0.485–0.508$), confirming that this divergence is an intrinsic biophysical property of the two targets rather "
        "than an artifact of sequence segmentation."
    )
    assert old_disc_3d in txt, "Missing old_disc_3d"
    txt = txt.replace(old_disc_3d, new_disc_3d)

    # 11. Discussion radicalness paragraph lead-in
    old_disc_rad = (
        "The statistically robust finding in our biological analysis is the radicalness divergence: BLOSUM62 amino acid radicalness — a zero-parameter, "
        "alignment-free feature — predicts resistance on UL23 (0.6380, 95% CI [0.5524, 0.7237]) but operates at chance on UL30 (0.5074–0.5494). "
        "This dissociation reflects the distinct resistance mechanisms of the two enzymes and holds under contiguous-block holdout with all individual UL23 blocks above 0.55."
    )
    new_disc_rad = (
        "The primary biological divergence observed in our benchmark is the radicalness dissociation: BLOSUM62 amino acid radicalness — a zero-parameter, "
        "alignment-free feature — predicts resistance on UL23 (0.6380, 95% CI [0.5524, 0.7237], clearing chance across all blocks) but operates at chance on "
        "UL30 (0.5074–0.5494). When tested formally at the block-mean level, wide fold-level variance leaves this difference descriptive and non-significant "
        "($p=0.431$ on ACV, $p=0.406$ on FOS), though secondary position-cluster resampling across the full cohort strictly excludes zero ($p=0.0020$ on ACV, $p < 0.0001$ on FOS)."
    )
    assert old_disc_rad in txt, "Missing old_disc_rad"
    txt = txt.replace(old_disc_rad, new_disc_rad)

    # 12. References additions
    old_refs = (
        "## References\n\n"
        "1. **Liu et al. (PRIME):** *PRIME: An evaluation framework for protein representation inference and generalization in viral mutation space.* *BMC Genomics* 2026; DOI: [10.1186/s12864-026-12976-5](https://doi.org/10.1186/s12864-026-12976-5).\n"
        "2. **Shankar et al. (HSV-1 Polymerase Cryo-EM):** *Viral DNA polymerase structures reveal mechanisms of antiviral drug resistance.* *Cell* 2024; 187(20):5572–5586.e15. PMID: [39197451](https://pubmed.ncbi.nlm.nih.gov/39197451/). DOI: [10.1016/j.cell.2024.07.048](https://doi.org/10.1016/j.cell.2024.07.048). PDB entries: `8V1Q`, `8V1R`, `8V1S`, `8EXX`, `8V1T`.\n"
        "3. **Safrin et al. (Foundational PRA Thresholds):** *Antiviral susceptibility of herpes simplex virus isolates from immunocompromised patients: analysis of 115 isolates from HIV-infected patients.* *Antimicrob Agents Chemother* 1994; 38(6):1246–1250. PMID: [8092821](https://pubmed.ncbi.nlm.nih.gov/8092821/).\n"
        "4. **Burrel et al. (Dye Uptake Standardization):** *Genotypic and phenotypic analysis of herpes simplex virus 1 clinical isolates resistant to acyclovir.* *Antimicrob Agents Chemother* 2010; 54(11):4833–4842. PMID: [20733037](https://pubmed.ncbi.nlm.nih.gov/20733037/).\n"
        "5. **Sauerbrei et al. (European Diagnostic Compilation):** *Characterization of herpes simplex virus isolates with ambiguous susceptibility to acyclovir and/or foscarnet.* *J Antimicrob Chemother* 2016; 71(1):6–16. PMID: [26433780](https://pubmed.ncbi.nlm.nih.gov/26433780/).\n"
        "6. **Sawhney et al. (Viral PLM Representation):** *Fine-tuning protein language models unlocks the potential of underrepresented viral proteomes.* *PeerJ* 2025; 13:e19919. PMID: [40985030](https://pubmed.ncbi.nlm.nih.gov/40985030/). DOI: [10.7717/peerj.19919](https://doi.org/10.7717/peerj.19919).\n"
        "7. **Gelman & Stern (Difference of Significance):** *The difference between \"significant\" and \"not significant\" is not itself statistically significant.* *The American Statistician* 2006; 60(4):328–331. DOI: [10.1198/000313006X152607](https://doi.org/10.1198/000313006X152607)."
    )
    new_refs = (
        "## References\n\n"
        "1. **Liu et al. (PRIME):** *PRIME: An evaluation framework for protein representation inference and generalization in viral mutation space.* *BMC Genomics* 2026; DOI: [10.1186/s12864-026-12976-5](https://doi.org/10.1186/s12864-026-12976-5).\n"
        "2. **Shankar et al. (HSV-1 Polymerase Cryo-EM):** *Viral DNA polymerase structures reveal mechanisms of antiviral drug resistance.* *Cell* 2024; 187(20):5572–5586.e15. PMID: [39197451](https://pubmed.ncbi.nlm.nih.gov/39197451/). DOI: [10.1016/j.cell.2024.07.048](https://doi.org/10.1016/j.cell.2024.07.048). PDB entries: `8V1Q`, `8V1R`, `8V1S`, `8EXX`, `8V1T`.\n"
        "3. **Charles et al. (HerpesDRG Database):** *HerpesDRG: a comprehensive resource for human herpesvirus antiviral drug resistance genotyping.* *BMC Bioinformatics* 2024; 25:283. PMID: [39192205](https://pubmed.ncbi.nlm.nih.gov/39192205/). DOI: [10.1186/s12859-024-05885-5](https://doi.org/10.1186/s12859-024-05885-5).\n"
        "4. **Dähne et al. (CMI Diagnostic Consensus Review):** *Herpes simplex virus and drug resistance—comprehensive update on resistance mutations and implications for clinical management: a narrative review.* *Clin Microbiol Infect* 2025; 31(9):1484–1490. PMID: [40349973](https://pubmed.ncbi.nlm.nih.gov/40349973/). DOI: [10.1016/j.cmi.2025.04.046](https://doi.org/10.1016/j.cmi.2025.04.046). Corrigendum: *Clin Microbiol Infect* 2026; DOI: [10.1016/j.cmi.2026.01.010](https://doi.org/10.1016/j.cmi.2026.01.010).\n"
        "5. **Bennett et al. (HSV-1 Thymidine Kinase Structure):** *Biochemical and structural characterization of (South)-methanocarbathymidine that specifically inhibits growth of herpes simplex virus.* *FEBS Lett* 1999; 443(2):121–125. PMID: [10025941](https://pubmed.ncbi.nlm.nih.gov/10025941/). DOI: [10.1016/S0014-5793(98)01701-1](https://doi.org/10.1016/S0014-5793(98)01701-1). PDB: `1OF1`.\n"
        "6. **Safrin et al. (Foundational PRA Thresholds):** *Antiviral susceptibility of herpes simplex virus isolates from immunocompromised patients: analysis of 115 isolates from HIV-infected patients.* *Antimicrob Agents Chemother* 1994; 38(6):1246–1250. PMID: [8092821](https://pubmed.ncbi.nlm.nih.gov/8092821/).\n"
        "7. **Burrel et al. (Dye Uptake Standardization):** *Genotypic and phenotypic analysis of herpes simplex virus 1 clinical isolates resistant to acyclovir.* *Antimicrob Agents Chemother* 2010; 54(11):4833–4842. PMID: [20733037](https://pubmed.ncbi.nlm.nih.gov/20733037/).\n"
        "8. **Sauerbrei et al. (European Diagnostic Compilation):** *Characterization of herpes simplex virus isolates with ambiguous susceptibility to acyclovir and/or foscarnet.* *J Antimicrob Chemother* 2016; 71(1):6–16. PMID: [26433780](https://pubmed.ncbi.nlm.nih.gov/26433780/).\n"
        "9. **Sawhney et al. (Viral PLM Representation):** *Fine-tuning protein language models unlocks the potential of underrepresented viral proteomes.* *PeerJ* 2025; 13:e19919. PMID: [40985030](https://pubmed.ncbi.nlm.nih.gov/40985030/). DOI: [10.7717/peerj.19919](https://doi.org/10.7717/peerj.19919).\n"
        "10. **Gelman & Stern (Difference of Significance):** *The difference between \"significant\" and \"not significant\" is not itself statistically significant.* *The American Statistician* 2006; 60(4):328–331. DOI: [10.1198/000313006X152607](https://doi.org/10.1198/000313006X152607)."
    )
    assert old_refs in txt, "Missing old_refs"
    txt = txt.replace(old_refs, new_refs)

    p.write_text(txt, encoding="utf-8")
    print("Successfully patched manuscript/manuscript_draft.md!")

    # Also mirror to conversation artifact
    artifact_path = pathlib.Path(r"C:\Users\HP\.gemini\antigravity\brain\7e2c0f5f-723c-45b3-9f56-f3ec4abc7b17\hsv_antiviral_resistance_research_paper.md")
    if artifact_path.parent.exists():
        artifact_path.write_text(txt, encoding="utf-8")
        print("Successfully mirrored to artifact directory!")

if __name__ == "__main__":
    main()
