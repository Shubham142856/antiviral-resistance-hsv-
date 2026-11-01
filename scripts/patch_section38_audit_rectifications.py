"""
patch_section38_audit_rectifications.py
Updates Section 3.8 and Discussion in manuscript/manuscript_draft.md to:
1. Provide the complete Table 7 for Supervised Dual-Buffered Domain Benchmark
2. Detail the sign resolution: b5 = +0.0368 is algebraically invariant to substitution signing ((-S)^2 = S^2)
3. Report corrected linear severity b2 = +0.3629 [+0.1235, +0.6831] on UL23
4. Frame tolerance-window null as 'not supported at this sample size (N_pos=62)'
5. Document Where vs What decomposition with exact numbers (0.8606, 0.4220, 0.8603)
6. Clarify replication scope: directionally consistent in n=1 additional pair (HCMV); VZV polymerase not evaluable
7. Synchronize the manuscript artifact
"""
import pathlib
import shutil

p = pathlib.Path('manuscript/manuscript_draft.md')
text = p.read_text(encoding='utf-8')

# 1. Update Discussion paragraph
old_disc = "Building on these findings, our pre-registered **Cross-Herpesvirus Replication** (§3.8) across orthologous targets in HerpesDRG confirms that this radicalness dissociation replicates in betaherpesvirus (HCMV UL97 kinase AUROC = 0.7005 [0.5943, 0.7960] vs. UL54 polymerase AUROC = 0.4377–0.5522) and alphaherpesvirus (VZV ORF36 kinase AUROC = 0.7611), demonstrating that target-dependent radicalness separation is indeed a family-wide empirical property across *Herpesviridae*."

new_disc = "Building on these findings, our cross-herpesvirus replication (§3.8) across orthologous targets in HerpesDRG demonstrates that the radicalness dissociation is **directionally consistent in one additional independent pair** (HCMV UL97 kinase AUROC = 0.7005 [0.5943, 0.7960] vs. UL54 polymerase AUROC = 0.4377–0.5522 across 3 drugs). On VZV, the kinase arm (ORF36 AUROC = 0.7611) exhibits high radicalness discrimination, but the polymerase arm (ORF28) had only 4 clinical resistant variants and was not evaluable. The empirical contrast is therefore demonstrated across two independent target pairs (HSV and HCMV), providing suggestive cross-species consistency rather than a proven family-wide law."

if old_disc in text:
    text = text.replace(old_disc, new_disc)
    print("Discussion updated.")
else:
    print("Warning: old_disc not found.")

# 2. Extract and replace Section 3.8
start_token = "### 3.8 Dual-Buffered Domain Benchmark, Pre-Registered Tolerance-Window Test, and Cross-Herpesvirus Replication"
end_token = "## 4. Discussion"

start_idx = text.find(start_token)
end_idx = text.find(end_token)

if start_idx != -1 and end_idx != -1:
    new_sec38 = """### 3.8 Dual-Buffered Domain Benchmark, Pre-Registered Tolerance-Window Test, and Cross-Herpesvirus Replication

To address spatial interpolation leakage at its structural root while rigorously testing the mechanistic hypotheses generated in earlier sections, we implemented an evaluation protocol encompassing dual-buffered domain holdout, a parametric tolerance-window test, exploratory variance decomposition, and cross-herpesvirus replication.

#### Dual-Buffered Leave-One-Domain-Out Protocol
Standard leave-one-domain-out holdout partitions proteins by structural domains (UL30: NTD 1–362, Exonuclease 363–593, Palm 594–766 & 826–956, Fingers 767–825, Thumb & C-term 957–1235; UL23: N-term 1–45, ATP Core 46–160, Nucleoside Pocket 161–240, Lid/Dimer 241–320, C-term 321–376 per Liu et al. 2006 and Shankar et al. 2024). To eliminate both sequence-adjacency and spatial-contact leakage across domain interfaces, we implemented a dual exclusion buffer: when holding out domain $D$, training residues within $\pm 5$ aa or $< 12.0$ \\AA in 3D space (cryo-EM holoenzyme `8V1Q` and crystal `1OF1`) are excluded from training.

Under this dual-buffered protocol (`logs/supervised_dual_buffer_benchmark.json`), the Position-Only Random Forest yields a pooled out-of-fold AUROC of **0.5185** on UL30-ACV, **0.3903** on UL30-FOS, and **0.3775** on UL23-ACV, successfully satisfying our pre-registered split-validity check ($[0.45, 0.55]$ on the primary target). Table 7 reports the complete comparative benchmark under this unified, leakage-neutral split.

**Table 7 | Supervised Predictors vs. Zero-Shot Baselines under Dual-Buffered Domain Holdout ($\pm 5$ aa, $< 12.0$ \\AA).** Pooled out-of-fold AUROCs evaluated across held-out structural domains.

| Target Pair | $N_{\\text{total}}$ ($N_{\\text{pos}}$) | Position-Only RF (Negative Control) | Zero-Shot ESM-2 650M | Zero-Shot ProtBERT 420M | Zero-Shot BLOSUM62 | Supervised Logistic Reg ($L_2$) | Supervised Random Forest |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **UL30-Aciclovir** | 332 (62) | **0.5185** | 0.8240 | 0.8398 | 0.4942 | **0.8688** | 0.8386 |
| **UL30-Foscarnet** | 311 (33) | **0.3903** | 0.8409 | 0.8545 | 0.4531 | **0.8930** | 0.8652 |
| **UL23-Aciclovir** | 286 (120) | **0.3775** | 0.7881 | 0.7243 | 0.6130 | 0.6347 | 0.6623 |

*Methodological Caveats:* Two operational limitations must be noted. First, the minimal passing buffer (12.0 \\AA) was the upper bound of our tested grid ($\\{8, 10, 12\\}$ \\AA), and the position-only RF AUROC (0.5185) sits near the upper threshold of the 0.55 neutral limit. Second, pooled out-of-fold AUROC concatenates predictions across domains whose underlying resistant base rates vary dramatically (from 0.0% in NTD to 88.9% in Fingers on UL30-ACV), which conflates within-domain ranking with across-domain calibration shift.

#### Pre-Registered Tolerance-Window Test (Negative Empirical Finding)
We pre-registered a formal test of the \"Where $\\times$ What\" tolerance-window hypothesis (`docs/WHERE_WHAT_PREREGISTRATION.md`, SHA-256 `ca4d2e7ff1...`, frozen at Git commit `ded62cb` before test execution). The hypothesis predicted that on the essential polymerase (UL30), resistance requires intermediate substitution severity at constrained positions, producing negative interaction curvature ($\\beta_5 < 0$ on $\\text{Where} \\times \\text{Severity}^2$).

Fitting $L_2$-penalized logistic models with 1,000 position-clustered bootstrap resamples yielded:
$$\\beta_5 = +0.0368, \\quad 95\\%\\text{ CI: } [-0.1868, +0.3039]$$
On the primary target (UL30-ACV), the point estimate is slightly positive and the 95% confidence interval broadly includes zero. Notably, because the interaction term is quadratic in severity ($\\text{Where} \\times \\text{Severity}^2$), $\\beta_5$ is algebraically invariant to whether severity is oriented as $\\Delta\\text{LLR}$ or $-\\Delta\\text{LLR}$ ($(-S)^2 = S^2$). Under our pre-registered decision rule, **the tolerance-window hypothesis is not supported by the clinical data at this sample size** ($N_{\\text{pos}}=62$ on UL30-ACV). We report this honestly as an unconfirmed prediction. However, because the 95% CI spans from $-0.19$ to $+0.30$, a moderate negative effect cannot be definitively excluded; at 62 resistant events, the test is statistically constrained, and masked-marginal $\\Delta\\text{LLR}$ represents only one proxy for mutation severity.

Conversely, on thymidine kinase (UL23-ACV), where resistance operates via monotonic catalytic loss-of-function, orienting severity in the pre-registered direction ($\\Delta\\text{LLR} \\ge 0$ for disruptive substitutions) yields a strictly positive linear coefficient:
$$\\beta_2 = +0.3629, \\quad 95\\%\\text{ CI: } [+0.1235, +0.6831]$$
with negligible curvature ($\\beta_3 = +0.0423$, 95% CI: $[-0.1485, +0.2463]$), confirming monotonic disruption on the kinase.

#### Exploratory Decomposition: Positional Constraint vs. Substitution Severity
To determine why zero-shot PLMs succeed on UL30 while amino acid radicalness fails, we conducted an exploratory post hoc variance decomposition separating positional constraint from within-residue substitution severity:
1. **Positional Constraint (`Where` alone):** Standardized positional constraint ($Z(-\\mu_{\\text{PLM-19}}) - Z(d_{\\text{active-site}})$) achieves an AUROC of **0.8606** on UL30-ACV and **0.9097** on UL30-FOS.
2. **Within-Residue Severity (`What` alone):** Relative mutation severity within position achieves an AUROC of **0.4220** on UL30-ACV and **0.3935** on UL30-FOS (below chance).
3. **Incremental Value:** Combining Where + Severity yields an AUROC of **0.8603** on UL30-ACV ($\\Delta = -0.0003$ over Where alone).
4. **Target Difference:** In sharp contrast, on UL23-ACV, severity alone achieves **0.5978** (and BLOSUM62 achieves **0.6130**), and adding severity to Where raises AUROC from 0.7265 to 0.7533 ($\\Delta = +0.0268$).

This establishes that **the zero-shot PLM signal on essential viral polymerase is predominantly driven by positional conservation and structural active-site localization, rather than substitution-specific chemical radicalness**. Because 98.4% of evaluated clinical UL30 positions are single-label (purely resistant or purely sensitive; only 2 mixed positions exist across the entire ACV cohort), position-level conservation provides virtually all discriminative power.

#### Cross-Herpesvirus Replication (HCMV and VZV)
To evaluate whether the empirical target-dependent radicalness dissociation generalizes beyond HSV-1, an analysis plan was specified prior to execution (`docs/CROSS_VIRUS_REPLICATION_PREREGISTRATION.md`, SHA-256 `429fd086...`, committed following the run) on the curated HerpesDRG-DB catalog (`logs/cross_virus_replication_results.json`):
- **Cytomegalovirus (HCMV) UL97 Kinase:** BLOSUM62 radicalness achieves an AUROC of **0.7005** (row bootstrap 95% CI: $[0.5943, 0.7960]$; position-cluster bootstrap 95% CI: $[0.6030, 0.7917]$ across 79 unique positions; $N=107$, 45 resistant, 62 sensitive), clearing chance.
- **HCMV UL54 DNA Polymerase:** BLOSUM62 radicalness operates near chance across all three therapeutic drugs: Ganciclovir AUROC = **0.5375** ($N=230$, 76 resistant), Foscarnet AUROC = **0.4377** ($N=160$, 50 resistant), and Cidofovir AUROC = **0.5522** ($N=149$, 57 resistant). All 95% bootstrap CIs span 0.50.
- **Varicella-Zoster Virus (VZV):** ORF36 thymidine kinase achieves **0.7611** (95% CI: $[0.5748, 0.9219]$, $N=28$, 18 resistant, 10 sensitive). However, on VZV ORF28 (DNA polymerase), only 4 clinical resistant missense variants exist in the published literature ($N_{\\text{res}}=4, N=24$, AUROC = 0.5500, 95% CI: $[0.2273, 0.9126]$); the polymerase counter-arm is therefore underpowered and **not evaluable** for contrast testing.

*Replication Scope & Limitations:* Across the two evaluable independent gene pairs (HSV-1 UL23/UL30 and HCMV UL97/UL54), the radicalness contrast matches directionally ($n=1$ replication pair). However, because UL97 loss impairs viral growth kinetics in cell culture (rendering it partially constrained rather than strictly dispensable) and clinical resistance mutations concentrate within specific kinase hotspots (~460, 520, 590–607), and VZV polymerase lacks sufficient resistant mutations to establish contrast, these results demonstrate directional consistency across two pairs rather than a proven family-wide law.

\n\n"""
    text = text[:start_idx] + new_sec38 + text[end_idx:]
    p.write_text(text, encoding='utf-8')
    print("Section 3.8 successfully replaced!")
    
    # Sync artifact
    artifact_path = pathlib.Path(r'C:\Users\HP\.gemini\antigravity\brain\7e2c0f5f-723c-45b3-9f56-f3ec4abc7b17\hsv_antiviral_resistance_research_paper.md')
    shutil.copyfile(p, artifact_path)
    print("Artifact successfully updated.")
else:
    print(f"Error: indices not found: {start_idx}, {end_idx}")
