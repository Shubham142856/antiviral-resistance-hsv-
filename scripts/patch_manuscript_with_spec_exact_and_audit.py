"""
patch_manuscript_with_spec_exact_and_audit.py
Updates Section 3.8 of manuscript/manuscript_draft.md with:
1. Side-by-side reporting of both tolerance runs (initial implementation vs. spec-exact rerun)
2. Position multiplicity and label purity data (86.9% single-variant, 2 mixed positions)
3. Both per-domain and pooled AUROCs for dual-buffer LODO
4. HCMV UL97 vs UL54 bootstrap contrast test (Delta = +0.1630, p = 0.0070)
5. Synchronize the research paper artifact
"""
import pathlib
import shutil

p = pathlib.Path('manuscript/manuscript_draft.md')
text = p.read_text(encoding='utf-8')

start_token = "### 3.8 Dual-Buffered Domain Benchmark, Pre-Registered Tolerance-Window Test, and Cross-Herpesvirus Replication"
end_token = "## 4. Discussion"

start_idx = text.find(start_token)
end_idx = text.find(end_token)

if start_idx != -1 and end_idx != -1:
    new_sec38 = """### 3.8 Dual-Buffered Domain Benchmark, Pre-Registered Tolerance-Window Test, and Cross-Herpesvirus Replication

To address spatial interpolation leakage at its structural root while rigorously testing the mechanistic hypotheses generated in earlier sections, we implemented an evaluation protocol encompassing dual-buffered domain holdout, a parametric tolerance-window test, exploratory variance decomposition, and cross-herpesvirus replication.

#### Dual-Buffered Leave-One-Domain-Out Protocol
Standard leave-one-domain-out holdout partitions proteins by structural domains (UL30: NTD 1–362, Exonuclease 363–593, Palm 594–766 & 826–956, Fingers 767–825, Thumb & C-term 957–1235; UL23: N-term 1–45, ATP Core 46–160, Nucleoside Pocket 161–240, Lid/Dimer 241–320, C-term 321–376 per Liu et al. 2006 and Shankar et al. 2024). To eliminate both sequence-adjacency and spatial-contact leakage across domain interfaces, we implemented a dual exclusion buffer: when holding out domain $D$, training residues within $\pm 5$ aa or $< 12.0$ \\AA in 3D space (cryo-EM holoenzyme `8V1Q` and crystal `1OF1`) are excluded from training.

Under this dual-buffered protocol (`logs/supervised_dual_buffer_benchmark.json`), the Position-Only Random Forest yields a pooled out-of-fold AUROC of **0.5185** on the primary target (UL30-ACV), successfully satisfying our pre-registered split-validity check ($[0.45, 0.55]$). Across evaluable individual domains, Position-Only RF operates consistently near chance (UL30-ACV domain mean: $0.5645$; UL30-FOS domain mean: $0.5147$; UL23-ACV domain mean: $0.4372$; Table 7). On UL30-FOS and UL23-ACV, pooled out-of-fold AUROCs drop to $0.3903$ and $0.3775$ because pooled metrics aggregate across domains with sharply divergent resistant base rates (e.g., UL30 NTD 0.0% vs. Fingers 88.9% resistant), conflating within-domain ranking with across-domain calibration shift.

**Table 7 | Supervised Predictors vs. Zero-Shot Baselines under Dual-Buffered Domain Holdout ($\pm 5$ aa, $< 12.0$ \\AA).** Performance reported as pooled out-of-fold AUROC [with evaluable per-domain mean in brackets].

| Target Pair | $N_{\\text{total}}$ ($N_{\\text{pos}}$) | Position-Only RF (Negative Control) | Zero-Shot ESM-2 650M | Zero-Shot ProtBERT 420M | Zero-Shot BLOSUM62 | Supervised Logistic Reg ($L_2$) | Supervised Random Forest |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **UL30-Aciclovir** | 332 (62) | **0.5185** [0.5645] | 0.8240 [0.7725] | 0.8398 [0.7912] | 0.4942 [0.5120] | **0.8688** [0.8354] | 0.8386 [0.8011] |
| **UL30-Foscarnet** | 311 (33) | **0.3903** [0.5147] | 0.8409 [0.8115] | 0.8545 [0.8320] | 0.4531 [0.4890] | **0.8930** [0.8741] | 0.8652 [0.8410] |
| **UL23-Aciclovir** | 286 (120) | **0.3775** [0.4372] | 0.7881 [0.7410] | 0.7243 [0.7055] | 0.6130 [0.6380] | 0.6347 [0.6610] | 0.6623 [0.6740] |

#### Pre-Registered Tolerance-Window Test & Implementation Reconciliation
We pre-registered a formal test of the \"Where $\\times$ What\" tolerance-window hypothesis (`docs/WHERE_WHAT_PREREGISTRATION.md`, SHA-256 `ca4d2e7ff1...`, frozen at Git commit `ded62cb` before test execution). The hypothesis predicted that on the essential polymerase (UL30), resistance requires intermediate substitution severity at constrained positions, producing negative interaction curvature ($\\beta_5 < 0$ on $\\text{Where} \\times \\text{Severity}^2$).

For full transparency, we report both the initial implementation run and the spec-exact rerun (`logs/where_what_spec_exact_results.json`):
1. **Initial Run (Implementation Deviation):** In the initial feature pipeline, raw masked marginals $\\mu_{\\text{raw}} = \\log P(\\text{mut}) - \\log P(\\text{wt}) \\le 0$ were used without negation, effectively inverting the evolutionary constraint component ($W = -Z(\\text{Constraint}) - Z(d_{\\text{site}})$). This suppressed full model fit AUROC to $0.5993$ and left $\\beta_1 = +0.1759$, with $\\beta_5 = +0.0368$ (95% CI: $[-0.1868, +0.3039]$).
2. **Spec-Exact Rerun (Pre-Registered Formulation):** Correcting the constraint feature to match the pre-registration specification ($W = Z(\\text{Constraint}) - Z(d_{\\text{site}})$ and $S = -\\text{severity}_{\\text{raw}}$):
   - Full model fit AUROC on UL30-ACV jumps to **0.8611** (and **0.9168** on UL30-FOS).
   - The linear positional constraint coefficient $\\beta_1$ is massive and strictly positive: $\\mathbf{\\beta_1 = +1.9372}$ (95% position-cluster bootstrap CI: $[+1.0555, +2.7415]$, $p < 0.0001$), confirming that structural active-site localization and evolutionary constraint dominate resistance probability.
   - The linear substitution severity coefficient is near zero: $\\beta_2 = -0.0799$ (95% CI: $[-0.4935, +0.6017]$).
   - The interaction curvature coefficient is:
     $$\\beta_5 = -0.0162, \\quad 95\\%\\text{ CI: } [-0.1876, +1.0402]$$
   While the point estimate is slightly negative ($-0.0162$), the 95% position-cluster bootstrap interval broadly spans zero (extending up to $+1.04$). Under our pre-registered decision rule (requiring the 95% CI upper bound to be $< 0$), **the tolerance-window hypothesis is not supported by the clinical data at this sample size** ($N_{\\text{pos}}=62$ on UL30-ACV).
3. **Monotonic Disruption on Thymidine Kinase:** In sharp contrast, on UL23-ACV (the dispensable kinase), the linear disruption coefficient is strictly positive: $\\mathbf{\\beta_2 = +0.4468}$ (95% CI: $[+0.1822, +0.7900]$, $p < 0.01$) with negligible curvature ($\\beta_3 = +0.0222$ [$-0.2491, +0.2513$]), confirming monotonic catalytic loss-of-function.

#### Exploratory Decomposition: The Dominance of Positional Prior
Evaluating the orthogonal components of the model reveals an essential structural finding:
1. **Positional Prior (`Where` alone):** A label-free, zero-parameter composite combining ESM-2 19-substitution positional constraint and 3D Euclidean distance to the drug-binding pocket achieves an AUROC of **0.8606** on UL30-ACV and **0.9097** on UL30-FOS.
2. **Within-Residue Severity (`What` alone):** Relative mutation severity within position achieves an AUROC of **0.4220** on UL30-ACV and **0.3935** on UL30-FOS.
3. **Combined Model:** Combining Where + Severity yields an AUROC of **0.8603** on UL30-ACV ($\\Delta = -0.0003$ over Where alone).

*Dataset Multiplicity and Identifiability Limit:* An audit of variant multiplicity demonstrates why within-position severity adds negligible predictive power: across the 332 rows of the UL30-ACV cohort, variants reside at **290 unique sequence positions**; **252 of these positions (86.9%) possess exactly one cataloged clinical variant**. Only 38 positions (13.1%) have $\\ge 2$ substitutions, and among these, **only 2 positions (0.7% of all positions) exhibit mixed resistant/sensitive labels**. On UL30-FOS, 235 of 271 positions (86.7%) have a single variant, and only 3 positions (1.1%) have mixed labels. Thus, clinical resistance catalogs in HSV-1 polymerase provide virtually no within-position phenotypic variance to identify substitution-level effects; the predictive task in current clinical data is almost entirely **hotspot (position) identification**, where a simple zero-parameter positional prior matches supervised models.

#### Cross-Herpesvirus Replication (HCMV and VZV)
To evaluate whether the empirical target-dependent radicalness dissociation generalizes beyond HSV-1, an analysis plan was specified prior to execution (`docs/CROSS_VIRUS_REPLICATION_PREREGISTRATION.md`, SHA-256 `429fd086...`, committed following the run) on the curated HerpesDRG-DB catalog (`logs/cross_virus_replication_results.json`):
- **Cytomegalovirus (HCMV) UL97 Kinase:** BLOSUM62 radicalness achieves an AUROC of **0.7005** (row bootstrap 95% CI: $[0.5943, 0.7960]$; position-cluster bootstrap 95% CI: $[0.6030, 0.7917]$ across 79 unique positions; $N=107$, 45 resistant, 62 sensitive), robustly clearing chance.
- **HCMV UL54 DNA Polymerase:** BLOSUM62 radicalness operates near chance across all three therapeutic drugs: Ganciclovir AUROC = **0.5375** ($N=230$, 76 resistant), Foscarnet AUROC = **0.4377** ($N=160$, 50 resistant), and Cidofovir AUROC = **0.5522** ($N=149$, 57 resistant). All 95% bootstrap CIs span 0.50.
- **Bootstrap Contrast Test:** Directly testing the rank-AUROC difference between kinase and polymerase in HCMV (UL97-GCV minus UL54-GCV) via 2,000 bootstrap resamples yields an observed $\\Delta = \\mathbf{+0.1630}$ (95% bootstrap CI: $[\\mathbf{+0.0341, +0.2861}]$, empirical $p = 0.0070$), confirming that the radicalness dissociation is statistically significant in HCMV.
- **Varicella-Zoster Virus (VZV):** ORF36 thymidine kinase achieves **0.7611** (95% CI: $[0.5748, 0.9219]$, $N=28$, 18 resistant, 10 sensitive). However, on VZV ORF28 (DNA polymerase), only 4 clinical resistant missense variants exist in the published literature ($N_{\\text{res}}=4, N=24$, AUROC = 0.5500, 95% CI: $[0.2273, 0.9126]$); the polymerase counter-arm is therefore underpowered and **not evaluable** for contrast testing.

*Replication Scope & Limitations:* With $n=1$ independent replication pair (HCMV UL97 vs. UL54), the radicalness contrast is directionally consistent with HSV-1. However, because UL97 loss impairs viral growth kinetics in cell culture (rendering it partially constrained rather than strictly dispensable), clinical resistance mutations concentrate within specific kinase hotspots (~460, 520, 590–607), and VZV polymerase lacks sufficient resistant mutations to establish contrast, these results represent directional consistency observed in one additional viral pair rather than a proven family-wide law.

\n\n"""
    text = text[:start_idx] + new_sec38 + text[end_idx:]
    p.write_text(text, encoding='utf-8')
    print("Section 3.8 successfully updated with audit rectifications.")
    
    # Sync artifact
    artifact_path = pathlib.Path(r'C:\Users\HP\.gemini\antigravity\brain\7e2c0f5f-723c-45b3-9f56-f3ec4abc7b17\hsv_antiviral_resistance_research_paper.md')
    shutil.copyfile(p, artifact_path)
    print("Artifact successfully synchronized.")
else:
    print(f"Error: indices not found: {start_idx}, {end_idx}")

