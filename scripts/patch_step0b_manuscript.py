"""
patch_step0b_manuscript.py
Updates manuscript/manuscript_draft.md with Step 0b recalibrated truths:
1. Reconciles Table 7 with max(0.50, pos_only) excess reference, FOS n=2 evaluable domains, distinct bootstrap multisets (35 and 3), and sign consistency.
2. Updates Abstract and Section 2.2 with 929 variant-drug records across 630 unique substitutions (and 90.1% overlap between UL30 ACV and FOS).
3. Adds the beta5 bimodality diagnostic in Section 3.8 (Palm domain influence: beta5 shifts to +0.4289 when Palm is omitted).
4. Synchronizes artifact.
"""
import pathlib
import shutil

p = pathlib.Path('manuscript/manuscript_draft.md')
text = p.read_text(encoding='utf-8')

# 1. Update Table 7 block
old_table7_marker = "**Table 7 | Dual-Buffered Domain Benchmark"
end_table7_marker = "#### Operational Split-Calibration: Excess AUROC over Position-Only Baseline"

start_idx = text.find(old_table7_marker)
end_idx = text.find(end_table7_marker)

new_table7_block = """**Table 7 | Dual-Buffered Domain Benchmark ($\pm 5$ aa, $< 12.0$ \AA): Recalibrated Predictor Performance and Excess AUROC.** 
Performance reported as pooled out-of-fold AUROC [with evaluable per-domain mean in brackets]. Evaluable domains strictly enforce the pre-registered $N_{\text{pos}} \ge 5$ rule (UL30-ACV: 4 domains; UL30-FOS: 2 domains [Palm and Fingers; Exo and Thumb omitted with $N_{\text{pos}} < 5$]; UL23-ACV: 4 domains). Excess AUROC is recalibrated against an empirical chance baseline: $\Delta_{\text{excess}} = \text{AUROC}_{\text{model}} - \max(0.50, \text{AUROC}_{\text{pos\_only}})$. 95% confidence intervals are computed via domain-clustered bootstrap (recording 35 distinct multisets for $n=4$ and 3 distinct multisets for $n=2$, reflecting the small-$n$ domain limit) and complemented by non-parametric sign consistency ($k/n$ domains with $\Delta > 0$).

| Target Pair | $N_{\text{evaluable}}$ ($n_{\text{dom}}$) | Position-Only RF Pooled [Domain Mean] | Split Validity Status ($[0.45, 0.55]$) | Zero-Shot ESM-2 650M Domain Mean | **ESM-2 Domain Excess $\Delta$ [95% CI] (Sign)** | Zero-Shot BLOSUM62 Domain Mean | **BLOSUM62 Domain Excess $\Delta$ [95% CI] (Sign)** | Where Composite Domain Mean | **Where Domain Excess $\Delta$ [95% CI] (Sign)** |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **UL30-Aciclovir** | 248 (4) | **0.5185** [0.5645] | **CONFIRMED** (Pooled in $[0.45, 0.55]$) | 0.7439 | **+0.1794** [+0.0281, +0.2664] (3/4) | 0.4906 | -0.0739 [-0.1589, +0.0112] (1/4) | 0.7801 | **+0.2157** [+0.1116, +0.3202] (4/4) |
| **UL30-Foscarnet** | 114 (2) | **0.5087** [0.5253] | **CONFIRMED** (In Palm & Fingers) | 0.5855 | +0.0602 [-0.0833, +0.2038] (1/2) | 0.5845 | +0.0592 [-0.1213, +0.2396] (1/2) | 0.7878 | **+0.2625** [+0.2500, +0.2750] (2/2) |
| **UL23-Aciclovir** | 254 (4) | **0.3775** [0.4372] | **OUTSIDE** (Domain rate shift: 14%–75%) | 0.6885 | **+0.1885** [+0.1046, +0.2725] (4/4) | 0.6726 | **+0.1726** [+0.1238, +0.2215] (4/4) | 0.5670 | +0.0670 [-0.0360, +0.1700] (2/4) |

*Calibration & Small-n Domain Notes:*
1. **Split Validity:** Operational split validity ($[0.45, 0.55]$ on Position-Only RF) is confirmed for UL30-Aciclovir pooled (0.5185) and within the 2 evaluable domains of UL30-Foscarnet (0.5253). For UL23-Aciclovir (0.3775), pooled position-only AUROC falls below 0.45 because cross-domain pooling aggregates across structural domains with divergent resistant rates (14% in C-term vs. 75% in Nucleoside Pocket). Clamping the negative control to $\max(0.50, \text{pos}_d)$ ensures models are never rewarded for beating a sub-chance baseline.
2. **Small-n Domain Bootstrap Discrete Bound:** Resampling with replacement across $n=4$ evaluable domains generates exactly $\binom{7}{4} = 35$ distinct multiset combinations; for $n=2$ domains, it generates only $\binom{3}{2} = 3$ combinations. Consequently, domain bootstrap confidence intervals are inherently discrete and descriptive; non-parametric sign consistency ($k/n$) provides the primary evidence of robust directional gain.
3. **Biological Dissociation Invariant:** With the sub-chance baseline clamped to 0.50, BLOSUM62 radicalness achieves a positive excess of **+0.1726** (4/4 domains positive) on thymidine kinase (UL23), while operating at or below chance on polymerase (UL30-ACV: -0.0739, only 1/4 domains positive).

"""

if start_idx != -1 and end_idx != -1:
    text = text[:start_idx] + new_table7_block + text[end_idx:]
    print("Table 7 successfully updated with Step 0b recalibrated metrics.")
else:
    print("Warning: Table 7 markers not found.")

# 2. Update beta5 bimodality diagnostic in Section 3.8
old_b5_text = "While the point estimate is slightly negative ($-0.0162$), the 95% position-cluster bootstrap interval broadly spans zero (extending up to $+1.04$). Under our pre-registered decision rule (requiring the 95% CI upper bound to be $< 0$), **the tolerance-window hypothesis is not supported by the clinical data at this sample size** ($N_{\\text{pos}}=62$ on UL30-ACV)."
new_b5_text = """While the point estimate is slightly negative ($-0.0162$), the 95% position-cluster bootstrap interval broadly spans zero (extending up to $+1.0402$). An inspection of the bootstrap distribution (`logs/beta5_bootstrap_distribution.json`) reveals a pronounced bimodal structure: 41.9% of resamples cluster around $\beta_5 \approx 0$, while 58.1% of resamples shift into a secondary positive mode ($\beta_5 \approx +0.40$ to $+0.80$, overall skewness $= +0.9960$). 

A leave-one-domain-out influence diagnostic (`logs/step0b_recalibrated_benchmarks.json`) identifies the structural driver: omitting the Palm domain ($N=111$, containing 32 of the 62 resistant events) shifts $\beta_5$ from $-0.0162$ to **$+0.4289$**, whereas omitting any other domain leaves $\beta_5$ near zero ($-0.012$ to $-0.086$). Because clinical resistant events are heavily concentrated in a single catalytic domain, $\beta_5$ is **statistically unidentifiable** in retrospective clinical cohorts. Under our pre-registered decision rule (requiring the 95% CI upper bound to be $< 0$), **the tolerance-window hypothesis is unidentifiable and neither confirmed nor refuted at current clinical sample sizes**."""

if old_b5_text in text:
    text = text.replace(old_b5_text, new_b5_text)
    print("Beta5 bimodality diagnostic successfully updated.")
else:
    print("Warning: old beta5 text not found directly.")

# 3. Update variant-drug records vs unique substitutions and position overlap
old_cohort_text = "The primary cohort comprises 929 rows:"
new_cohort_text = """The primary cohort comprises **929 single-site missense variant-drug records** across **630 unique amino acid substitutions** (286 unique substitutions in UL23; 344 unique substitutions in UL30):
- **UL23-Aciclovir:** 286 variant-drug records (120 Resistant, 166 Sensitive; 210 unique positions; 286 unique substitutions).
- **UL30-Aciclovir:** 332 variant-drug records (62 Resistant, 270 Sensitive; 290 unique positions; 332 unique substitutions).
- **UL30-Foscarnet:** 311 variant-drug records (33 Resistant, 278 Sensitive; 271 unique positions; 311 unique substitutions).

*Substantial Variant Overlap Across Polymerase Drugs:* UL30-ACV and UL30-FOS evaluate the **same enzyme** against two different inhibitor classes. Crucially, they share **264 sequence positions (91.0% of ACV positions; 97.4% of FOS positions)** and **299 exact substitutions (90.1% of ACV substitutions; 96.1% of FOS substitutions)**. Foscarnet therefore does not represent an independent biological validation arm, but rather a cross-drug comparison on an essentially identical variant set."""

if old_cohort_text in text:
    # Find up to next subheading
    sub_start = text.find(old_cohort_text)
    sub_end = text.find("The complete dataset is frozen and tracked by its SHA-256 digest")
    if sub_start != -1 and sub_end != -1:
        text = text[:sub_start] + new_cohort_text + "\n\n" + text[sub_end:]
        print("Cohort description successfully updated with overlap realities.")

p.write_text(text, encoding='utf-8')

# Sync artifact
artifact_path = pathlib.Path(r'C:\Users\HP\.gemini\antigravity\brain\7e2c0f5f-723c-45b3-9f56-f3ec4abc7b17\hsv_antiviral_resistance_research_paper.md')
shutil.copyfile(p, artifact_path)
print("Artifact successfully synchronized.")
