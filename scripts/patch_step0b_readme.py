"""
patch_step0b_readme.py
Updates README.md with Step 0b recalibrated truths.
"""
import pathlib

p = pathlib.Path('README.md')
text = p.read_text(encoding='utf-8')

# 1. Update cohort table caption and description
old_cohort = "| **Total Cohort** | - | - | **929** | **432** | **497** | 46.5% / 53.5% | **737 (79.3%)** |"
new_cohort = """| **Total Cohort** | - | - | **929 records** | **432** | **497** | 46.5% / 53.5% | **737 (79.3%)** |

*Variant Record Overlap:* The 929 records represent **630 unique amino acid substitutions** (286 in UL23, 344 in UL30). Across UL30, aciclovir ($N=332$) and foscarnet ($N=311$) share **264 sequence positions (91.0%)** and **299 exact substitutions (90.1%)**, testing the same enzyme against two drugs rather than independent biological replicates."""

if old_cohort in text:
    text = text.replace(old_cohort, new_cohort)
    print("README cohort description updated.")

# 2. Update Table 3 block in README
old_t3_marker = "### Table 3: Dual-Buffered Domain Holdout & Excess AUROC over Position-Only Baseline"
end_t3_marker = "### Variance Decomposition: \"Where\" (Positional Prior) vs. \"What\" (Mutation Severity)"

start_idx = text.find(old_t3_marker)
end_idx = text.find(end_t3_marker)

new_t3_block = """### Table 3: Recalibrated Dual-Buffered Domain Benchmark ($\pm 5$ aa, $< 12.0$ \AA)

Under a leave-one-domain-out protocol with a dual exclusion buffer, the Position-Only Random Forest operates at chance on the primary target (UL30-ACV pooled AUROC = **0.5185**). Evaluable domains strictly enforce the pre-registered $N_{\text{pos}} \ge 5$ rule (UL30-ACV: 4 domains; UL30-FOS: 2 domains [Palm and Fingers; Exo and Thumb omitted with $N_{\text{pos}} < 5$]; UL23-ACV: 4 domains).

To prevent sub-chance baselines from inflating model gains, **Excess AUROC is recalibrated against an empirical chance baseline**:
$$\Delta_{\text{excess}} = \text{AUROC}_{\text{model}} - \max(0.50, \text{AUROC}_{\text{pos\_only}})$$

95% confidence intervals are computed via domain-clustered bootstrap (recording 35 distinct multisets for $n=4$ and 3 distinct multisets for $n=2$, reflecting the small-$n$ domain limit) and complemented by non-parametric sign consistency ($k/n$ domains with $\Delta > 0$):

| Target Pair | $N_{\text{eval}}$ ($n_{\text{dom}}$) | Position-Only RF Pooled [Domain Mean] | Split Validity Status ($[0.45, 0.55]$) | Zero-Shot ESM-2 650M Domain Mean | **ESM-2 Domain Excess $\Delta$ [95% CI] (Sign)** | Zero-Shot BLOSUM62 Domain Mean | **BLOSUM62 Domain Excess $\Delta$ [95% CI] (Sign)** | Where Composite Domain Mean | **Where Domain Excess $\Delta$ [95% CI] (Sign)** |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **UL30-Aciclovir** | 248 (4) | **0.5185** [0.5645] | **CONFIRMED** (Pooled in $[0.45, 0.55]$) | 0.7439 | **+0.1794** [+0.0281, +0.2664] (3/4) | 0.4906 | -0.0739 [-0.1589, +0.0112] (1/4) | 0.7801 | **+0.2157** [+0.1116, +0.3202] (4/4) |
| **UL30-Foscarnet** | 114 (2) | **0.5087** [0.5253] | **CONFIRMED** (In Palm & Fingers) | 0.5855 | +0.0602 [-0.0833, +0.2038] (1/2) | 0.5845 | +0.0592 [-0.1213, +0.2396] (1/2) | 0.7878 | **+0.2625** [+0.2500, +0.2750] (2/2) |
| **UL23-Aciclovir** | 254 (4) | **0.3775** [0.4372] | **OUTSIDE** (Domain rate shift: 14%–75%) | 0.6885 | **+0.1885** [+0.1046, +0.2725] (4/4) | 0.6726 | **+0.1726** [+0.1238, +0.2215] (4/4) | 0.5670 | +0.0670 [-0.0360, +0.1700] (2/4) |

*Key Calibration Takeaway:* Clamping the negative control to $\max(0.50, \text{pos}_d)$ ensures models are never rewarded for beating a sub-chance baseline. BLOSUM62 achieves an excess of **+0.1726** (4/4 domains positive) on thymidine kinase (UL23), while remaining at or below chance on polymerase (UL30-ACV: -0.0739, 1/4 positive), confirming the biological dissociation under domain holdout.

---

"""

if start_idx != -1 and end_idx != -1:
    text = text[:start_idx] + new_t3_block + text[end_idx:]
    print("README Table 3 updated.")

# 3. Update beta5 diagnostic in README
old_b5_readme = "4. **Tolerance-Window Hypothesis Test:** Pre-registered test of non-monotonic curvature on UL30 ($\beta_5$) yielded a null result ($\beta_5 = -0.0162$, 95% CI: $[-0.1876, +1.0402]$, spanning zero; not supported under pre-registered stopping rule). In contrast, UL23 exhibits strictly monotonic disruption ($\beta_2 = +0.4468$, 95% CI: $[+0.1822, +0.7900]$, $p < 0.01$)."
new_b5_readme = """4. **Tolerance-Window Hypothesis Test & Bimodality Diagnostic:** Pre-registered test of non-monotonic curvature on UL30 ($\beta_5$) yielded $\beta_5 = -0.0162$ (95% CI: $[-0.1876, +1.0402]$, spanning zero). The bootstrap distribution exhibits a bimodal structure (overall skewness $= +0.9960$). A leave-one-domain-out diagnostic reveals that omitting the Palm domain ($N_{\text{pos}}=32/62$) flips $\beta_5$ to **$+0.4289$**, demonstrating that $\beta_5$ is statistically unidentifiable in retrospective clinical cohorts. The window hypothesis is neither confirmed nor refuted at current sample sizes. In contrast, UL23 exhibits strictly monotonic disruption ($\beta_2 = +0.4468$, 95% CI: $[+0.1822, +0.7900]$, $p < 0.01$)."""

if old_b5_readme in text:
    text = text.replace(old_b5_readme, new_b5_readme)
    print("README beta5 diagnostic updated.")

p.write_text(text, encoding='utf-8')
print("README.md successfully updated with Step 0b refinements.")
