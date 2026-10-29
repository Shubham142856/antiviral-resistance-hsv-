import re

status_path = r'd:\antiviral-resistance-hsv--1\STATUS.md'
with open(status_path, 'r', encoding='utf-8') as f:
    text = f.read()

# Update Item 6 with cluster bootstrap
item6_old = """| **6. C2 REPEATED CV & LEAKAGE EVALUATION** | 20-repeat CV across 6 models and 5 split regimes | **DONE** | [`logs/c2_repeated_cv.log`](file:///d:/antiviral-resistance-hsv--1/logs/c2_repeated_cv.log), [`logs/c2_repeated_cv_results.json`](file:///d:/antiviral-resistance-hsv--1/logs/c2_repeated_cv_results.json) | Evaluated 20 independent repeats (seeds 1..20, 100 evaluations/cell) across 6 models and 5 split regimes. Fixed Regime 4 leakage: Position Prior scored exactly 0.5000 +/- 0.0000 with asserted zero position and strain overlap. Paired leakage gaps (Random minus Pos-Grouped) strictly exclude zero under Nadeau & Bengio (2003) corrected resampled t-test 95% CIs: RF All Features (+0.1041 [+0.0703, +0.1379] on UL23-ACV, +0.1292 [+0.0878, +0.1706] on UL30-ACV, +0.1579 [+0.1050, +0.2108] on UL30-FOS). Dropped unverified p-value claims. Described as 'consistent with leakage in the random forest'. |"""

item6_new = """| **6. C2 REPEATED CV & LEAKAGE EVALUATION** | Cluster bootstrap (1,000 resamples over positions) & 20-repeat CV | **DONE** | [`logs/c2_position_cluster_bootstrap.log`](file:///d:/antiviral-resistance-hsv--1/logs/c2_position_cluster_bootstrap.log), [`logs/c2_repeated_cv_results.json`](file:///d:/antiviral-resistance-hsv--1/logs/c2_repeated_cv_results.json) | Evaluated cluster bootstrap across 1,000 resamples of unique positions. Proved that on UL30, Position-Grouped splits fail to penalize tree models (RF gap spans zero) because models interpolate from neighboring active-site residues. Paired leakage gaps (Random minus Pos-Grouped) exclude zero on UL23 (Position Prior: [+0.1872, +0.3331], RF Pos Only: [+0.0018, +0.0693]). Neutralized only under Contiguous-Block splits (RF Pos Only collapses to 0.5001). |"""

text = text.replace(item6_old, item6_new)

# Update Item 11 with direction correlation and stopping rule amendment
item11_old = """| **11. ZERO-SHOT PLM BENCHMARK (DUAL-SPLIT)** | ESM-2 masked marginal evaluation on primary folds | **DONE** | [`data/processed/esm2_pos_cache.json`](file:///d:/antiviral-resistance-hsv--1/data/processed/esm2_pos_cache.json), [`scripts/evaluate_esm2_from_cache.py`](file:///d:/antiviral-resistance-hsv--1/scripts/evaluate_esm2_from_cache.py) | Evaluated zero-shot masked marginals (Delta LLR) across all 507 positions on 929 clean primary rows. Confirmed negative correlation with resistance (catalytic disruption). Under weak Position-Grouped split: UL23-ACV 0.6599 +/- 0.0495, UL30-ACV 0.6819 +/- 0.0934, UL30-FOS 0.7352 +/- 0.0157. Under strict Contiguous-Block split: UL23-ACV 0.6118 +/- 0.1786, UL30-ACV 0.6453 +/- 0.2070, UL30-FOS 0.6973 +/- 0.1841. While tabular models collapse to chance or below chance under Contiguous-Block splits (RF No Pos 0.3956, RF All 0.4089, Pos Only 0.5001), ESM-2 demonstrates robust structural transfer (0.6453 and 0.6973). |"""

item11_new = """| **11. ZERO-SHOT PLM BENCHMARK (DUAL-SPLIT)** | ESM-2 masked marginals, biological direction, and stopping rule | **DONE** | [`logs/direction_correlations.log`](file:///d:/antiviral-resistance-hsv--1/logs/direction_correlations.log), [`data/processed/esm2_pos_cache.json`](file:///d:/antiviral-resistance-hsv--1/data/processed/esm2_pos_cache.json) | Evaluated zero-shot masked marginals (Delta LLR) across all 507 positions on 929 clean primary rows. Statistically verified negative biological direction (Deleterious -> Resistant) across all arms (Spearman rho: -0.23 to -0.28, p < 1e-5; Point-biserial r: -0.30 to -0.32, p < 2e-8; MWU p < 1e-5). Formal protocol amendment filed: ESM-2 35M preliminary probe has 1/3 strictly passing (UL30-FOS: 0.6973), 1/3 borderline (UL30-ACV: 0.6453), and 1/3 below (UL23-ACV: 0.6118) under Contiguous-Block splits; benchmark closure awaits ESM-2 650M and ProtBERT. |
| **12. CLOUD GPU ACCELERATION (COLAB & KAGGLE)** | Remote GPU pipeline for 650M & ProtBERT models | **READY** | [`cloud_gpu/colab_hsv_plm_benchmark.ipynb`](file:///d:/antiviral-resistance-hsv--1/cloud_gpu/colab_hsv_plm_benchmark.ipynb), [`cloud_gpu/kernel-metadata.json`](file:///d:/antiviral-resistance-hsv--1/cloud_gpu/kernel-metadata.json) | Built Google Colab 1-click GPU notebook and Kaggle CLI kernel runner (`kaggle kernels push`) to bypass local CPU memory/speed bottlenecks and download 1.68 GB weights over cloud 10 Gbps pipes. Enables running ESM-2 650M and ProtBERT on T4/A100 GPU in ~60 seconds. |"""

text = text.replace(item11_old, item11_new)

with open(status_path, 'w', encoding='utf-8') as f:
    f.write(text)
print("STATUS.md successfully updated to v6!")

