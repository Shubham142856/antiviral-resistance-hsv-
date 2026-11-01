import pathlib
import os

status_path = str(pathlib.Path(__file__).resolve().parent.parent / 'STATUS.md')
with open(status_path, 'r', encoding='utf-8') as f:
    lines = f.readlines()

new_lines = []
inserted = False

for l in lines:
    new_lines.append(l)
    if 'PRIME Precedent Restored' in l and not inserted:
        addition = (
            "   - **3D Spatial-Neighborhood Holdout Pre-Registered & Evaluated (2026-11-01):** "
            "Executed the empirical structural neighborhood holdout experiment (`scripts/run_spatial_3d_holdout_experiment.py`) "
            "using cryo-EM coordinates (`8V1Q` Chain A, 1097 CA atoms, 88.6% reference coverage, 80.8% unique benchmark position coverage) "
            "and crystal coordinates (`1OF1` Chain B, 313 CA atoms, 82.4% reference coverage, 81.9% unique benchmark position coverage). "
            "Pre-registered 5 deterministic spatial clusters via KMeans (`k=5, seed=42`) into `data/processed/primary_3d_spatial_folds.tsv` "
            "(SHA-256: `0ED51D37D558530EF8B2D0C991780A23936580920E62C7ADE46A426ADFA96F4C`, logged in `logs/spatial_3d_fold_assignments_hash.log`). "
            "Key findings serialized in `logs/spatial_3d_holdout_results.json`:\n"
            "     1. *Zero-Shot PLM Invariance to Spatial Holdout:* ESM-2 650M maintains stable generalization from 1D contiguous blocks to 3D spatial neighborhoods: "
            "UL30-ACV 0.7875 -> 0.7875 (resolved core: 0.7887 -> 0.7966); UL30-FOS 0.8809 -> 0.8932 (resolved core: 0.8793 -> 0.8977); "
            "UL23-ACV 0.7397 -> 0.7303 (resolved core: 0.7255 -> 0.7132). Inter-fold variance decreases under 3D spatial holdout because spatial partitioning yields balanced classes across all 5 folds, eliminating the single-class artifact of 1D Block 0.\n"
            "     2. *Target Dissociation Invariant in 3D Space:* Chemical dissimilarity alone (BLOSUM62 radicalness) consistently clears chance on the catalytic kinase "
            "(UL23-ACV: 1D 0.6380, 3D 0.6249; resolved core: 1D 0.6732, 3D 0.6501) but sits at chance on the allosteric polymerase "
            "(UL30-ACV: 1D 0.5494, 3D 0.4908; UL30-FOS: 1D 0.5074, 3D 0.4847).\n"
            "     3. *Evolutionary Baselines in 3D Space:* Multi-herpesvirus PSSM scores 0.7189 (UL30-ACV) and 0.8334 (UL30-FOS), while positional Shannon entropy scores 0.6806 and 0.8338. "
            "ESM-2 maintains its lead over both baselines (+0.069 to +0.107 on ACV; +0.060 on FOS).\n"
            "     4. *Publication Figure 5 Generated:* High-resolution 300-DPI multi-panel visualization saved to `manuscript/figures/figure5_spatial_3d_generalization.png` and artifact storage.\n"
        )
        new_lines.append(addition)
        inserted = True

with open(status_path, 'w', encoding='utf-8') as f:
    f.writelines(new_lines)

print("Updated STATUS.md cleanly with 3D Spatial Holdout results!")

