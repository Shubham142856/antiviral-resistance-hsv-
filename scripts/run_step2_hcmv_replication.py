"""
run_step2_hcmv_replication.py
==============================
Step 2: Cross-Herpesvirus Replication of Structural Proximity and Radicalness Dissociation.
Anchored strictly to frozen pre-registration: docs/STEP2_HCMV_REPLICATION_PREREGISTRATION.md
(Commit: da59f55)

Targets:
  - Primary Arm:   HCMV UL54 - Ganciclovir (Essential Replicative DNA Polymerase)
  - Secondary Arm: HCMV UL54 - Foscarnet (Essential Replicative DNA Polymerase)
  - Secondary Arm: HCMV UL54 - Cidofovir (Essential Replicative DNA Polymerase)
  - Counter-Arm:   HCMV UL97 - Ganciclovir (Dispensable Phosphotransferase / Kinase)

Evaluates:
  1. Positivity Check: M mixed studies, discordant pair count, top study share, distinct variants.
  2. S0 Within-Study Control (Primary):
     - Pair-weighted AUROC on same-study pairs
     - Study-clustered bootstrap (B=2000, seed 42) for 95% CIs
     - Study-equal unweighted mean of per-study AUROCs
     - Leave-One-Study-Out (LOSO) stability analysis
  3. Feature Sets:
     - BLOSUM62 Radicalness: - BLOSUM62(wt, mut)
     - Option A (Primary Where): Active-site 3D proximity on Swiss-Model P08546_swissmodel_8oj7.pdb
     - Option B (Sensitivity Where): Homology-projected active-site distance from HSV-1 8EXX
  4. Secondary Stratum S2: Quantitative numeric fold-change measurements only.
  5. Uncontrolled Naive Pooled D1.
  6. Empirical Contrast Tests:
     - Delta_BLOSUM: UL97-GCV vs UL54-GCV
     - Delta_Polymerase: UL54-FOS Proximity vs BLOSUM62

Outputs:
  - logs/step2_hcmv_replication_results.json
"""

import os
import re
import json
import warnings
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from Bio.PDB import PDBParser, MMCIFParser
from Bio.Align import PairwiseAligner, substitution_matrices

warnings.filterwarnings('ignore')

import sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

# 1. BLOSUM62 Matrix
AA_LIST = list("ARNDCQEGHILKMFPSTWYV")
BLOSUM62 = {
    'A': {'A': 4, 'R': -1, 'N': -2, 'D': -2, 'C': 0, 'Q': -1, 'E': -1, 'G': 0, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 0, 'W': -3, 'Y': -2, 'V': 0},
    'R': {'A': -1, 'R': 5, 'N': 0, 'D': -2, 'C': -3, 'Q': 1, 'E': 0, 'G': -2, 'H': 0, 'I': -3, 'L': -2, 'K': 2, 'M': -1, 'F': -3, 'P': -2, 'S': -1, 'T': -1, 'W': -3, 'Y': -2, 'V': -3},
    'N': {'A': -2, 'R': 0, 'N': 6, 'D': 1, 'C': -3, 'Q': 0, 'E': 0, 'G': 0, 'H': 1, 'I': -3, 'L': -3, 'K': 0, 'M': -2, 'F': -3, 'P': -2, 'S': 1, 'T': 0, 'W': -4, 'Y': -2, 'V': -3},
    'D': {'A': -2, 'R': -2, 'N': 1, 'D': 6, 'C': -3, 'Q': 0, 'E': 2, 'G': -1, 'H': -1, 'I': -3, 'L': -4, 'K': -1, 'M': -3, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -4, 'Y': -3, 'V': -3},
    'C': {'A': 0, 'R': -3, 'N': -3, 'D': -3, 'C': 9, 'Q': -3, 'E': -4, 'G': -3, 'H': -3, 'I': -1, 'L': -1, 'K': -3, 'M': -1, 'F': -2, 'P': -3, 'S': -1, 'T': -1, 'W': -2, 'Y': -2, 'V': -1},
    'Q': {'A': -1, 'R': 1, 'N': 0, 'D': 0, 'C': -3, 'Q': 5, 'E': 2, 'G': -2, 'H': 0, 'I': -3, 'L': -2, 'K': 1, 'M': 0, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -2, 'Y': -1, 'V': -2},
    'E': {'A': -1, 'R': 0, 'N': 0, 'D': 2, 'C': -4, 'Q': 2, 'E': 5, 'G': -2, 'H': 0, 'I': -3, 'L': -3, 'K': 1, 'M': -2, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2},
    'G': {'A': 0, 'R': -2, 'N': 0, 'D': -1, 'C': -3, 'Q': -2, 'E': -2, 'G': 6, 'H': -2, 'I': -4, 'L': -4, 'K': -2, 'M': -3, 'F': -3, 'P': -2, 'S': 0, 'T': -2, 'W': -2, 'Y': -3, 'V': -3},
    'H': {'A': -2, 'R': 0, 'N': 1, 'D': -1, 'C': -3, 'Q': 0, 'E': 0, 'G': -2, 'H': 8, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -1, 'P': -2, 'S': -1, 'T': -2, 'W': -2, 'Y': 2, 'V': -3},
    'I': {'A': -1, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -3, 'E': -3, 'G': -4, 'H': -3, 'I': 4, 'L': 2, 'K': -3, 'M': 1, 'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -3, 'Y': -1, 'V': 3},
    'L': {'A': -1, 'R': -2, 'N': -3, 'D': -4, 'C': -1, 'Q': -2, 'E': -3, 'G': -4, 'H': -3, 'I': 2, 'L': 4, 'K': -2, 'M': 2, 'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -2, 'Y': -1, 'V': 1},
    'K': {'A': -1, 'R': 2, 'N': 0, 'D': -1, 'C': -3, 'Q': 1, 'E': 1, 'G': -2, 'H': -1, 'I': -3, 'L': -2, 'K': 5, 'M': -1, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2},
    'M': {'A': -1, 'R': -1, 'N': -2, 'D': -3, 'C': -1, 'Q': 0, 'E': -2, 'G': -3, 'H': -2, 'I': 1, 'L': 2, 'K': -1, 'M': 5, 'F': 0, 'P': -2, 'S': -1, 'T': -1, 'W': -1, 'Y': -1, 'V': 1},
    'F': {'A': -2, 'R': -3, 'N': -3, 'D': -3, 'C': -2, 'Q': -3, 'E': -3, 'G': -3, 'H': -1, 'I': 0, 'L': 0, 'K': -3, 'M': 0, 'F': 6, 'P': -4, 'S': -2, 'T': -2, 'W': 1, 'Y': 3, 'V': -1},
    'P': {'A': -1, 'R': -2, 'N': -2, 'D': -1, 'C': -3, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -4, 'P': 7, 'S': -1, 'T': -1, 'W': -4, 'Y': -3, 'V': -2},
    'S': {'A': 1, 'R': -1, 'N': 1, 'D': 0, 'C': -1, 'Q': 0, 'E': 0, 'G': 0, 'H': -1, 'I': -2, 'L': -2, 'K': 0, 'M': -1, 'F': -2, 'P': -1, 'S': 4, 'T': 1, 'W': -3, 'Y': -2, 'V': -2},
    'T': {'A': 0, 'R': -1, 'N': 0, 'D': -1, 'C': -1, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 5, 'W': -2, 'Y': -2, 'V': 0},
    'W': {'A': -3, 'R': -3, 'N': -4, 'D': -4, 'C': -2, 'Q': -2, 'E': -3, 'G': -2, 'H': -2, 'I': -3, 'L': -2, 'K': -3, 'M': -1, 'F': 1, 'P': -4, 'S': -3, 'T': -2, 'W': 11, 'Y': 2, 'V': -3},
    'Y': {'A': -2, 'R': -2, 'N': -2, 'D': -3, 'C': -2, 'Q': -1, 'E': -2, 'G': -3, 'H': 2, 'I': -1, 'L': -1, 'K': -2, 'M': -1, 'F': 3, 'P': -3, 'S': -2, 'T': -2, 'W': 2, 'Y': 7, 'V': -1},
    'V': {'A': 0, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -2, 'E': -2, 'G': -3, 'H': -3, 'I': 3, 'L': 1, 'K': -2, 'M': 1, 'F': -1, 'P': -2, 'S': -2, 'T': 0, 'W': -3, 'Y': -1, 'V': 4}
}

# 2. Extract Structural Distance Features (Option A: Swiss-Model)
print("Parsing Option A: Swiss-Model HCMV UL54 structure (P08546_swissmodel_8oj7.pdb)...")
pdb_parser = PDBParser(QUIET=True)
s_sm = pdb_parser.get_structure('UL54', os.path.join(ROOT, 'data', 'structures', 'P08546_swissmodel_8oj7.pdb'))
m_sm = list(s_sm.get_models())[0]

with open(os.path.join(ROOT, 'logs', 'hcmv_ul54_active_site_mapping.json')) as f:
    mapping_data = json.load(f)
as_positions = set(mapping_data['hcmv_active_site_positions'])

as_atoms = []
ca_coords_sm = {}
for chain in m_sm:
    for res in chain:
        r_id = res.id[1]
        if r_id in as_positions:
            for a in res:
                as_atoms.append(a.get_coord())
        if 'CA' in res:
            ca_coords_sm[r_id] = res['CA'].get_coord()

as_atoms = np.array(as_atoms)
dists_option_a = {p: float(np.min(np.linalg.norm(as_atoms - c, axis=1))) for p, c in ca_coords_sm.items()}
max_dist_a = max(dists_option_a.values())
print(f"  Option A: {len(dists_option_a)} residue C-alpha coordinates mapped. Max distance: {max_dist_a:.2f} A.")

# 3. Extract Sensitivity Distances (Option B: HSV-1 8EXX Projected)
print("Parsing Option B: HSV-1 8EXX cryo-EM coordinates mapped via BLOSUM62 alignment...")
cif_parser = MMCIFParser(QUIET=True)
s_8exx = cif_parser.get_structure('8EXX', os.path.join(ROOT, 'data', 'structures', '8exx.cif'))
m_8exx = list(s_8exx.get_models())[0]

hsv_as_atoms = []
for res in m_8exx.get_residues():
    if res.get_resname() in ['PPF', 'MG']:
        for a in res.get_atoms(): hsv_as_atoms.append(a.get_coord())
    if res.get_parent().id == 'A' and res.id[1] in [717, 888, 889]:
        for a in res.get_atoms(): hsv_as_atoms.append(a.get_coord())
hsv_as_atoms = np.array(hsv_as_atoms)

hsv_ca = {}
for res in m_8exx['A'].get_residues():
    if 'CA' in res: hsv_ca[res.id[1]] = res['CA'].get_coord()

hsv_dists = {p: float(np.min(np.linalg.norm(hsv_as_atoms - c, axis=1))) for p, c in hsv_ca.items()}

# Pairwise alignment
from scripts.map_active_site_to_hcmv import get_gene_seq
seq_hsv = get_gene_seq('NC_001806.2', 'UL30')
seq_hcmv = get_gene_seq('NC_006273.2', 'UL54')
aligner = PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
aligner.open_gap_score = -10.0
aligner.extend_gap_score = -0.5
aligner.mode = 'global'
aln = aligner.align(seq_hsv, seq_hcmv)[0]

hsv_idx, hcmv_idx = 0, 0
hcmv_to_hsv = {}
for r_ch, q_ch in zip(aln[0], aln[1]):
    if r_ch != '-': hsv_idx += 1
    if q_ch != '-': hcmv_idx += 1
    if r_ch != '-' and q_ch != '-': hcmv_to_hsv[hcmv_idx] = hsv_idx

dists_option_b = {p: hsv_dists.get(hcmv_to_hsv.get(p), np.nan) for p in range(1, len(seq_hcmv)+1)}
max_dist_b = np.nanmax(list(dists_option_b.values()))
print(f"  Option B: {len([p for p, d in dists_option_b.items() if not np.isnan(d)])} residues mapped to HSV-1.")

# 4. Load and Filter HerpesDRG HCMV Data
print("\nLoading and phenotyping HerpesDRG HCMV cohort...")
raw_df = pd.read_csv(os.path.join(ROOT, 'data', 'raw', 'herpesdrg-db', 'herpesdrg-db.tsv'), sep='\t')
hcmv_df = raw_df[raw_df['virus'] == 'HCMV'].copy()

def parse_mut(s):
    m = re.match(r'^([A-Z])(\d+)([A-Z])$', str(s).strip())
    return (m.group(1), int(m.group(2)), m.group(3)) if m else (None, None, None)

p = [parse_mut(x) for x in hcmv_df['aa_change']]
hcmv_df['wt'] = [x[0] for x in p]
hcmv_df['pos'] = [x[1] for x in p]
hcmv_df['mut'] = [x[2] for x in p]
hcmv_df = hcmv_df[hcmv_df['wt'].notna()].copy()

def get_study_id(row):
    doi = str(row['ref_doi']).strip()
    if doi and doi not in ['nan', 'None', '']: return f"DOI:{doi}"
    t = str(row['ref_title']).strip()
    if t and t not in ['nan', 'None', '']: return f"TITLE:{t[:45]}"
    return "Unassigned"

hcmv_df['study_id'] = [get_study_id(r) for _, r in hcmv_df.iterrows()]

def is_res(val, thresh=2.0):
    s = str(val).strip()
    if s == 'Resistant': return 1
    if s in ['Polymorphism', 'Sensitive']: return 0
    try: return 1 if float(s) >= thresh else 0
    except: return None

def is_numeric_fc(val):
    if pd.isna(val): return False
    s = str(val).replace('>', '').replace('<', '').strip()
    try:
        float(s)
        return True
    except:
        return False

# Target configurations
targets = [
    {
        'key': 'HCMV_UL54_Ganciclovir',
        'gene': 'UL54',
        'drug': 'Ganciclovir',
        'role': 'DNA Polymerase (Essential Replicative Engine)',
        'arm_type': 'Primary Replication Arm'
    },
    {
        'key': 'HCMV_UL54_Foscarnet',
        'gene': 'UL54',
        'drug': 'Foscarnet',
        'role': 'DNA Polymerase (Essential Replicative Engine)',
        'arm_type': 'Secondary Descriptive Arm'
    },
    {
        'key': 'HCMV_UL54_Cidofovir',
        'gene': 'UL54',
        'drug': 'Cidofovir',
        'role': 'DNA Polymerase (Essential Replicative Engine)',
        'arm_type': 'Secondary Descriptive Arm'
    },
    {
        'key': 'HCMV_UL97_Ganciclovir',
        'gene': 'UL97',
        'drug': 'Ganciclovir',
        'role': 'Phosphotransferase (Dispensable Prodrug Activator)',
        'arm_type': 'Dispensable Counter-Arm'
    }
]

# Helper function for Mann-Whitney U AUROC on pairs
def get_study_pair_stats(df_sub, score_col):
    st_groups = df_sub.groupby('study_id')
    stats = {}
    for st_id, grp in st_groups:
        res_scores = grp.loc[grp['y'] == 1, score_col].values
        sens_scores = grp.loc[grp['y'] == 0, score_col].values
        n_r = len(res_scores)
        n_s = len(sens_scores)
        if n_r > 0 and n_s > 0:
            n_p = n_r * n_s
            diffs = res_scores[:, None] - sens_scores[None, :]
            concordant = float(np.sum(diffs > 0) + 0.5 * np.sum(diffs == 0))
            stats[st_id] = {
                'n_res': n_r,
                'n_sens': n_s,
                'n_pairs': n_p,
                'concordant': concordant,
                'auroc': concordant / n_p
            }
    return stats

def compute_pair_weighted_auroc(stats_dict):
    tot_p = sum(s['n_pairs'] for s in stats_dict.values())
    tot_c = sum(s['concordant'] for s in stats_dict.values())
    return (tot_c / tot_p) if tot_p > 0 else np.nan, tot_p

def compute_study_mean_auroc(stats_dict):
    if not stats_dict: return np.nan
    return float(np.mean([s['auroc'] for s in stats_dict.values()]))

# Shared Study-Clustered Bootstrap Engine
RNG = np.random.RandomState(42)
N_BOOT = 2000

results = {
    'protocol': {
        'document': 'docs/STEP2_HCMV_REPLICATION_PREREGISTRATION.md',
        'commit': 'da59f55',
        'n_boot': N_BOOT,
        'seed': 42
    },
    'targets': {},
    'contrasts': {}
}

# Collect target-specific dataframes for joint bootstrap
target_dfs = {}

print("\n" + "="*110)
print("STEP 2: HCMV REPLICATION BENCHMARK (WITHIN-STUDY S0 & EVIDENCE LADDER)")
print("="*110)

for tgt in targets:
    t_key = tgt['key']
    gene = tgt['gene']
    drug = tgt['drug']
    
    sub = hcmv_df[hcmv_df['gene'] == gene].copy()
    sub['y'] = [is_res(x) for x in sub[drug]]
    sub = sub[sub['y'].notna()].copy()
    sub['y'] = sub['y'].astype(int)
    sub['is_num_fc'] = [is_numeric_fc(x) for x in sub[drug]]
    
    # Deduplicate per (aa_change, study_id)
    sub = sub.drop_duplicates(subset=['aa_change', 'study_id']).copy().reset_index(drop=True)
    
    # Features
    sub['blosum_rad'] = [- BLOSUM62.get(r['wt'], {}).get(r['mut'], 0) for _, r in sub.iterrows()]
    if gene == 'UL54':
        def get_a(pos):
            v = dists_option_a.get(pos)
            return - float(v) if v is not None and not np.isnan(v) else - float(max_dist_a)
        def get_b(pos):
            v = dists_option_b.get(pos)
            return - float(v) if v is not None and not np.isnan(v) else - float(max_dist_b)
            
        sub['where_option_a'] = [get_a(r['pos']) for _, r in sub.iterrows()]
        sub['where_option_b'] = [get_b(r['pos']) for _, r in sub.iterrows()]
    
    target_dfs[t_key] = sub
    
    # --- A. Positivity Check ---
    st_counts = sub.groupby('study_id')['y'].agg(
        total='count',
        res=lambda x: int((x == 1).sum()),
        sens=lambda x: int((x == 0).sum())
    )
    mixed_st = st_counts[(st_counts['res'] > 0) & (st_counts['sens'] > 0)].copy()
    mixed_st['pairs'] = mixed_st['res'] * mixed_st['sens']
    tot_pairs = int(mixed_st['pairs'].sum())
    top_share = float((mixed_st['pairs'].max() / tot_pairs) * 100) if tot_pairs > 0 else 0.0
    
    res_in_mixed = sub[sub['study_id'].isin(mixed_st.index) & (sub['y'] == 1)]
    sens_in_mixed = sub[sub['study_id'].isin(mixed_st.index) & (sub['y'] == 0)]
    distinct_res = res_in_mixed['aa_change'].nunique()
    distinct_sens = sens_in_mixed['aa_change'].nunique()
    
    rule_pairs = (tot_pairs >= 100)
    rule_mixed = (len(mixed_st) >= 5)
    rule_share = (top_share <= 50.0)
    pos_pass = rule_pairs and rule_mixed and rule_share
    
    top5_studies = []
    for st_id, r in mixed_st.sort_values(by='pairs', ascending=False).head(5).iterrows():
        top5_studies.append({
            'study_id': st_id,
            'resistant': int(r['res']),
            'sensitive': int(r['sens']),
            'pairs': int(r['pairs']),
            'share_pct': round(float(r['pairs'] / tot_pairs * 100), 2)
        })
        
    positivity_info = {
        'total_cohort_rows': int(len(sub)),
        'unique_mutations_overall': int(sub['aa_change'].nunique()),
        'total_studies': int(sub['study_id'].nunique()),
        'mixed_study_count': int(len(mixed_st)),
        'total_same_study_pairs': tot_pairs,
        'top_1_study_share_pct': round(top_share, 2),
        'distinct_res_variants_in_mixed': distinct_res,
        'distinct_sens_variants_in_mixed': distinct_sens,
        'top_5_studies_by_pairs': top5_studies,
        'rules': {
            'pairs_ge_100': rule_pairs,
            'mixed_ge_5': rule_mixed,
            'top_share_le_50': rule_share,
            'overall_status': 'PASS' if pos_pass else 'FAIL'
        }
    }
    
    # --- B. S0 Evaluations ---
    features_to_eval = [('BLOSUM62 Radicalness', 'blosum_rad')]
    if gene == 'UL54':
        features_to_eval.append(('Where (Option A: Swiss-Model)', 'where_option_a'))
        features_to_eval.append(('Where (Option B: HSV-1 Mapped)', 'where_option_b'))
        
    feature_results = {}
    for f_label, f_col in features_to_eval:
        st_stats = get_study_pair_stats(sub, f_col)
        pw_auc, n_p = compute_pair_weighted_auroc(st_stats)
        sm_auc = compute_study_mean_auroc(st_stats)
        
        # Study-clustered bootstrap
        mixed_study_ids = np.array(list(st_stats.keys()))
        boot_pw = []
        boot_sm = []
        for _ in range(N_BOOT):
            boot_ids = RNG.choice(mixed_study_ids, size=len(mixed_study_ids), replace=True)
            b_pairs = sum(st_stats[s]['n_pairs'] for s in boot_ids)
            b_conc = sum(st_stats[s]['concordant'] for s in boot_ids)
            if b_pairs > 0:
                boot_pw.append(b_conc / b_pairs)
                boot_sm.append(np.mean([st_stats[s]['auroc'] for s in boot_ids]))
        
        ci_pw = [float(np.percentile(boot_pw, 2.5)), float(np.percentile(boot_pw, 97.5))]
        ci_sm = [float(np.percentile(boot_sm, 2.5)), float(np.percentile(boot_sm, 97.5))]
        
        # LOSO analysis
        loso_pw = []
        for s_omit in mixed_study_ids:
            rem = [s for s in mixed_study_ids if s != s_omit]
            r_pairs = sum(st_stats[s]['n_pairs'] for s in rem)
            r_conc = sum(st_stats[s]['concordant'] for s in rem)
            if r_pairs > 0:
                loso_pw.append(r_conc / r_pairs)
        loso_range = [float(min(loso_pw)), float(max(loso_pw))] if loso_pw else [np.nan, np.nan]
        
        # Pooled D1
        dedup_sub = sub.drop_duplicates(subset=['aa_change']).copy()
        y_d1 = dedup_sub['y'].values
        sc_d1 = dedup_sub[f_col].values
        d1_auc = float(roc_auc_score(y_d1, sc_d1))
        # Bootstrap on D1
        d1_boots = []
        n_d1 = len(y_d1)
        for _ in range(N_BOOT):
            b_idx = RNG.choice(n_d1, size=n_d1, replace=True)
            if len(np.unique(y_d1[b_idx])) > 1:
                d1_boots.append(roc_auc_score(y_d1[b_idx], sc_d1[b_idx]))
        ci_d1 = [float(np.percentile(d1_boots, 2.5)), float(np.percentile(d1_boots, 97.5))]
        
        # S2 Stratum (numeric fold change)
        s2_sub = sub[sub['is_num_fc']].copy()
        s2_stats = get_study_pair_stats(s2_sub, f_col)
        s2_pw, s2_pairs = compute_pair_weighted_auroc(s2_stats)
        s2_sm = compute_study_mean_auroc(s2_stats)
        s2_dedup = s2_sub.drop_duplicates(subset=['aa_change'])
        if len(np.unique(s2_dedup['y'])) > 1:
            s2_d1 = float(roc_auc_score(s2_dedup['y'], s2_dedup[f_col]))
        else:
            s2_d1 = np.nan
            
        # Decision status
        if f_label == 'Where (Option A: Swiss-Model)':
            if pw_auc >= 0.65 and ci_pw[0] > 0.55:
                decision = "PASS (Replicated Proximity)"
            elif pw_auc >= 0.55 and ci_pw[0] > 0.50:
                decision = "INCONCLUSIVE (Moderate Signal / Wide CI)"
            else:
                decision = "NOT SUPPORTED (Operates Near or Below Chance)"
        elif 'BLOSUM' in f_label:
            if gene == 'UL97':
                decision = "PASS (Replicated Kinase Radicalness)" if pw_auc >= 0.60 and ci_pw[0] > 0.50 else "NOT SUPPORTED"
            else:
                decision = "CONFIRMED (Polymerase Radicalness Chance/Suppressed)" if ci_pw[0] <= 0.55 else "ACTIVE"
        else:
            decision = "SENSITIVITY CONTROL"
            
        feature_results[f_label] = {
            's0_pair_weighted_auroc': round(pw_auc, 4),
            's0_pair_weighted_ci_95': [round(x, 4) for x in ci_pw],
            's0_study_mean_auroc': round(sm_auc, 4),
            's0_study_mean_ci_95': [round(x, 4) for x in ci_sm],
            's0_loso_range': [round(x, 4) for x in loso_range],
            's2_numeric_fc_pair_weighted_auroc': round(s2_pw, 4) if not np.isnan(s2_pw) else None,
            's2_numeric_fc_pairs': s2_pairs,
            's2_numeric_fc_pooled_auroc': round(s2_d1, 4) if not np.isnan(s2_d1) else None,
            'd1_pooled_auroc': round(d1_auc, 4),
            'd1_pooled_ci_95': [round(x, 4) for x in ci_d1],
            'decision': decision
        }
        
    results['targets'][t_key] = {
        'gene': gene,
        'drug': drug,
        'role': tgt['role'],
        'arm_type': tgt['arm_type'],
        'positivity': positivity_info,
        'features': feature_results
    }
    
    print(f"\nTarget: {gene} - {drug} ({tgt['arm_type']})")
    print(f"  Positivity: M={len(mixed_st)}, Pairs={tot_pairs}, Top-1 Study Share={top_share:.1f}% -> [{positivity_info['rules']['overall_status']}]")
    print(f"  {'Feature':<32} | {'S0 Pair-Weighted [95% CI]':<26} | {'Study-Mean':<12} | {'Pooled D1':<12} | {'Verdict'}")
    print("  " + "-" * 105)
    for f_label, fr in feature_results.items():
        ci_str = f"[{fr['s0_pair_weighted_ci_95'][0]:.4f}, {fr['s0_pair_weighted_ci_95'][1]:.4f}]"
        pw_str = f"{fr['s0_pair_weighted_auroc']:.4f} {ci_str}"
        print(f"  {f_label:<32} | {pw_str:<26} | {fr['s0_study_mean_auroc']:.4f}       | {fr['d1_pooled_auroc']:.4f}       | {fr['decision']}")


# --- C. Contrast Tests ---
print("\n" + "="*110)
print("PRE-REGISTERED CONTRAST TESTS")
print("="*110)

# Contrast 1: Kinase (UL97-GCV) vs Polymerase (UL54-GCV) BLOSUM62 Radicalness
df_ul97 = target_dfs['HCMV_UL97_Ganciclovir']
df_ul54 = target_dfs['HCMV_UL54_Ganciclovir']

st_97 = get_study_pair_stats(df_ul97, 'blosum_rad')
st_54 = get_study_pair_stats(df_ul54, 'blosum_rad')

obs_pw_97, _ = compute_pair_weighted_auroc(st_97)
obs_pw_54, _ = compute_pair_weighted_auroc(st_54)
delta_obs_gcv = obs_pw_97 - obs_pw_54

studies_97 = np.array(list(st_97.keys()))
studies_54 = np.array(list(st_54.keys()))

boot_deltas_gcv = []
for _ in range(N_BOOT):
    b97 = RNG.choice(studies_97, size=len(studies_97), replace=True)
    b54 = RNG.choice(studies_54, size=len(studies_54), replace=True)
    
    p97 = sum(st_97[s]['n_pairs'] for s in b97)
    c97 = sum(st_97[s]['concordant'] for s in b97)
    
    p54 = sum(st_54[s]['n_pairs'] for s in b54)
    c54 = sum(st_54[s]['concordant'] for s in b54)
    
    if p97 > 0 and p54 > 0:
        boot_deltas_gcv.append((c97 / p97) - (c54 / p54))

ci_delta_gcv = [float(np.percentile(boot_deltas_gcv, 2.5)), float(np.percentile(boot_deltas_gcv, 97.5))]
p_emp_gcv = float(np.mean(np.array(boot_deltas_gcv) <= 0))

print(f"Contrast 1: BLOSUM62 Radicalness Dissociation (UL97-GCV vs UL54-GCV)")
print(f"  UL97-GCV S0: {obs_pw_97:.4f} | UL54-GCV S0: {obs_pw_54:.4f}")
print(f"  Observed Delta: +{delta_obs_gcv:.4f} [95% CI: {ci_delta_gcv[0]:.4f}, {ci_delta_gcv[1]:.4f}] | Empirical p = {p_emp_gcv:.4f}")

# Contrast 2: Polymerase (UL54-FOS) Where vs BLOSUM62
df_ul54_fos = target_dfs['HCMV_UL54_Foscarnet']
st_fos_where = get_study_pair_stats(df_ul54_fos, 'where_option_a')
st_fos_blosum = get_study_pair_stats(df_ul54_fos, 'blosum_rad')

obs_fos_where, _ = compute_pair_weighted_auroc(st_fos_where)
obs_fos_blosum, _ = compute_pair_weighted_auroc(st_fos_blosum)
delta_obs_fos = obs_fos_where - obs_fos_blosum

studies_fos = np.array(list(st_fos_where.keys()))
boot_deltas_fos = []
for _ in range(N_BOOT):
    b_ids = RNG.choice(studies_fos, size=len(studies_fos), replace=True)
    pw = sum(st_fos_where[s]['n_pairs'] for s in b_ids)
    cw = sum(st_fos_where[s]['concordant'] for s in b_ids)
    pb = sum(st_fos_blosum[s]['n_pairs'] for s in b_ids)
    cb = sum(st_fos_blosum[s]['concordant'] for s in b_ids)
    if pw > 0 and pb > 0:
        boot_deltas_fos.append((cw / pw) - (cb / pb))

ci_delta_fos = [float(np.percentile(boot_deltas_fos, 2.5)), float(np.percentile(boot_deltas_fos, 97.5))]
p_emp_fos = float(np.mean(np.array(boot_deltas_fos) <= 0))

print(f"\nContrast 2: Polymerase Biophysical Dissociation on UL54-Foscarnet (Where Proximity vs BLOSUM62 Radicalness)")
print(f"  UL54-FOS Proximity S0: {obs_fos_where:.4f} | UL54-FOS BLOSUM62 S0: {obs_fos_blosum:.4f}")
print(f"  Observed Delta: +{delta_obs_fos:.4f} [95% CI: {ci_delta_fos[0]:.4f}, {ci_delta_fos[1]:.4f}] | Empirical p = {p_emp_fos:.4f}")

results['contrasts']['kinase_vs_polymerase_gcv_blosum'] = {
    'ul97_gcv_s0': round(obs_pw_97, 4),
    'ul54_gcv_s0': round(obs_pw_54, 4),
    'delta_observed': round(delta_obs_gcv, 4),
    'delta_ci_95': [round(x, 4) for x in ci_delta_gcv],
    'empirical_p': round(p_emp_gcv, 4),
    'verdict': 'CONFIRMED (Positive Contrast, Kinase > Polymerase)' if ci_delta_gcv[0] > 0 else 'DESCRIPTIVE'
}

results['contrasts']['polymerase_fos_where_vs_blosum'] = {
    'ul54_fos_proximity_s0': round(obs_fos_where, 4),
    'ul54_fos_blosum_s0': round(obs_fos_blosum, 4),
    'delta_observed': round(delta_obs_fos, 4),
    'delta_ci_95': [round(x, 4) for x in ci_delta_fos],
    'empirical_p': round(p_emp_fos, 4),
    'verdict': 'CONFIRMED (Massive Proximity Advantage over Radicalness)' if ci_delta_fos[0] > 0 else 'NOT SUPPORTED'
}

out_log = os.path.join(ROOT, 'logs', 'step2_hcmv_replication_results.json')
with open(out_log, 'w') as f:
    json.dump(results, f, indent=2)

print("\n" + "="*110)
print(f"Step 2 replication results successfully serialized to {out_log}")
print("="*110)
