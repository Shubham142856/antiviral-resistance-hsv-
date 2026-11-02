"""
build_kaggle_kernel.py
Prepares a self-contained Kaggle kernel script to evaluate ESM-2 650M and 35M
positional constraint on HCMV UL54 and compute the pre-registered composite Where score.
"""

import os
import json
import re
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, MMCIFParser
from Bio.Align import PairwiseAligner, substitution_matrices

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
import sys
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

# 1. UL54 Sequence
from scripts.map_active_site_to_hcmv import get_gene_seq
seq_ul54 = get_gene_seq('NC_006273.2', 'UL54')
seq_ul30 = get_gene_seq('NC_001806.2', 'UL30')

# 2. Option A distances
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

# 3. Option B distances
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

aligner = PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
aligner.open_gap_score = -10.0
aligner.extend_gap_score = -0.5
aligner.mode = 'global'
aln = aligner.align(seq_ul30, seq_ul54)[0]

hsv_idx, hcmv_idx = 0, 0
hcmv_to_hsv = {}
for r_ch, q_ch in zip(aln[0], aln[1]):
    if r_ch != '-': hsv_idx += 1
    if q_ch != '-': hcmv_idx += 1
    if r_ch != '-' and q_ch != '-': hcmv_to_hsv[hcmv_idx] = hsv_idx

dists_option_b = {p: hsv_dists.get(hcmv_to_hsv.get(p), np.nan) for p in range(1, len(seq_ul54)+1)}
max_dist_b = np.nanmax(list(dists_option_b.values()))

# Clean up nan in dists_option_b for serialization
dists_option_b_clean = {k: (v if not np.isnan(v) else None) for k, v in dists_option_b.items()}

# 4. HerpesDRG UL54 rows
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

ul54_df = hcmv_df[hcmv_df['gene'] == 'UL54'].copy()

rows_data = []
for _, r in ul54_df.iterrows():
    rows_data.append({
        'aa_change': str(r['aa_change']),
        'pos': int(r['pos']),
        'wt': str(r['wt']),
        'mut': str(r['mut']),
        'study_id': str(r['study_id']),
        'Ganciclovir': str(r['Ganciclovir']),
        'Foscarnet': str(r['Foscarnet']),
        'Cidofovir': str(r['Cidofovir'])
    })

print(f"Serialized UL54 data: {len(rows_data)} rows, {len(dists_option_a)} Option A dists, {len(dists_option_b)} Option B dists")

# Ensure cloud_gpu/hcmv_composite_kernel exists
out_dir = os.path.join(ROOT, 'cloud_gpu', 'hcmv_composite_kernel')
os.makedirs(out_dir, exist_ok=True)

# Write metadata
metadata = {
    "id": "shubhamraikwar142/hsv-plm-gpu-benchmark",
    "title": "HSV PLM GPU Benchmark",
    "code_file": "run_hcmv_composite.py",
    "language": "python",
    "kernel_type": "script",
    "is_private": "true",
    "enable_gpu": "true",
    "enable_tpu": "false",
    "enable_internet": "true",
    "dataset_sources": [],
    "competition_sources": [],
    "kernel_sources": [],
    "model_sources": []
}

with open(os.path.join(out_dir, 'kernel-metadata.json'), 'w') as f:
    json.dump(metadata, f, indent=2)

# Write runner python script
runner_path = os.path.join(out_dir, 'run_hcmv_composite.py')

with open(runner_path, 'w', encoding='utf-8') as f:
    f.write('''"""
run_hcmv_composite.py
Standalone Cloud GPU runner for HCMV UL54 Step 2 Replication:
Evaluates ESM-2 650M (and 35M) masked marginals on UL54, computes the
pre-registered composite Where score (Z(Constraint) - Z(d_site)),
and runs S0 within-study pair-weighted AUROC with B=2000 clustered bootstrap.
"""

import os
import sys
import json
import time
import numpy as np
import pandas as pd
import torch
from sklearn.metrics import roc_auc_score
from transformers import AutoTokenizer, AutoModelForMaskedLM

print("="*80)
print("HCMV UL54 PRE-REGISTERED COMPOSITE BENCHMARK (KAGGLE GPU)")
print("="*80)
print(f"CUDA Available: {torch.cuda.is_available()}")
if torch.cuda.is_available():
    print(f"Device: {torch.cuda.get_device_name(0)}")

# 1. UL54 Sequence & Constants
SEQ_UL54 = ''' + repr(seq_ul54) + '''
DISTS_OPTION_A = ''' + repr(dists_option_a) + '''
MAX_DIST_A = ''' + repr(max_dist_a) + '''
DISTS_OPTION_B = ''' + repr(dists_option_b_clean) + '''
MAX_DIST_B = ''' + repr(max_dist_b) + '''
UL54_ROWS = ''' + repr(rows_data) + '''

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

df_all = pd.DataFrame(UL54_ROWS)
unique_positions = sorted(df_all['pos'].unique())
print(f"Total rows: {len(df_all)}, Unique mutated positions in UL54: {len(unique_positions)}")

# Helper to compute ESM-2 scores for a given model
def compute_esm2_for_positions(model_id, device='cuda' if torch.cuda.is_available() else 'cpu'):
    print(f"\\nLoading {model_id} on {device}...")
    t0 = time.time()
    tokenizer = AutoTokenizer.from_pretrained(model_id)
    model = AutoModelForMaskedLM.from_pretrained(model_id).to(device)
    model.eval()
    print(f"Model loaded in {time.time()-t0:.1f}s")
    
    canonical_aas = 'ACDEFGHIKLMNPQRSTVWY'
    aa_tokens = {a: tokenizer.convert_tokens_to_ids(a) for a in canonical_aas}
    mask_token_id = tokenizer.mask_token_id
    
    full_seq = SEQ_UL54
    pos_cache = {}
    mu_dict = {}
    
    print(f"Computing masked marginals across {len(unique_positions)} positions...")
    t_start = time.time()
    with torch.no_grad():
        for i, pos in enumerate(unique_positions):
            p_idx = pos - 1 # 0-based
            if p_idx < 0 or p_idx >= len(full_seq):
                continue
            wt_aa = full_seq[p_idx]
            
            # Centered 1020-residue window (ESM max length is 1024)
            if len(full_seq) > 1020:
                p_start = max(0, p_idx - 510)
                p_end = min(len(full_seq), p_start + 1020)
                if p_end - p_start < 1020 and p_start > 0:
                    p_start = max(0, p_end - 1020)
                sub_seq = full_seq[p_start:p_end]
                local_p_idx = p_idx - p_start
            else:
                sub_seq = full_seq
                local_p_idx = p_idx
                
            masked_seq = sub_seq[:local_p_idx] + "<mask>" + sub_seq[local_p_idx+1:]
            inputs = tokenizer(masked_seq, return_tensors='pt').to(device)
            
            input_ids = inputs['input_ids'][0]
            mask_indices = (input_ids == mask_token_id).nonzero(as_tuple=True)[0]
            if len(mask_indices) == 0:
                continue
            mask_token_pos = mask_indices[0].item()
            
            logits = model(**inputs).logits[0, mask_token_pos]
            log_probs = torch.log_softmax(logits, dim=-1)
            
            wt_token_id = aa_tokens.get(wt_aa)
            wt_lp = log_probs[wt_token_id].item() if wt_token_id else 0.0
            
            scores = {a: log_probs[aa_tokens[a]].item() - wt_lp for a in canonical_aas if a in aa_tokens}
            pos_cache[pos] = scores
            
            # mu_plm_19: mean of non-wildtype delta LLR
            non_wt = [v for k, v in scores.items() if k != wt_aa]
            mu_dict[pos] = float(np.mean(non_wt)) if len(non_wt) > 0 else 0.0
            
            if (i + 1) % 50 == 0 or (i + 1) == len(unique_positions):
                print(f"  Processed {i+1}/{len(unique_positions)} positions in {time.time()-t_start:.1f}s...")
                
    return pos_cache, mu_dict

# Compute with ESM-2 650M
cache_650m, mu_650m = compute_esm2_for_positions("facebook/esm2_t33_650M_UR50D")

# Compute with ESM-2 35M as baseline
cache_35m, mu_35m = compute_esm2_for_positions("facebook/esm2_t12_35M_UR50D")

# Save caches to disk
with open("esm2_650m_ul54_cache.json", "w") as f:
    json.dump({str(k): v for k, v in cache_650m.items()}, f)
with open("esm2_35m_ul54_cache.json", "w") as f:
    json.dump({str(k): v for k, v in cache_35m.items()}, f)
print("Positional caches saved to disk.")

# Evaluation functions
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

# Shared Bootstrap Engine
RNG = np.random.RandomState(42)
N_BOOT = 2000

drugs = ['Ganciclovir', 'Foscarnet', 'Cidofovir']
results = {
    'protocol': {
        'seed': 42,
        'n_boot': N_BOOT,
        'models': ['ESM-2 650M (facebook/esm2_t33_650M_UR50D)', 'ESM-2 35M (facebook/esm2_t12_35M_UR50D)']
    },
    'targets': {}
}

for drug in drugs:
    t_key = f"HCMV_UL54_{drug}"
    sub = df_all.copy()
    sub['y'] = [is_res(x) for x in sub[drug]]
    sub = sub[sub['y'].notna()].copy()
    sub['y'] = sub['y'].astype(int)
    sub['is_num_fc'] = [is_numeric_fc(x) for x in sub[drug]]
    
    # Deduplicate per (aa_change, study_id) exactly as Step 2
    sub = sub.drop_duplicates(subset=['aa_change', 'study_id']).copy().reset_index(drop=True)
    
    # Features
    sub['blosum_rad'] = [- BLOSUM62.get(r['wt'], {}).get(r['mut'], 0) for _, r in sub.iterrows()]
    
    # Distance Option A & B
    def get_a(pos):
        v = DISTS_OPTION_A.get(pos)
        return float(v) if v is not None and not np.isnan(v) else float(MAX_DIST_A)
    def get_b(pos):
        v = DISTS_OPTION_B.get(pos)
        return float(v) if v is not None and not np.isnan(v) else float(MAX_DIST_B)
        
    sub['dist_a'] = [get_a(r['pos']) for _, r in sub.iterrows()]
    sub['dist_b'] = [get_b(r['pos']) for _, r in sub.iterrows()]
    sub['dist_only_a'] = - sub['dist_a'] # As executed in Step 2: closer = higher
    sub['dist_only_b'] = - sub['dist_b']
    
    # Positional constraints mu_plm_19
    sub['mu_650m'] = [mu_650m.get(r['pos'], 0.0) for _, r in sub.iterrows()]
    sub['mu_35m'] = [mu_35m.get(r['pos'], 0.0) for _, r in sub.iterrows()]
    
    # Spec-Exact Z-scoring for UL54 cohort
    # z_mu = - (mu - mean) / std (higher = more constrained / less tolerant)
    sub['z_mu_650m'] = - (sub['mu_650m'] - sub['mu_650m'].mean()) / (sub['mu_650m'].std() if sub['mu_650m'].std() > 0 else 1.0)
    sub['z_mu_35m'] = - (sub['mu_35m'] - sub['mu_35m'].mean()) / (sub['mu_35m'].std() if sub['mu_35m'].std() > 0 else 1.0)
    
    # z_dist = (dist - mean) / std (closer = lower z_dist)
    sub['z_dist_a'] = (sub['dist_a'] - sub['dist_a'].mean()) / (sub['dist_a'].std() if sub['dist_a'].std() > 0 else 1.0)
    sub['z_dist_b'] = (sub['dist_b'] - sub['dist_b'].mean()) / (sub['dist_b'].std() if sub['dist_b'].std() > 0 else 1.0)
    
    # Composite Where: (z_mu - z_dist) / std
    sub['where_comp_650m_a'] = (sub['z_mu_650m'] - sub['z_dist_a']) / np.std(sub['z_mu_650m'] - sub['z_dist_a'])
    sub['where_comp_35m_a'] = (sub['z_mu_35m'] - sub['z_dist_a']) / np.std(sub['z_mu_35m'] - sub['z_dist_a'])
    sub['where_comp_650m_b'] = (sub['z_mu_650m'] - sub['z_dist_b']) / np.std(sub['z_mu_650m'] - sub['z_dist_b'])
    sub['constraint_only_650m'] = sub['z_mu_650m']
    
    feature_cols = {
        'Distance-Only (Option A: As Executed)': 'dist_only_a',
        'Composite Where (ESM-2 650M + Option A: As Pre-Registered)': 'where_comp_650m_a',
        'Composite Where (ESM-2 35M + Option A)': 'where_comp_35m_a',
        'Constraint-Only (ESM-2 650M Positional)': 'constraint_only_650m',
        'Distance-Only (Option B)': 'dist_only_b',
        'Composite Where (ESM-2 650M + Option B)': 'where_comp_650m_b',
        'BLOSUM62 Radicalness': 'blosum_rad'
    }
    
    tgt_res = {}
    
    # Compute stats for each feature
    unique_studies = sub['study_id'].unique()
    
    for f_label, f_col in feature_cols.items():
        base_stats = get_study_pair_stats(sub, f_col)
        pw_auc, n_p = compute_pair_weighted_auroc(base_stats)
        sm_auc = compute_study_mean_auroc(base_stats)
        
        # LOSO
        loso_vals = []
        for s_ex in base_stats.keys():
            loso_stats = {k: v for k, v in base_stats.items() if k != s_ex}
            l_pw, _ = compute_pair_weighted_auroc(loso_stats)
            if not np.isnan(l_pw): loso_vals.append(l_pw)
        loso_range = [round(float(min(loso_vals)), 4), round(float(max(loso_vals)), 4)] if loso_vals else [None, None]
        
        # S2
        sub_s2 = sub[sub['is_num_fc']].copy()
        s2_stats = get_study_pair_stats(sub_s2, f_col)
        s2_pw, s2_p = compute_pair_weighted_auroc(s2_stats)
        s2_pooled = float(roc_auc_score(sub_s2['y'], sub_s2[f_col])) if sub_s2['y'].nunique() > 1 else np.nan
        
        # D1 pooled
        d1_auc = float(roc_auc_score(sub['y'], sub[f_col])) if sub['y'].nunique() > 1 else np.nan
        
        tgt_res[f_label] = {
            'col': f_col,
            's0_pair_weighted': pw_auc,
            's0_study_mean': sm_auc,
            'loso_range': loso_range,
            's2_pair_weighted': s2_pw,
            's2_pairs': s2_p,
            's2_pooled': s2_pooled,
            'd1_pooled': d1_auc,
            'boot_pw': [],
            'boot_sm': [],
            'boot_d1': []
        }
        
    # Run joint bootstrap across studies
    study_dfs = {s_id: sub[sub['study_id'] == s_id] for s_id in unique_studies}
    
    for b in range(N_BOOT):
        sampled_studies = RNG.choice(unique_studies, size=len(unique_studies), replace=True)
        boot_df = pd.concat([study_dfs[s] for s in sampled_studies], ignore_index=True)
        
        for f_label, f_col in feature_cols.items():
            b_stats = get_study_pair_stats(boot_df, f_col)
            b_pw, _ = compute_pair_weighted_auroc(b_stats)
            b_sm = compute_study_mean_auroc(b_stats)
            b_d1 = float(roc_auc_score(boot_df['y'], boot_df[f_col])) if boot_df['y'].nunique() > 1 else np.nan
            
            if not np.isnan(b_pw): tgt_res[f_label]['boot_pw'].append(b_pw)
            if not np.isnan(b_sm): tgt_res[f_label]['boot_sm'].append(b_sm)
            if not np.isnan(b_d1): tgt_res[f_label]['boot_d1'].append(b_d1)
            
    # Finalize CIs
    final_features = {}
    for f_label, f_data in tgt_res.items():
        pw_ci = [round(float(np.percentile(f_data['boot_pw'], 2.5)), 4), round(float(np.percentile(f_data['boot_pw'], 97.5)), 4)] if f_data['boot_pw'] else [None, None]
        sm_ci = [round(float(np.percentile(f_data['boot_sm'], 2.5)), 4), round(float(np.percentile(f_data['boot_sm'], 97.5)), 4)] if f_data['boot_sm'] else [None, None]
        d1_ci = [round(float(np.percentile(f_data['boot_d1'], 2.5)), 4), round(float(np.percentile(f_data['boot_d1'], 97.5)), 4)] if f_data['boot_d1'] else [None, None]
        
        final_features[f_label] = {
            's0_pair_weighted_auroc': round(float(f_data['s0_pair_weighted']), 4),
            's0_pair_weighted_ci_95': pw_ci,
            's0_study_mean_auroc': round(float(f_data['s0_study_mean']), 4),
            's0_study_mean_ci_95': sm_ci,
            's0_loso_range': f_data['loso_range'],
            's2_numeric_fc_pair_weighted_auroc': round(float(f_data['s2_pair_weighted']), 4) if not np.isnan(f_data['s2_pair_weighted']) else None,
            's2_numeric_fc_pairs': f_data['s2_pairs'],
            's2_numeric_fc_pooled_auroc': round(float(f_data['s2_pooled']), 4) if not np.isnan(f_data['s2_pooled']) else None,
            'd1_pooled_auroc': round(float(f_data['d1_pooled']), 4) if not np.isnan(f_data['d1_pooled']) else None,
            'd1_pooled_ci_95': d1_ci
        }
        
    results['targets'][t_key] = {
        'drug': drug,
        'features': final_features
    }

# Save results to JSON
out_json_path = "step2_hcmv_composite_results.json"
with open(out_json_path, "w") as f:
    json.dump(results, f, indent=2)
print(f"\\nSaved replication benchmark results to {out_json_path}")

# Print Comparison Table
print("\\n" + "="*110)
print("SIDE-BY-SIDE REPLICATION RESULTS: AS EXECUTED (DISTANCE-ONLY) vs AS PRE-REGISTERED (COMPOSITE)")
print("="*110)
print(f"{'Target Arm':<20} | {'Score Definition':<55} | {'S0 Pair-Weighted [95% CI]':<30} | {'Study-Mean':<12}")
print("-" * 125)
for drug in drugs:
    t_k = f"HCMV_UL54_{drug}"
    feats = results['targets'][t_k]['features']
    for f_lbl in [
        'Distance-Only (Option A: As Executed)',
        'Composite Where (ESM-2 650M + Option A: As Pre-Registered)',
        'Composite Where (ESM-2 35M + Option A)',
        'Constraint-Only (ESM-2 650M Positional)',
        'BLOSUM62 Radicalness'
    ]:
        v = feats[f_lbl]
        pw_str = f"{v['s0_pair_weighted_auroc']:.4f} [{v['s0_pair_weighted_ci_95'][0]:.4f}, {v['s0_pair_weighted_ci_95'][1]:.4f}]"
        sm_str = f"{v['s0_study_mean_auroc']:.4f}"
        print(f"{f'UL54-{drug}':<20} | {f_lbl:<55} | {pw_str:<30} | {sm_str:<12}")
    print("-" * 125)

print("\\nBenchmark execution complete!")
''')

print(f"Generated {runner_path} successfully!")
