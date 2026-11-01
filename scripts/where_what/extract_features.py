"""
extract_features.py
Extracts orthogonal Where (positional evolutionary constraint + catalytic site distance)
and What (within-residue centered severity + chemical radicalness) features.
Anchored to pre-registration: docs/WHERE_WHAT_PREREGISTRATION.md
"""

import json
import math
import numpy as np
import pandas as pd
from Bio.PDB.MMCIFParser import MMCIFParser

# 1. Load primary benchmark cohort
df = pd.read_csv('data/processed/primary_position_folds.tsv', sep='\t')
print(f"Loaded primary benchmark cohort: {df.shape[0]} rows across {df['gene'].nunique()} genes.")

# 2. Load ESM-2 position cache (all 20 amino acid scores per position)
with open('data/processed/esm2_pos_cache.json') as f:
    esm2_cache = json.load(f)

# BLOSUM62 matrix
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

# 3. Parse 3D Structures and extract active site / drug distances
print("Parsing experimental structures for active site coordinates...")
parser = MMCIFParser(QUIET=True)

# UL30 active site coordinates from 8EXX (Chain A + PPF ligand / MG ions / catalytic Asp717, Asp888, Glu889)
s_8exx = parser.get_structure('8EXX', 'data/structures/8exx.cif')
m_8exx = list(s_8exx.get_models())[0]

# Collect target site atoms in 8EXX
ul30_site_coords = []
for res in m_8exx.get_residues():
    resname = res.get_resname()
    # PPF ligand (foscarnet) and MG catalytic ions
    if resname in ['PPF', 'MG']:
        for a in res.get_atoms():
            ul30_site_coords.append(a.get_coord())
    # Catalytic triad: Asp717, Asp888, Glu889 in Chain A
    if res.get_parent().id == 'A' and res.id[1] in [717, 888, 889]:
        for a in res.get_atoms():
            ul30_site_coords.append(a.get_coord())

ul30_site_coords = np.array(ul30_site_coords)
print(f"UL30 (8EXX) active-site reference atoms: {len(ul30_site_coords)}")

# Map UL30 Chain A CA coordinates
ul30_ca_coords = {}
for res in m_8exx['A'].get_residues():
    if 'CA' in res:
        ul30_ca_coords[res.id[1]] = res['CA'].get_coord()

# UL23 active site coordinates from 1OF1 (Chain B + SCT ligand / Arg163, Glu225)
s_1of1 = parser.get_structure('1OF1', 'data/structures/1of1.cif')
m_1of1 = list(s_1of1.get_models())[0]

ul23_site_coords = []
for res in m_1of1.get_residues():
    if res.get_resname() == 'SCT':
        for a in res.get_atoms():
            ul23_site_coords.append(a.get_coord())
    if res.get_parent().id == 'B' and res.id[1] in [163, 225]:
        for a in res.get_atoms():
            ul23_site_coords.append(a.get_coord())

ul23_site_coords = np.array(ul23_site_coords)
print(f"UL23 (1OF1) active-site reference atoms: {len(ul23_site_coords)}")

ul23_ca_coords = {}
for res in m_1of1['B'].get_residues():
    if 'CA' in res:
        ul23_ca_coords[res.id[1]] = res['CA'].get_coord()

# Distance helper
def get_site_distance(gene, pos):
    if gene == 'UL30':
        if pos in ul30_ca_coords and len(ul30_site_coords) > 0:
            ca = ul30_ca_coords[pos]
            dists = np.linalg.norm(ul30_site_coords - ca, axis=1)
            return float(np.min(dists))
    elif gene == 'UL23':
        if pos in ul23_ca_coords and len(ul23_site_coords) > 0:
            ca = ul23_ca_coords[pos]
            dists = np.linalg.norm(ul23_site_coords - ca, axis=1)
            return float(np.min(dists))
    return np.nan

# Domain mapping helpers
def assign_ul30_domains(pos):
    if pos <= 362:
        return 'NTD', 'NTD (1-362)'
    elif pos <= 593:
        return 'Exo', 'Exo (363-593)'
    elif pos <= 766:
        return 'Palm', 'Palm_1 (594-766)'
    elif pos <= 825:
        return 'Fingers', 'Fingers (767-825)'
    elif pos <= 956:
        return 'Palm', 'Palm_2 (826-956)'
    else:
        return 'Thumb_Cterm', 'Thumb_Cterm (957-1235)'

def assign_ul23_domains(pos):
    if pos < 46:
        return 'N-term', 'N-term (1-45)'
    elif pos <= 160:
        return 'ATP_Core', 'ATP_Core (46-160)'
    elif pos <= 240:
        return 'Nucleoside_Pocket', 'Nucleoside_Pocket (161-240)'
    elif pos <= 320:
        return 'Lid_Dimer', 'Lid_Dimer (241-320)'
    else:
        return 'C-term', 'C-term (321-376)'

# Compute features for each variant
mu_plm_list = []
sigma_plm_list = []
sev_plm_list = []
sev_blosum_list = []
site_dist_list = []
dom5_list = []
dom_detail_list = []

for idx, r in df.iterrows():
    gene = r['gene']
    pos = int(r['pos'])
    wt = r['wt']
    mut = r['mut']
    
    # Domain
    if gene == 'UL30':
        d5, dd = assign_ul30_domains(pos)
    else:
        d5, dd = assign_ul23_domains(pos)
    dom5_list.append(d5)
    dom_detail_list.append(dd)
    
    # ESM-2 19-substitution positional constraint
    k = f"{gene}:{pos}"
    if k in esm2_cache:
        scores = esm2_cache[k]
        # non-wildtype values
        non_wt_vals = [scores[aa] for aa in AA_LIST if aa != wt and aa in scores]
        mu = float(np.mean(non_wt_vals)) if len(non_wt_vals) > 0 else 0.0
        sig = float(np.std(non_wt_vals)) if len(non_wt_vals) > 0 else 1.0
        mut_val = float(scores.get(mut, mu))
        sev = mut_val - mu
    else:
        mu, sig, sev = np.nan, np.nan, np.nan
        
    mu_plm_list.append(mu)
    sigma_plm_list.append(sig)
    sev_plm_list.append(sev)
    
    # BLOSUM62 radicalness (negative score: higher = more radical)
    b_val = -BLOSUM62.get(wt, {}).get(mut, 0)
    sev_blosum_list.append(b_val)
    
    # Active site 3D distance
    dist = get_site_distance(gene, pos)
    site_dist_list.append(dist)

df['domain_5'] = dom5_list
df['domain_detail'] = dom_detail_list
df['mu_plm_19'] = mu_plm_list
df['sigma_plm_19'] = sigma_plm_list
df['severity_plm'] = sev_plm_list
df['severity_blosum'] = sev_blosum_list
df['dist_to_site_angstrom'] = site_dist_list

# Impute missing distances with gene median for unresolved positions
for g in ['UL23', 'UL30']:
    mask = (df['gene'] == g)
    med_dist = df.loc[mask, 'dist_to_site_angstrom'].median()
    df.loc[mask & df['dist_to_site_angstrom'].isna(), 'dist_to_site_angstrom'] = med_dist

# Standardize Where components per gene
for g in ['UL23', 'UL30']:
    mask = (df['gene'] == g)
    # Z(mu_plm)
    mu_mean = df.loc[mask, 'mu_plm_19'].mean()
    mu_std = df.loc[mask, 'mu_plm_19'].std()
    z_mu = (df.loc[mask, 'mu_plm_19'] - mu_mean) / mu_std
    
    # Z(dist)
    d_mean = df.loc[mask, 'dist_to_site_angstrom'].mean()
    d_std = df.loc[mask, 'dist_to_site_angstrom'].std()
    z_dist = (df.loc[mask, 'dist_to_site_angstrom'] - d_mean) / d_std
    
    # Composite Where: high constraint (high mu_plm) AND close to site (low dist)
    df.loc[mask, 'where_composite'] = z_mu - z_dist
    
    # Standardize severities per gene
    df.loc[mask, 'z_severity_plm'] = (df.loc[mask, 'severity_plm'] - df.loc[mask, 'severity_plm'].mean()) / df.loc[mask, 'severity_plm'].std()
    df.loc[mask, 'z_severity_blosum'] = (df.loc[mask, 'severity_blosum'] - df.loc[mask, 'severity_blosum'].mean()) / df.loc[mask, 'severity_blosum'].std()

out_path = 'data/processed/where_what_features.tsv'
df.to_csv(out_path, sep='\t', index=False)
print(f"Saved Where x What features to {out_path} ({df.shape[0]} rows, {df.shape[1]} columns).")
print("\nFeature Summary by Target Pair:")
for p in [('UL23', 'Aciclovir'), ('UL30', 'Aciclovir'), ('UL30', 'Foscarnet')]:
    sub = df[(df['gene'] == p[0]) & (df['drug'] == p[1])]
    print(f"{p[0]}-{p[1]}: N={len(sub)}, Pos={sub['target_binary'].sum():.0f} | Where: [{sub['where_composite'].min():.2f}, {sub['where_composite'].max():.2f}] | SevPLM: [{sub['z_severity_plm'].min():.2f}, {sub['z_severity_plm'].max():.2f}]")
