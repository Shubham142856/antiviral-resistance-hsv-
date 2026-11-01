import pathlib
import os
import hashlib
import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser
from sklearn.cluster import KMeans
from download_and_audit_structures import get_hsv1_reference_sequences, parse_chain_ca

def main():
    parser = MMCIFParser(QUIET=True)
    ref_seqs = get_hsv1_reference_sequences()
    
    # 1. Load structures
    s_8v1q = parser.get_structure('8V1Q', str(pathlib.Path(__file__).resolve().parent.parent / 'data/structures/8v1q.cif'))
    ca_ul30_list = parse_chain_ca(s_8v1q, 'A')
    
    s_1of1 = parser.get_structure('1OF1', str(pathlib.Path(__file__).resolve().parent.parent / 'data/structures/1of1.cif'))
    ca_ul23_list = parse_chain_ca(s_1of1, 'B')
    
    ca_ul30 = {r['resseq']: (r['x'], r['y'], r['z']) for r in ca_ul30_list}
    ca_ul23 = {r['resseq']: (r['x'], r['y'], r['z']) for r in ca_ul23_list}
    
    # 2. Load benchmark dataset
    df = pd.read_csv(str(pathlib.Path(__file__).resolve().parent.parent / 'data/processed/primary_position_folds.tsv'), sep='\t')
    
    # Check spatial clustering for each gene
    coords_dict = {'UL30': ca_ul30, 'UL23': ca_ul23}
    
    # We will build spatial fold assignments
    df['has_3d_coord'] = 0
    df['coord_x'] = np.nan
    df['coord_y'] = np.nan
    df['coord_z'] = np.nan
    df['spatial_cluster_id'] = -1
    
    for gene in ['UL23', 'UL30']:
        gene_coords = coords_dict[gene]
        # Get all positions in the protein that have coordinates
        resolved_positions = sorted(list(gene_coords.keys()))
        X_resolved = np.array([gene_coords[p] for p in resolved_positions])
        
        # Fit KMeans with k=5, seed=42
        kmeans = KMeans(n_clusters=5, random_state=42, n_init=10)
        labels = kmeans.fit_predict(X_resolved)
        pos_to_cluster = {p: int(lbl) for p, lbl in zip(resolved_positions, labels)}
        
        # For unmodeled positions, find the nearest sequence position that is resolved
        ref_len = len(ref_seqs[gene])
        all_positions = list(range(1, ref_len + 1))
        full_pos_to_cluster = {}
        for p in all_positions:
            if p in pos_to_cluster:
                full_pos_to_cluster[p] = pos_to_cluster[p]
            else:
                # Nearest sequence neighbor with resolved coords
                nearest_p = min(resolved_positions, key=lambda x: abs(x - p))
                full_pos_to_cluster[p] = pos_to_cluster[nearest_p]
                
        # Assign to df
        mask = (df['gene'] == gene)
        for idx in df[mask].index:
            p = df.loc[idx, 'pos']
            if p in gene_coords:
                df.loc[idx, 'has_3d_coord'] = 1
                df.loc[idx, 'coord_x'] = gene_coords[p][0]
                df.loc[idx, 'coord_y'] = gene_coords[p][1]
                df.loc[idx, 'coord_z'] = gene_coords[p][2]
            df.loc[idx, 'spatial_cluster_id'] = full_pos_to_cluster.get(p, -1)
            
    # Inspect cluster distributions per (gene, drug)
    print("=== Spatial Cluster Breakdown (Full Cohort) ===")
    for (g, d), grp in df.groupby(['gene', 'drug']):
        print(f"\n--- {g}-{d} (Total N = {len(grp)}) ---")
        ct = grp.groupby(['spatial_cluster_id', 'target_binary']).size().unstack(fill_value=0)
        print(ct)
        # Check if any cluster has single class
        single_class = [cid for cid in ct.index if 0 in ct.loc[cid] and ct.loc[cid, 0] == 0 or 1 in ct.loc[cid] and ct.loc[cid, 1] == 0]
        if single_class:
            print(f"Warning: Single class clusters in {g}-{d}: {single_class}")
        else:
            print(f"All 5 clusters have BOTH positive and negative classes! Excellent for AUROC calculation.")
            
    print("\n=== Spatial Cluster Breakdown (Structurally Resolved Core) ===")
    df_res = df[df['has_3d_coord'] == 1]
    for (g, d), grp in df_res.groupby(['gene', 'drug']):
        print(f"\n--- {g}-{d} Resolved Core (N = {len(grp)}) ---")
        ct = grp.groupby(['spatial_cluster_id', 'target_binary']).size().unstack(fill_value=0)
        print(ct)
        single_class = [cid for cid in ct.index if 0 in ct.loc[cid] and ct.loc[cid, 0] == 0 or 1 in ct.loc[cid] and ct.loc[cid, 1] == 0]
        if single_class:
            print(f"Warning: Single class clusters in {g}-{d} core: {single_class}")
        else:
            print(f"All 5 clusters have BOTH positive and negative classes in resolved core!")

    # Check 1D sequence block vs 3D spatial cluster alignment
    from compute_conservation_baseline import assign_contiguous_block
    df['seq_block_id'] = df.apply(lambda r: assign_contiguous_block(r['gene'], r['pos']), axis=1)
    
    print("\n=== Cross-Tabulation: 1D Sequence Block vs 3D Spatial Cluster ===")
    for (g, d), grp in df.groupby(['gene', 'drug']):
        print(f"\n--- {g}-{d} (1D Block vs 3D Cluster) ---")
        ct = pd.crosstab(grp['seq_block_id'], grp['spatial_cluster_id'], margins=True)
        print(ct)

if __name__ == '__main__':
    main()

