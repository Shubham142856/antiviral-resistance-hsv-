import pathlib
import os
import sys
import urllib.request
import gzip
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1

# 1. Canonical Reference Sequences from NC_001806.2
def get_hsv1_reference_sequences():
    ref_path = str(pathlib.Path(__file__).resolve().parent.parent / 'herpesdrg/inst/ref/NC_001806.2.fasta')
    ref = SeqIO.read(ref_path, 'fasta')
    # UL23 (TK)
    ul23_nt = ref.seq[46672:47803].reverse_complement()
    ul23_aa = str(ul23_nt.translate())[:-1]
    # UL30 (Pol)
    ul30_nt = ref.seq[62806:66514]
    ul30_aa = str(ul30_nt.translate())[:-1]
    return {'UL23': ul23_aa, 'UL30': ul30_aa}

def download_cif(pdb_id, dest_dir):
    pdb_id = pdb_id.lower()
    cif_path = os.path.join(dest_dir, f"{pdb_id}.cif")
    if os.path.exists(cif_path) and os.path.getsize(cif_path) > 1000:
        print(f"[{pdb_id.upper()}] Already downloaded at {cif_path}")
        return cif_path
    
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    print(f"[{pdb_id.upper()}] Downloading from {url}...")
    urllib.request.urlretrieve(url, cif_path)
    print(f"[{pdb_id.upper()}] Saved to {cif_path} ({os.path.getsize(cif_path):,} bytes)")
    return cif_path

def parse_chain_ca(structure, chain_id):
    chain = structure[0][chain_id]
    residues = []
    for res in chain:
        # Check standard amino acid residue
        res_name = res.get_resname()
        try:
            one_letter = seq1(res_name)
        except Exception:
            continue
        if one_letter == 'X' or not one_letter:
            continue
        # Check for CA atom
        if 'CA' in res:
            ca = res['CA']
            coord = ca.get_coord() # numpy array (3,)
            res_id = res.get_id() # (' ', resseq, icode)
            residues.append({
                'resseq': res_id[1],
                'icode': res_id[2],
                'resname': res_name,
                'aa': one_letter,
                'x': float(coord[0]),
                'y': float(coord[1]),
                'z': float(coord[2])
            })
    return residues

def align_and_map_coords(ref_seq, chain_residues):
    """
    Map chain residues to 1-based reference sequence positions using exact subsequence matching or alignment.
    """
    struct_seq = "".join([r['aa'] for r in chain_residues])
    print(f"  Structure sequence length (CA atoms): {len(struct_seq)}")
    
    # Check if resseq directly matches reference coordinates
    direct_match = 0
    direct_mismatch = 0
    mapping = {}
    
    for r in chain_residues:
        pos = r['resseq']
        if 1 <= pos <= len(ref_seq):
            ref_aa = ref_seq[pos - 1]
            if ref_aa == r['aa']:
                direct_match += 1
                mapping[pos] = (r['x'], r['y'], r['z'], r['aa'])
            else:
                direct_mismatch += 1
                
    if direct_match > len(chain_residues) * 0.8:
        print(f"  Direct numbering matches reference! Match: {direct_match}, Mismatch: {direct_mismatch}")
        return mapping
        
    # Otherwise perform local/global alignment
    from Bio import pairwise2
    alignments = pairwise2.align.globalms(ref_seq, struct_seq, 2, -1, -5, -0.5)
    best = alignments[0]
    ref_aln, struct_aln = best[0], best[1]
    
    mapping = {}
    ref_idx = 0
    struct_idx = 0
    matched = 0
    for a_ref, a_str in zip(ref_aln, struct_aln):
        curr_ref = ref_idx + 1 if a_ref != '-' else None
        curr_str = chain_residues[struct_idx] if a_str != '-' else None
        
        if a_ref != '-' and a_str != '-':
            if a_ref == a_str:
                matched += 1
                mapping[curr_ref] = (curr_str['x'], curr_str['y'], curr_str['z'], curr_str['aa'])
            else:
                # Still map if high confidence alignment
                mapping[curr_ref] = (curr_str['x'], curr_str['y'], curr_str['z'], curr_str['aa'])
                
        if a_ref != '-':
            ref_idx += 1
        if a_str != '-':
            struct_idx += 1
            
    print(f"  Aligned to reference! Aligned matches: {matched}/{len(struct_seq)} ({matched/len(struct_seq)*100:.1f}%)")
    return mapping

def main():
    dest_dir = str(pathlib.Path(__file__).resolve().parent.parent / 'data/structures')
    ref_seqs = get_hsv1_reference_sequences()
    print(f"Reference UL23 length: {len(ref_seqs['UL23'])} aa")
    print(f"Reference UL30 length: {len(ref_seqs['UL30'])} aa")
    
    # Load benchmark dataset to see position coverage
    folds_path = str(pathlib.Path(__file__).resolve().parent.parent / 'data/processed/primary_position_folds.tsv')
    df_folds = pd.read_csv(folds_path, sep='\t')
    
    target_pdbs = {
        'UL30': ['8V1Q', '8EXX', '8V1T'],
        'UL23': ['1OF1', '2C09', '1VTK']
    }
    
    parser = MMCIFParser(QUIET=True)
    
    coverage_results = []
    
    for gene, pdbs in target_pdbs.items():
        ref_seq = ref_seqs[gene]
        gene_df = df_folds[df_folds['gene'] == gene]
        gene_unique_pos = set(gene_df['pos'].unique())
        print(f"\n==========================================")
        print(f"Auditing structures for {gene} (Unique benchmark positions: {len(gene_unique_pos)})")
        print(f"==========================================")
        
        for pdb_id in pdbs:
            try:
                cif_file = download_cif(pdb_id, dest_dir)
                structure = parser.get_structure(pdb_id, cif_file)
                
                # Find the best chain
                best_chain = None
                best_len = 0
                for chain in structure[0]:
                    ca_count = sum(1 for res in chain if 'CA' in res)
                    if ca_count > best_len:
                        best_len = ca_count
                        best_chain = chain.id
                        
                print(f"\n[{pdb_id}] Best chain: {best_chain} (CA atoms: {best_len})")
                chain_res = parse_chain_ca(structure, best_chain)
                pos_mapping = align_and_map_coords(ref_seq, chain_res)
                
                # Check benchmark position coverage
                covered_benchmark_pos = gene_unique_pos.intersection(set(pos_mapping.keys()))
                missing_benchmark_pos = gene_unique_pos - set(pos_mapping.keys())
                
                # Check row coverage
                covered_rows = gene_df[gene_df['pos'].isin(pos_mapping.keys())]
                row_coverage_pct = len(covered_rows) / len(gene_df) * 100
                pos_coverage_pct = len(covered_benchmark_pos) / len(gene_unique_pos) * 100
                
                print(f"[{pdb_id}] Total reference positions resolved: {len(pos_mapping)}/{len(ref_seq)} ({len(pos_mapping)/len(ref_seq)*100:.1f}%)")
                print(f"[{pdb_id}] Benchmark unique position coverage: {len(covered_benchmark_pos)}/{len(gene_unique_pos)} ({pos_coverage_pct:.1f}%)")
                print(f"[{pdb_id}] Benchmark row coverage: {len(covered_rows)}/{len(gene_df)} ({row_coverage_pct:.1f}%)")
                if missing_benchmark_pos:
                    print(f"[{pdb_id}] Missing benchmark positions ({len(missing_benchmark_pos)}): {sorted(list(missing_benchmark_pos))}")
                    
                coverage_results.append({
                    'gene': gene,
                    'pdb_id': pdb_id,
                    'chain': best_chain,
                    'structure_len': len(chain_res),
                    'ref_len': len(ref_seq),
                    'ref_covered': len(pos_mapping),
                    'benchmark_unique_pos': len(gene_unique_pos),
                    'benchmark_covered_pos': len(covered_benchmark_pos),
                    'benchmark_missing_pos_count': len(missing_benchmark_pos),
                    'benchmark_pos_pct': pos_coverage_pct,
                    'benchmark_rows_total': len(gene_df),
                    'benchmark_rows_covered': len(covered_rows),
                    'benchmark_rows_pct': row_coverage_pct,
                    'missing_positions': str(sorted(list(missing_benchmark_pos)))
                })
            except Exception as e:
                print(f"Error processing {pdb_id}: {e}")
                
    cov_df = pd.DataFrame(coverage_results)
    out_csv = str(pathlib.Path(__file__).resolve().parent.parent / 'logs/structure_coverage_audit.csv')
    cov_df.to_csv(out_csv, index=False)
    print(f"\nSaved coverage audit results to: {out_csv}")
    print("\nSummary Table:")
    print(cov_df[['gene', 'pdb_id', 'chain', 'ref_covered', 'benchmark_covered_pos', 'benchmark_pos_pct', 'benchmark_rows_pct']].to_string(index=False))

if __name__ == '__main__':
    main()

