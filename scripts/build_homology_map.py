import pathlib
"""
Build Pairwise Sequence Alignment and Homology Map for HSV-1 vs HSV-2:
- UL23 (Thymidine Kinase: HSV-1 376 aa vs HSV-2 376 aa)
- UL30 (DNA Polymerase: HSV-1 1235 aa vs HSV-2 1240 aa)
Produces canonical homologous cluster IDs so homologous positions share identical cross-validation folds.
Saves to data/processed/hsv1_hsv2_homology_map.tsv and logs/homology_mapping.log.
"""
from Bio import SeqIO
from Bio.Align import PairwiseAligner
import pandas as pd

def main():
    ref1_path = str(pathlib.Path(__file__).resolve().parent.parent / 'herpesdrg/inst/ref/NC_001806.2.fasta')
    ref2_path = str(pathlib.Path(__file__).resolve().parent.parent / 'herpesdrg/inst/ref/NC_001798.2.fasta')
    
    ref1 = SeqIO.read(ref1_path, 'fasta')
    ref2 = SeqIO.read(ref2_path, 'fasta')
    
    # UL23 coordinates
    ul23_1_nt = ref1.seq[46672:47803].reverse_complement()
    ul23_2_nt = ref2.seq[46908:48039].reverse_complement()
    ul23_1_aa = str(ul23_1_nt.translate())[:-1]
    ul23_2_aa = str(ul23_2_nt.translate())[:-1]
    
    # UL30 coordinates
    ul30_1_nt = ref1.seq[62806:66514]
    ul30_2_nt = ref2.seq[63299:67022]
    ul30_1_aa = str(ul30_1_nt.translate())[:-1]
    ul30_2_aa = str(ul30_2_nt.translate())[:-1]
    
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    
    records = []
    log_lines = ["=== HSV-1 vs HSV-2 PAIRWISE HOMOLOGY MAPPING LOG ===", ""]
    
    for gene, s1, s2 in [('UL23', ul23_1_aa, ul23_2_aa), ('UL30', ul30_1_aa, ul30_2_aa)]:
        aln = aligner.align(s1, s2)[0]
        aligned_s1, aligned_s2 = aln[0], aln[1]
        
        log_lines.append(f"Gene: {gene}")
        log_lines.append(f"  HSV-1 Length: {len(s1)} aa | HSV-2 Length: {len(s2)} aa")
        log_lines.append(f"  Alignment Length: {len(aligned_s1)} columns | Score: {aln.score}")
        
        idx1 = 0
        idx2 = 0
        cluster_counter = 1
        
        for col_idx in range(len(aligned_s1)):
            c1 = aligned_s1[col_idx]
            c2 = aligned_s2[col_idx]
            
            p1 = None
            p2 = None
            
            if c1 != '-':
                idx1 += 1
                p1 = idx1
            if c2 != '-':
                idx2 += 1
                p2 = idx2
                
            cluster_id = f"{gene}_hom_pos_{cluster_counter}"
            cluster_counter += 1
            
            records.append({
                'gene': gene,
                'alignment_col': col_idx + 1,
                'homology_cluster_id': cluster_id,
                'hsv1_pos': p1,
                'hsv1_aa': c1 if c1 != '-' else None,
                'hsv2_pos': p2,
                'hsv2_aa': c2 if c2 != '-' else None,
                'is_identical': 1 if c1 == c2 and c1 != '-' else 0,
                'is_gap': 1 if c1 == '-' or c2 == '-' else 0
            })
            
        # Summary stats
        gene_recs = [r for r in records if r['gene'] == gene]
        gaps = sum(r['is_gap'] for r in gene_recs)
        ident = sum(r['is_identical'] for r in gene_recs)
        log_lines.append(f"  Identical Residues: {ident}/{len(s1)} ({ident/len(s1)*100:.2f}%)")
        log_lines.append(f"  Gapped Columns: {gaps}")
        log_lines.append("")
        
    df = pd.DataFrame(records)
    out_tsv = str(pathlib.Path(__file__).resolve().parent.parent / 'data/processed/hsv1_hsv2_homology_map.tsv')
    df.to_csv(out_tsv, sep='\t', index=False)
    
    out_log = str(pathlib.Path(__file__).resolve().parent.parent / 'logs/homology_mapping.log')
    with open(out_log, 'w', encoding='utf-8') as f:
        f.write("\n".join(log_lines))
        
    print("\n".join(log_lines))
    print(f"Saved homology map to: {out_tsv}")
    print(f"Saved log to: {out_log}")

if __name__ == '__main__':
    main()

