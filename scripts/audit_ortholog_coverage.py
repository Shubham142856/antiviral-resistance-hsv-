"""
Ortholog Alignment Inspection & Coverage Audit
==============================================
Audits the reference-anchored ortholog alignments used in Track B1:
  UL23: HSV-1, HSV-2, VZV (ORF36)
  UL30: HSV-1, HSV-2, VZV (ORF28), HCMV (UL54), HHV-6B (U38)

Measures:
  1. Overall pairwise sequence identity and alignment length
  2. Per-position ortholog coverage (how many orthologs have non-gap residues)
  3. Coverage specifically at the 929 evaluated mutation positions
  4. Sensitivity analysis: impact of excluding poorly-covered positions on conservation AUROC
"""

import os
import sqlite3
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner, substitution_matrices

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PRIMARY_TSV = os.path.join(ROOT, "data", "processed", "primary_position_folds.tsv")

def get_gene_seq(ref, gene_name):
    db_path = os.path.join(ROOT, "herpesdrg", "inst", "ref", f"{ref}.sqlite")
    fasta_path = os.path.join(ROOT, "herpesdrg", "inst", "ref", f"{ref}.fasta")
    genome = str(SeqIO.read(fasta_path, 'fasta').seq)
    conn = sqlite3.connect(db_path)
    query = """
    SELECT t._tx_id, t.tx_name, t.tx_strand, t.tx_start, t.tx_end, c.cds_start, c.cds_end
    FROM transcript t
    JOIN splicing s ON t._tx_id = s._tx_id
    JOIN cds c ON s._cds_id = c._cds_id
    WHERE t.tx_name = ? OR t.tx_name LIKE ?
    ORDER BY s.exon_rank
    """
    rows = conn.execute(query, (gene_name, f'%{gene_name}%')).fetchall()
    if not rows:
        query2 = """
        SELECT t._tx_id, t.tx_name, t.tx_strand, t.tx_start, t.tx_end, c.cds_start, c.cds_end
        FROM gene g
        JOIN transcript t ON g._tx_id = t._tx_id
        JOIN splicing s ON t._tx_id = s._tx_id
        JOIN cds c ON s._cds_id = c._cds_id
        WHERE g.gene_id = ? OR g.gene_id LIKE ?
        ORDER BY s.exon_rank
        """
        rows = conn.execute(query2, (gene_name, f'%{gene_name}%')).fetchall()
    strand = rows[0][2]
    cds_dna = "".join([genome[r[5]-1:r[6]] for r in rows])
    s = Seq(cds_dna)
    if strand == '-':
        s = s.reverse_complement()
    return str(s.translate(to_stop=True))

orthologs = {
    'UL23': {
        'HSV-1': get_gene_seq('NC_001806.2', 'UL23'),
        'HSV-2': get_gene_seq('NC_001798.2', 'UL23'),
        'VZV':   get_gene_seq('NC_001348.1', 'ORF36'),
    },
    'UL30': {
        'HSV-1':  get_gene_seq('NC_001806.2', 'UL30'),
        'HSV-2':  get_gene_seq('NC_001798.2', 'UL30'),
        'VZV':    get_gene_seq('NC_001348.1', 'ORF28'),
        'HCMV':   get_gene_seq('NC_006273.2', 'UL54'),
        'HHV-6B': get_gene_seq('NC_000898.1', 'U38'),
    }
}

aligner = PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -10.0
aligner.extend_gap_score = -0.5
aligner.mode = 'global'

def align_pair(ref_seq, qry_seq):
    aln = aligner.align(ref_seq, qry_seq)[0]
    r_aln, q_aln = aln[0], aln[1]
    
    # Calculate global alignment identity
    matches = sum(1 for r, q in zip(r_aln, q_aln) if r == q and r != '-')
    aligned_len = len(r_aln)
    identity = matches / max(len(ref_seq), len(qry_seq))
    
    # Map ref pos -> query char
    mapping = {}
    r_idx = 0
    for r, q in zip(r_aln, q_aln):
        if r != '-':
            r_idx += 1
            mapping[r_idx] = q
    return identity, mapping, r_aln, q_aln

print("=== 1. GLOBAL PAIRWISE ALIGNMENT METRICS TO HSV-1 ===")
alignments = {}
for g in ['UL23', 'UL30']:
    ref = orthologs[g]['HSV-1']
    print(f"\nTarget Gene: {g} (HSV-1 reference length: {len(ref)} aa)")
    alignments[g] = {}
    for o_name, o_seq in orthologs[g].items():
        if o_name == 'HSV-1':
            continue
        ident, mapping, r_aln, q_aln = align_pair(ref, o_seq)
        gaps_in_qry = sum(1 for c in mapping.values() if c == '-')
        coverage = (len(ref) - gaps_in_qry) / len(ref)
        alignments[g][o_name] = mapping
        print(f"  {o_name:8s} ({len(o_seq)} aa) | Identity: {ident*100:5.2f}% | HSV-1 Coverage: {coverage*100:5.2f}% ({len(ref)-gaps_in_qry}/{len(ref)} non-gap pos)")

print("\n=== 2. COVERAGE AUDIT AT EVALUATED BENCHMARK POSITIONS ===")
df = pd.read_csv(PRIMARY_TSV, sep='\t')

for g in ['UL23', 'UL30']:
    sub_df = df[df['gene'] == g]
    unique_pos = sorted(sub_df['pos'].unique())
    n_pos = len(unique_pos)
    ref = orthologs[g]['HSV-1']
    o_names = [k for k in orthologs[g].keys() if k != 'HSV-1']
    
    print(f"\nGene {g}: {n_pos} unique evaluated positions in benchmark")
    
    # Check coverage at each evaluated position
    cov_counts = {len(o_names): 0}
    for k in range(len(o_names)):
        cov_counts[k] = 0
        
    poorly_covered = []
    for p in unique_pos:
        non_gaps = sum(1 for o in o_names if alignments[g][o].get(p, '-') != '-')
        cov_counts[non_gaps] = cov_counts.get(non_gaps, 0) + 1
        if non_gaps < len(o_names) // 2 + 1:
            poorly_covered.append(p)
            
    print(f"  Ortholog coverage distribution across evaluated positions (total orthologs: {len(o_names)}):")
    for non_gaps, cnt in sorted(cov_counts.items(), reverse=True):
        print(f"    {non_gaps}/{len(o_names)} orthologs aligned: {cnt:3d} positions ({cnt*100/n_pos:5.1f}%)")
    
    print(f"  Positions with < 50% ortholog coverage: {len(poorly_covered)} / {n_pos}")
    if poorly_covered:
        print(f"    Positions: {poorly_covered}")
