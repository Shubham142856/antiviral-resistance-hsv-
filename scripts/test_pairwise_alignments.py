from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
from test_extract_orthologs import get_gene_seq

aligner = PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -10.0
aligner.extend_gap_score = -0.5
aligner.mode = 'global'

hsv1_ul23 = get_gene_seq('NC_001806.2', 'UL23')
hsv2_ul23 = get_gene_seq('NC_001798.2', 'UL23')
vzv_ul23  = get_gene_seq('NC_001348.1', 'ORF36')

print("\n--- UL23 Alignments to HSV-1 ---")
for name, seq in [('HSV-2', hsv2_ul23), ('VZV', vzv_ul23)]:
    alignments = aligner.align(hsv1_ul23, seq)
    best = alignments[0]
    print(f'{name} score: {best.score:.1f}')
    print(f'  Target (HSV-1) length: {len(hsv1_ul23)}')
    print(f'  Query ({name}) length: {len(seq)}')
    # Check aligned length
    hsv1_aligned, other_aligned = best[0], best[1]
    print(f'  Alignment length: {len(hsv1_aligned)}')

hsv1_ul30 = get_gene_seq('NC_001806.2', 'UL30')
hsv2_ul30 = get_gene_seq('NC_001798.2', 'UL30')
vzv_ul30  = get_gene_seq('NC_001348.1', 'ORF28')
hcmv_ul30 = get_gene_seq('NC_006273.2', 'UL54')
hhv6_ul30 = get_gene_seq('NC_000898.1', 'U38')

print("\n--- UL30 Alignments to HSV-1 ---")
for name, seq in [('HSV-2', hsv2_ul30), ('VZV', vzv_ul30), ('HCMV', hcmv_ul30), ('HHV-6B', hhv6_ul30)]:
    alignments = aligner.align(hsv1_ul30, seq)
    best = alignments[0]
    print(f'{name} score: {best.score:.1f}')

