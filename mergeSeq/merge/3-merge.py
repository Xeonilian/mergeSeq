from Bio import Align
from Bio.Seq import Seq

# Read in the sequences from the input files
seq1 = Seq('TGCAGTCGAGCGCACTCTTCGGAGTGAGCGGCGGACGGGTTAGTAACGCGTGGGAACGTGCCCAGATCTAAGGAATAGCCACTGGAAACGGTGAGTAATACCTTATACGCCCTTCGGGGGAAAGATTTATCGGATTTGGATCGGCCCGCGTTAGATTAGGTAGTTGGTGGGGTAATGGCCTACCAAGCCTACGATCTATAGCTGGTTTTAGAGGATGATCAGCAACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATCTTGGACAATGGGCGCAAGCCTGATCCAGCCATGCCGCGTGAGTGATGAAGGCCCTAGGGTCGTAAAGCTCTTTCGCCAGGGATGATAATGACAGTACCTGGTAAAGAAACCCCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGGGTTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCGCGTAGGCGGACTGGAAAGTTGGGGGTGAAATCCCGGGGCTCAACCCCGGAACTGCCTCCAAAACTCCCAGTCTTGAGTTCGAGAGAGGTGAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAGGAACACCAGTGGCGAAGGCGGCTCACTGGCTCGATACTGACGCTGAGGTGCGAAAGTGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACACCGTAAACGATGAATGCCAGTCGTCGGGTAGCATGCTATTCGGTGACACACCTAACGGATTAAGCATTCCGCCTGG')
seq2 = Seq('ATTCGGAGGAACACCAGTGGCGAAGGCGGCTCACTGGCTCGATACTGACGCTGAGGTGCGAAAGTGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACACCGTAAACGATGAATGCCAGTCGTCGGGTAGCATGCTATTCGGTGACACACCTAACGGATTAAGCATTCCGCCTGGGGAGTACGGTCGCAAGATTAAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGCAGAACCTTACCAACCCTTGACATCCTCGGACCGGCCCAGAGATGGGTCTTTCACTTCGGTGACCGAGTGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTCGGTTAAGTCCGGCAACGAGCGCAACCCACATCCTTAGTTGCCAGCAGTTCGGCTGGGCACTCTAGGGAAACTGCCCGTGATAAGCGGGAGGAAGGTGTGGATGACGTCAAGTCCTCATGGCCCTTACGGGTTGGGCTACACACGTGCTACAATGGCAGTGACAATGGGTTAATCCCAAAAAACTGTCTCAGTTCGGATTGTCGTCTGCAACTCGGCGGCATGAAGTCGGAATCGCTAGTAATCGCGTAACAGCATGACGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTTGGGTTTACCCGAAGGCCGTGCGCCAACCTTTGGAGGC')

# Create a pairwise aligner object
aligner = Align.PairwiseAligner()

# Set the aligner to perform local alignments
aligner.mode = 'local'
import argparse

parser = argparse.ArgumentParser(description='Set the desired parameters')
parser.add_argument('--match_score', type=float, default=2, help='match score')
parser.add_argument('--mismatch_score', type=float, default=-2, help='mismatch score')
parser.add_argument('--gap_score', type=float, default=-1, help='gap score')
parser.add_argument('--open_gap_score', type=int, default=-5, help='open gap score')
parser.add_argument('--extend_gap_score', type=float, default=-0.5, help='extend gap score')
args = parser.parse_args()



# Set the desired parameters
aligner.match_score = args.match_score
aligner.mismatch_score = args.mismatch_score
aligner.gap_score = args.gap_score
aligner.open_gap_score = args.open_gap_score
aligner.extend_gap_score = args.extend_gap_score

# Display the current parameter settings to stdout
print(f"Current parameter settings: match_score={args.match_score}, mismatch_score={args.mismatch_score}, gap_score={args.gap_score}, open_gap_score={args.open_gap_score}, extend_gap_score={args.extend_gap_score}\n")


# Align seq1 and seq2_rc using the local alignment mode
alignments = aligner.align(seq1, seq2)
print(len(alignments))
# Print the alignment score and the aligned sequences
for alignment in alignments:
    print("Score:", alignment.score)
    print(alignment)
    print(len(alignment))
    

merged_seq = ''
for alignment in alignments:
    merged_seq += alignment.query
    merged_seq += alignment.target[alignment.target_begin:alignment.target_end]
merged_seq = Seq(merged_seq)
print("Merged sequence:", merged_seq, len(merged_seq))




