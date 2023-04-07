
# Importing necessary modules from Biopython
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO

# Defining the input file and output file names
input_file = "0005_31320101300621_(5)_[merge].fasta"
output_file = "output.csv"

# Performing BLAST locally and saving the output in XML format
blastn_cline = NcbiblastnCommandline(query=input_file, db="16S", out=output_file, outfmt="10 qseqid sseqid pident length bitscore", num_alignments=1 )
# Performing BLAST locally and saving the output in XML format
try:
    stdout, stderr = blastn_cline()
except Exception as e:
    stdout = None
    print("An error occurred:", e)

if stdout:
    with open(output_file, 'w') as out_file:
        out_file.write(stdout)
