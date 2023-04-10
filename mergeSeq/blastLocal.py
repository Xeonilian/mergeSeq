# Importing necessary modules from Biopython
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
import pandas as pd

# Defining the input file and output file names
def blastLocal(input_file="batch_blast.fasta", output_file="batch_blast.csv"):
    
    if os.path.isfile(input_file):
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

def parseTax(blast_res ="batch_blast.csv", tax_file =""):
    if os.path.isfile(blast_res) and os.path.isfile(tax_file):
        res = pd.read(blast_res)
        tax = pd.read(tax_file)