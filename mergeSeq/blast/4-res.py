
from Bio import SeqIO
import pandas as pd
from datetime import datetime
import gzip
import os

os.chdir("bacteria.16SrRNA-20230406")
res = pd.read_csv("tax-typestrain-20230406.csv")
i=0
date = datetime.now().strftime("%Y%m%d")
# open the fasta file

records = SeqIO.parse("bacteria.16SrRNA.fna", "fasta")
# create an empty list to store the type strain records
typestrain_records = []
# iterate over the records
for record in records:
    
    # check if the sequence id is in the type strain dataframe
    if res.loc[res["saccver"] == record.id, "typestrain"].any(): #如果找不到就是false,如果找到就是True
        # if it is, add the record to the list
        typestrain_records.append(record)
    else:
        i += 1
           
print(i)
# write the type strain records to a new fasta file
with open(f"16S-type-{date}.fna", "w") as f:
    SeqIO.write(typestrain_records, f, "fasta")