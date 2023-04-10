# 构成 tax-ncbidata-<date> 文件
# 筛选非模式菌株后形成blastdb
#  https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz 
#  https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.gbff.gz
# 构成 tax-ncbidata-<date> 文件

from Bio import SeqIO
import pandas as pd
from datetime import datetime
import gzip
import os


date = datetime.now().strftime("%Y%m%d")
# create an empty dataframe with the required columns
res = pd.DataFrame(columns=["saccver", "sname", "staxonomy", "sstrain", "typestrain"])
# initialize the index
i = 0

# create a folder named "bacteria.16SrRNA-{date}"

folder_name = f"bacteria.16SrRNA-{date}"
if not os.path.exists(folder_name):
    os.mkdir(folder_name)

os.chdir(folder_name)

with gzip.open( "bacteria.16SrRNA.gbff.gz", "rt") as con:
    # read the file line by line
    for line in con:
        # check if the line is empty
        if len(line) == 0:
            break
        # check if the line starts with "//"
        if line.startswith("//"):
            i += 1
        # check if the line contains "VERSION"
        if "VERSION" in line:
            res.loc[i, "saccver"] = line.split("VERSION")[-1].strip()
        # check if the line contains "strain="
        if "strain=" in line:
            res.loc[i, "sstrain"] = line.split('strain="')[1].split('"')[0]
        # check if the line contains "ORGANISM"
        if "ORGANISM" in line:
            name = line.split("ORGANISM")[-1].strip().split(" ")[0:2]
            name = " ".join(name)
            res.loc[i, "sname"] = name
            line1 = con.readline()
            line2 = con.readline()
            line=line1.strip()+" "+line2.strip()
            tax = line[:-1]
            res.loc[i, "staxonomy"] = tax
        # check if the line contains "type strain of"
        if "type strain of" in line:
            res.loc[i, "typestrain"] = True

# close the gbff file
con.close()

# fill the NaN values in the "typestrain" column with False
res["typestrain"].fillna(False, inplace=True)
# filter the dataframe to only include type strains
res = res[res["typestrain"]]
# write the filtered dataframe to a csv file
res.to_csv(f"tax-typestrain-{date}.csv", index=False)

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
