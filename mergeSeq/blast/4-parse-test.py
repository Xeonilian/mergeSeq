import pandas as pd
import os
os.chdir("bacteria.16SrRNA-20230406")
res = pd.DataFrame(columns=["saccver", "sname", "staxonomy", "sstrain", "typestrain"])
i=0
with open( "gbff-test.gbff", "rt") as con:
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

res