from Bio import SeqIO
from datetime import datetime
import pandas as pd
import gzip
import os
import pandas as pd

"""
    当16S rRNA基因序列数据库需要更新时进行。
"""
## TODO 测试函数成功
def cleanGbff(dbPath=""):
    """
    Input: the path of database
        - folder_name: tax-ncbidata-<date> 
        - https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz 
        - https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.gbff.gz
        - fna.gz has to be unzipped
        dfPath：the folder path that contains database files
    Output: 
        tax-typestrain-{date}.csv
        pandas dataframe
    """
    date = datetime.now().strftime("%Y%m%d")
    # create an empty dataframe with the required columns
    res = pd.DataFrame(columns=["saccver", "sname", "staxonomy", "sstrain", "typestrain"])
    # create a folder named "bacteria.16SrRNA-{date}"
    # folder_name = f"bacteria.16SrRNA-{date}"
    # if not os.path.exists(folder_name):
    #    os.mkdir(folder_name)
    if os.path.exists(dbPath):
        os.chdir(dbPath)
    else:
        Print("Wrong path")
        return -1
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
    return res

def cleanFna(res, dbfile, dbPath=""):
    """
    Input: the path of database
        - folder_name: tax-ncbidata-<date> 
        - https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz 
        - https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.gbff.gz
        - fna.gz has to be unzipped

    Output:

    """
    date = datetime.now().strftime("%Y%m%d")
    if not isinstance(res, isinstance(res, pd.DataFrame)):
        if os.path.exists(dbPath):
            os.chdir(dbPath)
            res = pd.read.csv(dbfile)
        else:
            Print("Wrong path")
            return -1
    # open the fasta file
    records = SeqIO.parse("bacteria.16SrRNA.fna", "fasta")
    # create an empty list to store the type strain records
    typestrain_records = []
    # initialize the index
    i = 0
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
    return typestrain_records


def ftpGetDb():
    """
    能运行，但是非常的慢 
    ## DEBUG 解决ftp下载NCBI速度慢
    """
    import os
    from datetime import datetime
    import ftplib
    import time


    date = datetime.now().strftime("%Y%m%d")

    # create a folder named "bacteria.16SrRNA-{date}"

    folder_name = f"bacteria.16SrRNA-{date}"
    if not os.path.exists(folder_name):
        os.mkdir(folder_name)

    # download the file from the given url and save it in the folder
    # connect to the FTP server
    ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
    ftp.login()
    print("Successfully connected to FTP server")
    ftp.cwd("/refseq/TargetedLoci/Bacteria")
    ftp.maxline = 1048576 # set buffer size to 1 MB
    filename = "bacteria.16SrRNA.fna.gz"
    local_filename = os.path.join(folder_name, filename)

    ftp.retrbinary(f"RETR {filename}", open(local_filename, "wb").write)
    ftp.quit()

def main(folder):
    res = cleanGbff(folder) # 同时会导出一个文件
    if not res == -1:
        cleanFna(res, dbPath=folder) #同时会导出一个fna文件

if __name__ == "__main__":
    import argparse
    # 创建参数解析器
    parser = argparse.ArgumentParser(description='Process some variables.')
    # 添加参数
    parser.add_argument('--folder', type=str, default="", help='the folder containing seq files')    
    
    # 解析参数
    args = parser.parse_args()
    # 运行main函数
    main(args.folder)