import pandas as pd
import os
from Bio import SeqIO
from datetime import datetime
import gzip
import shutil
import glob
import argparse

"""
    当16S rRNA基因序列数据库需要更新时，从NCBI数据库下载bacteria.16SrRNA.fna和bacteria.16SrRNA.gbff数据。筛选typestain
"""
date = datetime.now().strftime("%Y%m%d")
def cleanGbff(db_path=""):
    """
    清理数据库并保存为csv和db文件
    Args:
        db_path: the folder path that contains database files
    Returns:
        res: 如果成功pandas.dataframe，如果失败-1
    Output: 
        tax-typestrain-{date}.csv
        tax-typestrain-{date}.db
    """
    res = pd.DataFrame(columns=["saccver", "sname", "staxonomy", "sstrain", "typestrain"])

    with gzip.open("bacteria.16SrRNA.gbff.gz", "rt") as con:
        i = 0
        for line in con:
            if len(line) == 0:
                break
            if line.startswith("//"):
                i += 1
            if "VERSION" in line:
                res.loc[i, "saccver"] = line.split("VERSION")[-1].strip()
            if "strain=" in line:
                res.loc[i, "sstrain"] = line.split('strain="')[1].split('"')[0]
            if "ORGANISM" in line:
                name = line.split("ORGANISM")[-1].strip().split(" ")[0:2]
                name = " ".join(name)
                res.loc[i, "sname"] = name
                line1 = con.readline()
                line2 = con.readline()
                line = line1.strip() + " " + line2.strip()
                tax = line[:-1]
                res.loc[i, "staxonomy"] = tax
            if "type strain of" in line:
                res.loc[i, "typestrain"] = True
    con.close()

    res["typestrain"].fillna(False, inplace=True)
    res = res[res["typestrain"]]
    res.to_csv(f"tax-typestrain-{date}.csv", index=False)

    import sqlite3
    conn = sqlite3.connect(f"tax-typestrain-{date}.db")
    res.to_sql('mytable', conn, index=False)
    conn.close()
    return res

def cleanFna(res="", db_path=""):
    """
    Args:
        res: 直接输入pandas.DataFrame。
        gb_file: 通过存有gb_file的路径，读取pandas.DataFrame。
        db_path: 工作文件夹的路径，其中包括fna文件。
    Returns:
        typestrain_records:以seqRecord的生成器返回。
    Output:
        16S-type-{date}.fna: 仅有type strain的fna文件。
    """
    date = datetime.now().strftime("%Y%m%d")
    if os.path.exists(db_path):
        os.chdir(db_path)
        if not isinstance(res, pd.DataFrame):
            file = glob.glob("tax-typestrain-*.csv")[0]
            res = pd.read_csv(os.path.join(db_path,file))
    else:
        print("Wrong path")

    # open the compressed file in read mode
    with gzip.open('bacteria.16SrRNA.fna.gz', 'rb') as f_in:
        with open('bacteria.16SrRNA.fna', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    records = SeqIO.parse("bacteria.16SrRNA.fna", "fasta")
    typestrain_records = []
    i = 0
    for record in records:
        if res.loc[res["saccver"] == record.id, "typestrain"].any():
            typestrain_records.append(record)
        else:
            i += 1
    print(i)
    with open(f"16S-type-{date}.fna", "w") as f:
        SeqIO.write(typestrain_records, f, "fasta")
    return typestrain_records

def getDb(db_path="", fna=True, gbff=True):
    """
    从NCBI的ftp数据库下载用request下载16SrRNA 基因数据库信息
    Args:
        db_path:用来存放下载数据的路径。
        fna:是否下载fna。
        gbff: 是否下载gbff。
    Output:
        - bacteria.16SrRNA.fna.gz 
        - bacteria.16SrRNA.gbff.gz 
        - folder: bacteria.16SrRNA-{date}
    Note:
        - 使用request下载
    """
    # download the file from the given url and save it in the folder
    import requests
    if fna:
        print("Try connect fna")
        url = "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz"
        filename = "bacteria.16SrRNA.fna.gz"
        response = requests.get(url, stream = True)
        print("connected")
        if response.status_code == 200:
            
            # Do something with the response content
            with open(os.path.join(db_path, filename), "wb") as f:
                f.write(response.content)
            print(f"downloaded {filename}")
        else:
            print("Error: Could not retrieve fna from NCBI")
    if gbff:
        url = "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.gbff.gz"
        filename = "bacteria.16SrRNA.gbff.gz"
        print("Try connect gbff")
        response = requests.get(url, stream = True)
        print("connected")
        if response.status_code == 200:   
            with open(os.path.join(db_path, filename), "wb") as f:
                total_size = int(response.headers.get('content-length', 0))
                chunk_size = 8192
                downloaded_size = 0 # initialize downloaded_size
                for chunk in response.iter_content(chunk_size=chunk_size):
                    f.write(chunk)
                    downloaded_size += len(chunk)
                    progress = downloaded_size / total_size * 100
                    print(f"Downloaded {downloaded_size} bytes ({progress:.2f}%)")
            print(f"Downloaded {filename}")
        else:
            print("Error: Could not retrieve gbff from NCBI")

def main():
    """
    运行获取数据库，清理数据得到只有模式菌株的数据。
    Args: 
        db_path: 用来存放NCBI数据库的文件夹的路径。
        db_folder: 存放数据的文件夹名，如果改参数不为空，且文件夹不存在时新建文件夹。
        db: 是否下载数据库文件，默认False。
        fna: 是否下载数据库中的fna，当db为True时有效，默认True。
        gbff: 是否下载数据库中gbff文件，当db为True时有效，默认为True。
    OutPut:
        - tax-typestrain-{date}.csv
        - 16S-type-{date}.fna
        - bacteria.16SrRNA.fna
        - bacteria.16SrRNA.fna.gz <optional>
        - bacteria.16SrRNA.gbff.gz <optional>
    """
    # 创建参数解析器
    parser = argparse.ArgumentParser()
    # 添加参数
    parser.add_argument('--db_path', type=str, default="", help='the path used to store db')    
    parser.add_argument('--db_folder', type=str, default="", help='the folder used to store db')
    parser.add_argument('--db', type=bool, default="False", help='whether download database from NCBI')
    parser.add_argument('--fna', type=bool, default="True", help='download fna file')
    parser.add_argument('--gbff', type=bool, default="True", help='download gbff file')
    parser.add_argument('--test', type=bool, default="False", help=' do test')
    # 解析参数
    args = parser.parse_args()
    # 测试用 完成调试
    if args.test:
        args.db_path="E:\\Desktop\\python-learn-2023\\4-mergeSeq\\mergeSeq\\examples"
        args.db_folder="16SDB-20230412"
        args.db=False
        args.fna=False
        args.gbff=False
    

    if os.path.exists(args.db_path): # 存在指定路径
        if not args.db_folder=="": #新建一个文件夹，并设置为当前路径
            os.chdir(args.db_path)
            if os.path.exists(args.db_folder):
                print("Folder already exists")
                db_path = os.path.join(args.db_path, args.db_folder)
            else:
                os.mkdir(args.db_folder)
                print("Folder created")
                db_path = os.path.join(args.db_path, args.db_folder)
            os.chdir(db_path)
        
        if args.db: #下载到
            print("download db from NCBI")
            try:
                getDb(db_path=db_path, fna=args.fna, gbff=args.gbff)
            except:
                print("error: not downloaded!")
        else:
            print("skip db download")
        if os.path.exists("bacteria.16SrRNA.gbff.gz") and os.path.exists("bacteria.16SrRNA.fna.gz"):
            if any(filename.startswith("tax-typestrain") for filename in os.listdir(db_path)):
                print("tax-typestrain file import \nstart fna clean")
                cleanFna(db_path=db_path)
                print("fna clean") 
            else:
                print("start gbff clean")
                res = cleanGbff(db_path)
                print("gbff clean \nstart fna clean")
                cleanFna(res, db_path)
                print("fna clean")            
        else:
            print("Need fna file and gbff file.\nSet --db True can download files from NCBI")
    else:
        if args.db_path=="":
            print("db_path is requeried, can set using --db_path <path>")
        else:
            print("Error: wrong path")          

if __name__ == "__main__":  
    main()
