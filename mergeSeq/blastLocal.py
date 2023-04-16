from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
import pandas as pd
import argparse
import os
from io import StringIO
"""
    1. 本地16S数据库blast
    2. 成对序列的blast

"""


def blastPair(subject_file, query_file, out, outfmt):
    """
    比对两条序列 
    Args:
    output:
        blast result default name is "pari_blast.txt" in outfmt 10
    """
    if os.path.isfile(subject_file) and os.path.isfile(query_file):
        # Performing BLAST locally and saving the output in XML format
        blastn_cline = NcbiblastnCommandline(query=query_file,
                                             subject=subject_file,  outfmt=outfmt)
        # Performing BLAST locally and saving the output in XML format
        try:
            stdout, stderr = blastn_cline()
        except Exception as e:
            stdout = None
            print("An error occurred:", e)
        if stdout and out:
            with open(out, 'w') as out_file:
                out_file.write(stdout)
        else:
            print(stdout)
    else:
        print("Error: wrong path!")

def blastLocal(input_file, out, db, db_folder):
    """
    将序列和本地数据库进行blast
    """
    
    if os.path.isfile(input_file):
        if os.path.exists(db_folder):
            os.chdir(db_folder)
            blastn_cline = NcbiblastnCommandline(query=input_file, db=db,  outfmt="10 qseqid sseqid pident length bitscore", num_alignments=1)
        # Performing BLAST locally and saving the output in table format
        try:
            stdout, stderr = blastn_cline()
        except Exception as e:
            stdout = None
            print("An error occurred:", e)

        if stdout:
            # Converting stdout to pandas dataframe
            blast_res = pd.read_csv(StringIO(stdout), sep=',', header=None)
            blast_res.columns = ['qseqid', 'saccver', 'pident', 'length', 'bitscore']
            return blast_res
            if out:
            # Joining the current working directory and output_file to create a path
                path = os.path.join(os.getcwd(), out)
                # Saving the stdout into the joined path
                with open(path, 'w') as out_file:
                    out_file.write(stdout)

def parseTax(blast_res, tax_file):
    """
    解析blast结果中的分类地位
        Args:
        blast_res: pandas DataFrame containing the following columns:
            - qseqid: query sequence ID
            - saccver: subject sequence accession.version
            - pident: percentage of identical matches
            - length: alignment length
            - bitscore: bit score
        tax_file: path to a file containing taxonomy information for the subject sequences, with the following columns:
            - saccver: subject sequence accession.version
            - sname: subject sequence name
            - staxonomy: subject sequence taxonomy
            - sstrain: subject sequence strain
            - typestrain: whether the subject sequence is a type strain
        Returns:
        tax_res: pandas DataFrame containing the following columns:
            - qseqid: query sequence ID
            - accver: subject sequence accession.version
            - taxonomy: taxonomy information for the subject sequence
            - name: subject sequence name
            - strain: subject sequence strain
            - similarity: percentage of identical matches
            - length:  alignment length
            - bitscore: bit score
            - db-version: version of the subject sequence database
    """
    
    tax_res = blast_res
    import sqlite3
    # Check if the tax_file exists
    if os.path.isfile(tax_file):
        db_name = os.path.splitext(os.path.basename(tax_file))[0]
        conn = sqlite3.connect(tax_file)
        # Query the database for rows where the 'saccver' column matches any value in the 'saccver' column of matching_rows
        matching_rows = pd.read_sql_query(f"SELECT * FROM '{db_name}' WHERE saccver IN ({','.join(['?']*len(blast_res))})", 
            conn, params=blast_res['saccver'].tolist())
 
        conn.close()
        # Merge the matching_rows DataFrame with the tax_res DataFrame on the 'saccver' column
        if matching_rows.empty:
            print("No matching!")
            return -1
        else:
            tax_res = pd.merge(tax_res, matching_rows, on='saccver')
            tax_res = tax_res.drop(columns='typestrain')
            tax_res['db-version'] = db_name
            return tax_res
    else:
        print("Error: wrong path!")

def __getDescription__(tax_res_row):
    '''
    将一行分类比对信息生产description追加内容
        Args:
            tax_res_row: <pandas.dataframe> 一行数据
        Returns:
            des_res: <list> 包括追加的description和文件名。
            - description 追加内容：|taxonomy=xx; xx; xx; xx; xx; name=Desulfolithobacter dissulfuricans strain=GF1; similarity=xx.x%; db=16S-xxxxxx。
            - sequence file name: name similarity  31321032400942_(23)_[xxx.xxx 99.9].fasta。
    '''
    tax = tax_res_row.iloc[0].tolist()
    des_string = f"|taxonomy={tax[6]}; name={tax[5]}; strain={tax[7]}; similarity={tax[2]:.0f}%; db-version={tax[8]}"
    # Generate the fasta file name with taxonomic information
    fasta_file = f"{tax[0]}_[{tax[5]} {tax[2]}].fasta"
    return [des_string, fasta_file]

def __parseDescription__(seq_description):
    '''
    解析序列的描述，确认是否比对过，以及比对的数据库版本
    '''
    if "db-version=" in seq_description:
        return seq_description.split("db-version=")[-1]
    else:
        return False

def __updateFastafiles__(tax_res, input_file, redo):
    """
    增加或更新description之后输出fasta文件
        Args:
           tax_res: <pandas.Dataframe> 保存有序列分类地位信息的表格
           input_file: <str> 对应序列的fasta文件
           redo: <bool> 是否替换已有的description 
        Output:
            序列名更新文件名，序列中追加或者替换现有的数据库
    """

    dir_path = os.path.dirname(input_file)
    for record in SeqIO.parse(input_file, "fasta"):  
        if __parseDescription__(record.description) : # 如果已经解析过taxonomy
            if redo: # 如果需要更新数据
                if not __parseDescription__(record.description) == tax_res[tax_res['qseqid']==record.id]['db-version'][0]:
                    # 如果数据库没有重复且需要更新，去掉已有的信息，追加新的信息
                    record.description = record.description.split("|taxonomy=")[0]+__getDescription__(tax_res[tax_res['qseqid']==record.id])[0]
                    print(f"updated {record.id}")
                else: #数据库重复
                    print(f"no need to update {record.id}")
            else:
                print(f"pass {record.id}")
        else:
            record.description += __getDescription__(tax_res[tax_res['qseqid']==record.id])[0]              
            print(f"indentified {record.id}")
        with open(os.path.join(dir_path,__getDescription__(tax_res[tax_res['qseqid']==record.id])[1]), "w") as f:
            SeqIO.write(record, f, "fasta")

def main():
    """
    运行以本地16S数据库为模板的blast，解析分类地位，并将分类地位写入，fasta文件的description中
        Args:
            mode: <str>
            subject_file:
            query_file:
            db: <str> 数据库所在路径，并解析数据库的version \\TODO: 以相对引用的方式引用example中的数据库
            input_file: <str> blast的输入结果 \\TODO:以相对引用方式引用sequence中的
            out: <str> 导出的用于分类地位解析的文件。
            outfmt:<str>
            tax_file: <str> 解析后保存结果的文件。
            taxonomy: <bool> 是否检解析分类地位。
            redo: <bool> 强制重新比对，默认False。
            single: <bool> 如果多条序列比对是否保存为一个文件，默认True，保存为单个文件。
        Output:    
            sequence file: blasted_sequence.fasta，或者多条序列文件，文件名根据blast结果命名。
        Note:
            - description 追加内容：|taxonomy=xx; xx; xx; xx; xx; name=Desulfolithobacter dissulfuricans strain=GF1; similarity=xx.x%; db=16S-xxxxxx。
            - sequence file name: name similarity  31321032400942_(23)_[xxx.xxx 99.9].fasta。
            - 数据库的版本依靠输入文件的文件名区分，不代表使用的数据库本身内容是否不同。
    """
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--mode', default="db", type=str, choices=['pairwise', 'db'], 
        help='Choose between pair-wise BLAST or batch BLAST and parse results')
    parser.add_argument(
        '--outfmt', type=str, default="10", help='Output format for BLAST results')
    parser.add_argument(
        '--out', type=str, help='Output file name for BLAST results')
    parser.add_argument(
        '--subject_file', type=str, help='Path to the subject file for pair-wise BLAST')
    parser.add_argument(
        '--query_file', type=str, help='Path to the query file for pair-wise BLAST or batch BLAST')
    parser.add_argument(
        '--db', type=str, default="16S", help='the name of database')
    parser.add_argument(
        '--db_folder', type=str,default="",help="Path to folder having local database")
    parser.add_argument(
        '--input_file', type=str, default="output.fasta", help='Path to the input file')
    parser.add_argument(
        '--tax_file', type=str, help='Path to the tax file, the taxonomy of type strain')
    parser.add_argument(
        '--taxonomy', type=bool, default=True, help='Whether to parse taxonomy')
    parser.add_argument(
        '--redo', type=bool, default=False, help='Whether to redo the search')
    
    parser.add_argument(
        '--test', type=bool, default=False, help='run text setting')
    args = parser.parse_args()
    # 测试
    # args.test = True #测试完成 \\Todo: 修改为相对路径    
    if args.test:
        args.db_folder ="E:\\Desktop\\python-learn-2023\\4-mergeSeq\\mergeSeq\\examples\\16SDB-20230412"
    #    args.input_file ="E:\\Desktop\\python-learn-2023\\4-mergeSeq\\mergeSeq\\examples\\sequences\\output.fasta"
        args.input_file ="E:\\Desktop\\python-learn-2023\\4-mergeSeq\\mergeSeq\\examples\\sequences\\31320101400400_(11)_[merge]_[Mameliella sediminis 100.0].fasta"
        args.tax_file = "E:\\Desktop\\python-learn-2023\\4-mergeSeq\\mergeSeq\\examples\\16SDB-20230412\\tax-typestrain-20230413.db"
        args.redo = True  
    # 运行主体部分
    if args.mode == 'db':
        blast_res = blastLocal(input_file=args.input_file, out=args.out,  db=args.db, db_folder=args.db_folder)
        if args.taxonomy:
            tax_res = parseTax(blast_res=blast_res, tax_file=args.tax_file)
            __updateFastafiles__(tax_res=tax_res,input_file=args.input_file, redo=args.redo)
    elif args.mode == "pairwise":    
        blastPair(args.subject_file, args.query_file, args.outfmt)
    else:
        print("--mode is wrong")

if __name__ == "__main__":
    main()
