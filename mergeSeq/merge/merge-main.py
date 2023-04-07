import os
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import defaultdict

################## 函数定义 ###################
def ab1IO(path_var: list):
    '''
    功能：导入ab1文件
    '''
    for record in SeqIO.parse(path_var['file_name'], "abi"):
        record.primer = path_var['primer']
        record.sample = path_var['sample']
        record.type = path_var['file_type']
        record.description += f"{record.type}|"
        print(f"import ab1:{record.sample} {record.primer}")
        return record

def fasIO(path_var: list):
    '''
    功能：导入fasta文件
    '''
    try:  # 试了试文件对不对
        next(SeqIO.parse(path_var['file_name'], "fasta"))
        for record in SeqIO.parse(path_var['file_name'], "fasta"):  # 如果是fasta文件
            record.primer = path_var['primer']
            record.sample = path_var['sample']
            record.type = path_var['file_type']
            record.description += f"{record.type}|"
            print(f"import fasta:{record.sample} {record.primer}")
            return record
    except Exception as e:
        with open(path_var['file_name'], "r") as f:
            file_contents = f.read()
        record = SeqRecord(seq=Seq(file_contents),
                           id=path_var['name'],
                           name=path_var['name'],
                           description="add >")
        record.primer = path_var['primer']
        record.sample = path_var['sample']
        record.type = path_var['file_type']
        record.description += f"{record.type}|"
        
        print(f"import seq:{record.sample} {record.primer}")
        return record

def pathPaser(file_name: str, primers:list):
    '''
    功能：将文件名解析成sample，primer；确定primer在不在primers中，也就是保证数据成对才导入
    '''
    path_var = {}
    try:
        path_var['file_type'] = os.path.splitext(os.path.basename(file_name))[-1]
        path_var['name'] = os.path.splitext(os.path.basename(file_name))[0]
        match = re.search(r"(.*)_\[(\w+)\]", path_var['name'])
        if match.group(2) not in primers:

            path_var['check'] = False
        else:
            path_var['primer'] = match.group(2)
            path_var['sample'] = match.group(1)
            path_var['file_name'] = file_name
            path_var['check'] = True
    except:

        path_var['check'] = False
    return path_var

def checkRecordDict(record_dict, path_var: list):
    # 判断在record_dict中path_var中保存的sample和primer是不是存在
    if path_var['sample'] in record_dict:
        for record in record_dict[path_var['sample']]:
            if record.primer == path_var['primer']:
                print(f"pass: {path_var['file_name']}")
                return True
    return False

def checkFiles(file_names, primers):
    """
    功能：生成一个字典存放导入的结果，根据文件名，以及待拼接的成对引物 
    """ 
    record_dict = defaultdict(list)
    # print(type(record_dict))
    for file_name in file_names:
        # print(file_name)
        path_var = pathPaser(file_name, primers)
        if path_var['check'] == True:
            record_dict[path_var['sample']] = []
    return record_dict

def showRecordDict(record_dict):
    """
    查看数据结构
    """
    for sample in record_dict.keys():
        print(sample, len(record_dict[sample]))
        for record in record_dict[sample]:
            print("  ", record.primer, record.type, len(record.seq))

def checkPair(record_dict):
    """
    清除不成对的数据
    """
    record_dict_check = {k: v for k, v in record_dict.items() if isinstance(v, list) and len(v) == 2}
    return record_dict_check

def Phred2Index(record, check_len=0, thredhold=30, break_num=3):
    # 将Phred转化为index
    index = [-1, -1]  # 如果是空的，后面赋值就会有问题，这里有-1+1=0也没有办法通过最后的长度检测
    # 根据阈值转变为boolean值
    Phred_list = record.letter_annotations["phred_quality"]
    if not Phred_list:
        return index   
    Phred_list = [(x-thredhold) > 0 for x in Phred_list]
    # 找到大于阈值的起点
    for i, num in enumerate(Phred_list):
        if num:
            index[0] = i
            break
    Phred_list = Phred_list[index[0]:]
    for i, num in enumerate(Phred_list):
        if not num:
            if not any(Phred_list[i:i+break_num]):
                index[1] = i
                break
    # print(not check_len)
    if not check_len:  # 如果 = 0 这里就是True
        return index
    else:
        if index[1]-index[0] > check_len:
            return index
        else:
            return [index[0], -1]
 
def Ab12Fas(record):
    # ab1转化为fasta  
    if record.type == ".ab1":
        record_fas = SeqRecord(record.seq, id=record.id,
                               description=record.description, 
                               name=record.name)
        record_fas.sample = record.sample
        record_fas.primer = record.primer
        record_fas.type = record.type
        return record_fas
    else:
        return record

def trimFasta(record, index=[20, 750]):
    # 基于index剪切fasta数据
    record_trim = record
    if len(record.seq) > index[1] > index[0]:
        record_trim.seq = record.seq[index[0]:index[1]]
        record_trim.description += f"trimmed {index[0]}:{index[1]}|"
        record_trim.index = index
    else:
        record.index = [-1,-1]
    return record_trim

def swapSeq(seq1, seq2):
    '''
    序列顺序交换
    '''
    seq = seq1
    seq1 = seq2
    seq2 = seq
    return seq1, seq2

def findLong(aligned_segs):
    """
    alignedSeges 三层 tuple，从Bio.Align.PairwiseAlignment对象的aligned获得
    返回其中最长的片段的对应的target的结束位置，query的结束位置
    """
    marks = []
    seg_index = []
    if len(aligned_segs[0]) ==1:
        marks.append(aligned_segs[0][0][1])
        marks.append(aligned_segs[1][0][1])
        
    else:
        for i in range(0,len(aligned_segs)):
            seg_index.append(aligned_segs[0][i][1]-aligned_segs[0][i][0])
        pos = seg_index.index(max(seg_index))
        marks.append(aligned_segs[0][pos][1])
        marks.append(aligned_segs[1][pos][1])
    return marks

def mergeSeq(seq1, seq2):
    """
    输入2个Seq对象，第二条反补，然后拼接。
    返回Seq对象的序列，返回相似度，联配长度构成的list
    """
    # 反补
    seq2 = seq2.reverse_complement()
    # 建立联配
    aligner = Align.PairwiseAligner()
    # Set the aligner to perform local alignments
    aligner.mode = 'local'
    aligner.match_score = 2
    aligner.mismatch_score = -2
    aligner.gap_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    alignments = aligner.align(seq1, seq2)

    # 解析结果
    if len(alignments) == 1:
        res = alignments[0]
        marks = findLong(res.aligned)
        seqM = res.target[:marks[0]]+res.query[marks[1]:]
        similarity = res.score/res.shape[1]/2*100
        aligned_len = res.shape[1]
        return [seqM, similarity, aligned_len, len(seqM)]
    else:
        return len(alignments)

def main(file_names, primers=['27F', '1492R'], mod="auto", check_len=0, thredhold=30, break_num=3, fastacut=[20, 750], verbose=True，out="none"):
    """
    功能：
        1. 输入文件，优先导入ab1，生成record_dict
        2. 清理序列字典，选择自动或手动将每条序列剪切到需要的范围生成merge_dict: {sample:[record1,record2],}
        3. 联配，遍历每个key，其中成对的序列进行拼接，拼接结果追加到当前的key中，得到merge_dict，选择导出序列
    参数
        primers默认为27F和1492R
        mod: auto or mannual
        check_len: 剪切后短于这个数值的不要
        threhold：phred的阈值，目前20和40都不对
        break_num：连续几个错误认为不再要
        fastacut：fasta序列的默认剪切范围
        verbose:是不是要显示数据处理的全过程
        out:导出数据：none-不导出，merge-只导出merge结果，all-全部导出为1个文件
    输出：
        merge_dict{sample:[record1,record2，record3],}
    """
    # 1 输入阶段--------
    # 文件排序
    file_names.sort(key=lambda x: (not x.endswith('ab1'), x.endswith('ab1') and x))
    # 解析文件名构建保存数据的字典
    record_dict = checkFiles(file_names, primers)
    # 遍历文件导入
    for file_name in file_names:
        # 解析文件名 如果无法解析到primer信息，就停止导入
        path_var = pathPaser(file_name, primers)
        if path_var['check'] == True:
            # 判断该文件是否存在
            if not checkRecordDict(record_dict, path_var):
                if path_var['file_type'] == ".ab1":  
                    record = ab1IO(path_var)
                    record_dict[path_var['sample']].append(record)
                else:
                    record = fasIO(path_var)
                    record_dict[path_var['sample']].append(record)
    if verbose:
        showRecordDict(record_dict)
    # 2 剪切阶段----------
    record_dict_clean = defaultdict(list)
    record_dict_clean = checkPair(record_dict)
    # 保存结果用字典
    merge_dict = {k: [] for k in record_dict_clean.keys()}
    if mod == "auto":
        # Index生成
        for sample in record_dict_clean.keys():
            for record in record_dict_clean[sample]:
                if record.type == ".ab1":                  
                    ab1cut = Phred2Index(record, check_len, thredhold, break_num)
                    merge_dict[sample].append(trimFasta(Ab12Fas(record), index=ab1cut))
                else:
                    merge_dict[sample].append(trimFasta(record, index=fastacut))
        if verbose:
            print("-----------cut-----------")
            showRecordDict(merge_dict)  
    # 3 拼接阶段-----------
    for sample in merge_dict.keys():
        seqF = merge_dict[sample][0]
        seqR = merge_dict[sample][1]
        if "F" in seqR.primer and "R" in seqF.primer:
            seqF, seqR = swapSeq(seqF, seqR)
        seqM = mergeSeq(seqF.seq, seqR.seq)  # 需要转化为Seq对象
        description = f"len={str(seqM[3])}bp|{str(seqM[1])}% aligned of {str(seqM[2])}bp|{primers[0]} {primers[1]}"
        merged_seq = SeqRecord(seq=seqM[0], name=sample+"_[merge]", id=sample+"_[merge]",
                               description=description)
        merged_seq.primer = "merge"
        merged_seq.type = ".fasta"
        merge_dict[sample].append(merged_seq)
        # 拼接完成后是否输出summary
        if verbose:
            rint("-----------merge-----------")
            showRecordDict(merge_dict)
        # 拼接完后是否导出数据
        if out == "none":
            return merge_dict
        elif out == "merge":
            SeqIO.write(merged_seq, merged_seq.name+merged_seq.type, "fasta")
        elif out == "all":
            with open(f"{sample}.fasta", "a") as output_file:
                # Loop over the sequences you want to append
                for seq in merge_dict[sample]:
                    # Write each sequence to the output file
                    SeqIO.write(seq, output_file, "fasta")
    return merge_dict