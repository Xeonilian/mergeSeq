import os
import re
from Bio import SeqIO, Align
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import defaultdict
from tkinter import filedialog, Tk

def ab1IO(path_var: list, verbose=True):
    '''
    导入ab1文件
    Args:
        path_var: <list> 解析后的路径信息。
        verbose: <bool> 是否打印导入的过程，默认True。
    Returns:
        record: <SeqRecord> 导入的记录。
    '''
    for record in SeqIO.parse(path_var['path'], "abi"):
        record.primer = path_var['primer']
        record.sample = path_var['sample']
        record.type = path_var['file_type']
        record.description += f"{record.type}|"
        if verbose:
            print(f"import ab1:{record.sample} {record.primer}")
        return record

def fasIO(path_var: list, verbose=True):
    '''
    导入1条fasta文件，如果文件格式不符合">"开头，对文件进行修正后导入。
    
    Args:
        path_var: <list> 解析后的路径信息，包括文件名、样本名、引物名和文件类型。
        verbose: <bool> 是否打印导入的过程，默认True。
        
    Returns:
        record: <SeqRecord> 导入的记录。
        
    Note: 
        修正时未检测是不是满足DNA序列的特征。仅导入第一条序列。
    '''
    try:  # 试了试文件是否是fasta文件
        next(SeqIO.parse(path_var['path'], "fasta"))
        for record in SeqIO.parse(path_var['path'], "fasta"):  # 如果是fasta文件
            record.primer = path_var['primer']
            record.sample = path_var['sample']
            record.type = path_var['file_type']
            record.description += f"{record.type}|"
            if verbose:
                print(f"import fasta:{record.sample} {record.primer}")
            return record
    except Exception as e:
        with open(path_var['path'], "r") as f:
            file_contents = f.read()
        record = SeqRecord(seq=Seq(file_contents),
                           id=path_var['name'],
                           name=path_var['name'],
                           description="add >")
        record.primer = path_var['primer']
        record.sample = path_var['sample']
        record.type = path_var['file_type']
        record.description += f"{record.type}|"
        
        if verbose:
            print(f"import seq:{record.sample} {record.primer}")
        return record

def pathPaser(path: str, primers:list,  regex_para:list, verbose=True) -> dict:
    '''
    根据引物对信息解析文件名信息。
    
    Args:
        path: <str> 文件的路径。
        primers: <list> 引物列表。
        regex_para:<list> primer所属group，sample所属的group，正则表达式，默认为生工序列的解析方式。
        
    Returns:
        path_var: <dict> 包含文件类型、文件名、样本名、引物名和'check'键的值。
            如果引物不在列表中，则'check'键的值为False。
            如果引物在列表中，则'check'键的值为True。      
    '''
    path_var = {}
    try:
        path_var['file_type'] = os.path.splitext(os.path.basename(path))[-1]
        path_var['name'] = os.path.splitext(os.path.basename(path))[0]
        match = re.search(regex_para[2], path_var['name'])
        if match.group(regex_para[0]) not in primers:
            path_var['check'] = False
            if verbose:
                print(f"pass: {path}")
        else:
            path_var['primer'] = match.group(regex_para[0])
            path_var['sample'] = match.group(regex_para[1])
            path_var['path'] = path
            path_var['check'] = True
    except:
        path_var['check'] = False
        if verbose:
            print(f"pass: {path}")
    return path_var

def checkRecordDict(record_dict: defaultdict, path_var: dict, verbose=True) -> bool:
    """
    判断在record_dict中path_var中保存的sample和primer是否存在。

    Args:
        record_dict: <defaultdict> 包含样本名和对应的序列记录列表的字典。
        path_var: <dict> 包含文件类型、文件名、样本名、引物名和'check'键的值。
            如果引物不在列表中，则'check'键的值为False。
            如果引物在列表中，则'check'键的值为True。
        verbose: <bool> 是否打印信息。

    Returns:
        <bool> 如果在record_dict中存在相同的sample和primer，则返回True，否则返回False。
    """
    if path_var['sample'] in record_dict:
        for record in record_dict[path_var['sample']]:
            if record.primer == path_var['primer']:
                if verbose:
                    print(f"pass: {path_var['path']}")
                return True
    return False

def checkFiles(paths, primers, regex_para) ->defaultdict:
    """
    根据文件名和待拼接的成对引物，生成一个字典存放导入的结果。

    Args:
        paths: <list> 文件名列表。
        primers: <list> 引物列表。
        regex_para:<list> primer所属group，sample所属的group，正则表达式。

    Returns:
        record_dict: <defaultdict> 包含样本名和对应的序列记录列表的字典。
    """ 
    record_dict = defaultdict(list)
    for path in paths:
        path_var = pathPaser(path, primers, regex_para)
        if path_var['check'] == True:
            record_dict[path_var['sample']] = []
    return record_dict

def showRecordDict(record_dict):
    """
    打印record_dict中每个样本的记录数、引物、文件类型和序列长度。

    Args:
        record_dict: <defaultdict> 包含样本名和对应的序列记录列表的字典。

    Returns:
        None
    """
    for sample in record_dict.keys():
        print(sample, len(record_dict[sample]))
        for record in record_dict[sample]:
            print("  ", record.primer, record.type, len(record.seq))

def checkPair(record_dict):
    """
    清除不成对的数据

    Args:
        record_dict: <defaultdict> 包含样本名和对应的序列记录列表的字典。

    Returns:
        None     
    """
    record_dict_check = {k: v for k, v in record_dict.items() if isinstance(v, list) and len(v) == 2}
    return record_dict_check

def phred2Index(record, check_len=0, thredhold=30, break_num=3):
    """
    将Phred转化为index

    Args:
        record: <SeqRecord> 序列记录。
        check_len: <int> 检查序列长度。
        thredhold: <int> 阈值。
        break_num: <int> 连续低于阈值的碱基数。

    Returns:
        index: <list> 包含起始和终止位置的列表。
    
    Note: //Fix 如果thredhold调整为20或40效果与预期不同 
    """
    index = [-1, -1]  # 如果是空的，后面赋值就会有问题，这里有-1+1=0也没有办法通过最后的长度检测
    # 根据阈值转变为boolean值
    Phred_list = record.letter_annotations["phred_quality"]
    if not Phred_list:
        return index   
    Phred_list_bool = [(x-thredhold) > 0 for x in Phred_list]
    # 找到大于阈值的起点
    for i, num in enumerate(Phred_list_bool):
        if num:
            index[0] = i
            break
    Phred_list = Phred_list[index[0]:]
    Phred_list_bool = Phred_list_bool[index[0]:]
    j=0
    for j, num in enumerate(Phred_list_bool):
        if j < (len(Phred_list)-1)/2: #过半之后开始找连续错误
            continue
        if not num: #出现了False
            if not any(Phred_list_bool[j:j+break_num]):
                index[1] = j+i
                break
    if not check_len:
        return index
    else:
        if index[1]-index[0] > check_len:
            return index
        else:
            return [index[0], -1]
 
def ab12Fas(record):
    """
    将ab1文件转化为fasta文件。

    Args:
        record: <SeqRecord> 序列记录。

    Returns:
        record_fas: <SeqRecord> fasta格式的序列记录。
    """
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
    """
    根据索引修剪fasta数据。

    Args:
        record: <SeqRecord> 序列记录。
        index: <list> 包含两个整数的列表，表示修剪的起始和结束索引。默认为[20, 750]。

    Returns:
        record_trim: <SeqRecord> 修剪后的fasta格式序列记录。
    """
    record_trim = record
    if len(record.seq) > index[1] > index[0]:
        record_trim.seq = record.seq[index[0]:index[1]]
        record_trim.description += f"trimmed {index[0]}:{index[1]}|"
        record_trim.index = index
    else:
        print(f"Seq length shorter than {index[1]}")
        record.index = [-1,-1]
    return record_trim

def swapSeq(seq1, seq2):
    '''
    交换两个序列的顺序。

    Args:
        seq1: <Seq> 第一个序列。
        seq2: <Seq> 第二个序列。

    Returns:
        seq1: <Seq> 第二个序列。
        seq2: <Seq> 第一个序列。
    '''

    seq = seq1
    seq1 = seq2
    seq2 = seq
    return seq1, seq2

def findLong(aligned_segs):
    """
    找到alignedSeges中最长的片段，它是从Bio.Align.PairwiseAlignment对象的aligned属性中获得的三层元组。
    返回目标序列和查询序列中最长片段的结束位置。

    Args:
        aligned_segs: <tuple> 从Bio.Align.PairwiseAlignment对象的aligned属性中获得的三层元组。

    Returns:
        marks: <list> 由两个整数组成的列表，表示目标序列和查询序列中最长片段的结束位置。
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

def mergeSeq(seq1, seq2, verbose=True):
    """
    输入两个Seq对象，反向互补第二个对象，然后将它们连接起来。
    返回一个包含序列、相似度、对齐长度和合并序列长度的列表。

    Args:
        seq1: <Seq> 第一个序列。
        seq2: <Seq> 第二个序列。
        verbose: <bool>是否输出拼接结果。

    Returns:
        res: <list> 包含合并序列、相似度、对齐长度和合并序列长度的列表。
    """


    # 反向互补
    seq2 = seq2.reverse_complement()
    # 构建比对
    aligner = Align.PairwiseAligner()
    # 将比对器设置为执行局部比对
    aligner.mode = 'local'
    aligner.match_score = 2
    aligner.mismatch_score = -2
    aligner.gap_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    # 进行比对
    alignments = aligner.align(seq1, seq2)

    # 解析结果
    if len(alignments) == 1:
        res = alignments[0]
        # 找到最长的片段
        marks = findLong(res.aligned)
        # 合并序列
        seqM = res.target[:marks[0]]+res.query[marks[1]:]
        # 计算相似度
        similarity = res.score/res.shape[1]/2*100
        # 计算对齐长度
        aligned_len = res.shape[1]
        # 返回结果
        return [seqM, similarity, aligned_len, len(seqM)]
    else:
        if verbose:
            print(f"   Failed {len(alignments)} alignments!")
        # 返回比对结果的数量
        return len(alignments)

def setRegex(path="", path_regex=""):
    """
    检测路径通过一个正在表达式能否提取到样品名和引物，返回primer，sample对应的group，以及通过测试的正则表达式
    Args:
        path: <str> 通过input()输入，需要解析的路径
        path_regex: <str> 用来解析的正则表达式
    Returns:
        regex_para: <list> [primer.index, sample.index, path_regex] 如果失败输出-1。
    """
    import os
    import re
    if not path:
        path=input("input path for parsing:")
    if not path_regex:
        path_regex= input("input path regex:")
    name_extension=os.path.splitext(os.path.basename(path))[-1]
    name=os.path.splitext(os.path.basename(path))[0]
    print(name)
    print("regex:",path_regex)
    match = re.search(path_regex,name)
    if not match:
        path_regex = input("Try another regex:")
        if path_regex == chr(27) or path_regex =="":
            return -1
        else:
            setRegex(path=path,path_regex=path_regex)
    else:
        for i in range(0,len(match.groups())+1): #match.groups() 和match.group()是不一样的
            print(f"Group {i} is {match.group(i)}")
        try:
            primer_index = int(input("Which group is primer? no fit = -1:")) 
            sample_index = int(input("Which group is sample? no fit = -1:")) 
        except:
            print("int is needed, -1 mean no fit and try new regex!")
            return -1
        if primer_index + sample_index>0:
            return [primer_index, sample_index, path_regex]
        else: 
            path_regex = input("Try another regex:")
            if path_regex == chr(27) or path_regex =="":
                return -1
            else:
                setRegex(path=path, path_regex=path_regex)


def main(folder="", primers=['27F', '1492R'], mod="auto", 
        check_len=0, thredhold=30, break_num=3, fastacut=[20, 750], 
        verbose=True, out="none", set_regex=False, regex_para=[]):
    """
    Functions:
        1. 输入文件，优先导入ab1，生成record_dict
        2. 清理序列字典，选择自动或手动将每条序列剪切到需要的范围生成merge_dict: {sample:[record1,record2],}
        3. 联配，遍历每个key，其中成对的序列进行拼接，拼接结果追加到当前的key中，得到merge_dict，选择导出序列
    Args:
        primers:<list>默认为27F和1492R
        mod:<str> auto or mannual
        check_len: <int> 剪切后短于这个数值的不要。
        threhold: <int> phred的阈值，目前20和40都不对
        break_num: <int>连续几个错误认为不再要
        fastacut:<list> fasta序列的默认剪切范围，默认20-750bp。
        verbose: <bool> 是不是要显示数据处理的全过程
        out:导出数据：none-不导出，single-一个文件一条序列，one-结果导出为1个文件
        set_regex: <bool>是否要设置文件名的解析正则表达式。
        regex_para: <list> primer所属group，sample所属的group，正则表达式。
    Returns:
        merge_dict: <dict> 数据结构是{sample:[record1,record2,record3],}。
    """
    # 1 输入阶段--------
    # 判断folder是不是一个文件夹，如果是，提取里面全部的文件名保存到paths
    if verbose:
        print("-----------input-----------")
    if folder=="":
        # 通过对话框输入   
        root = Tk()
        # 隐藏Tk窗口
        root.withdraw()
        # 打开文件对话框
        paths = filedialog.askopenfilenames()
        paths = list(paths)
    else:    
        if os.path.isdir(folder):
            paths = os.listdir(folder)
            paths = [os.path.join(folder, path) for path in paths]   
        else:
            print("Error: Input is not a folder.")
            return 
   
    # 文件排序 
    paths.sort(key=lambda x: (not x.endswith('ab1'), x.endswith('ab1') and x))
    if set_regex:
        #//Fix: 显示前10个导入的路径，选择其中1个进行解析，输出
        regex_para= setRegex()
        if regex_para == -1:
            print("Failed: setting path_regex")
            return
    # 解析文件名构建保存数据的字典
    if len(regex_para)==3 and isinstance(regex_para, list):
        record_dict = checkFiles(paths, primers, regex_para)
    else:
        print("need regex_para")
        return
    if len(record_dict)==0:
        print("no qualified files in the path.")
        return
    # 遍历文件导入
    for path in paths:
        # 解析文件名 如果无法解析到primer信息，就停止导入
        path_var = pathPaser(path, primers, regex_para)
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
                    ab1cut = phred2Index(record, check_len, thredhold, break_num)
                    merge_dict[sample].append(trimFasta(ab12Fas(record), index=ab1cut))
                else:
                    merge_dict[sample].append(trimFasta(record, index=fastacut))
        if verbose:
            print("-----------cut-----------")
            showRecordDict(merge_dict)  
    # 3 拼接阶段-----------
    if out =="one" and os.path.exists(os.path.join(folder,"output.fasta")):
        os.remove(os.path.join(folder,"output.fasta"))
        print("Warning: the output file is deleted!")

    for sample in merge_dict.keys():
        seqF = merge_dict[sample][0]
        seqR = merge_dict[sample][1]
        if "F" in seqR.primer and "R" in seqF.primer:
            seqF, seqR = swapSeq(seqF, seqR)
        seqM = mergeSeq(seqF.seq, seqR.seq)  # 需要转化为Seq对象
        if isinstance(seqM, int):
            continue
        description = f"len={str(seqM[3])}bp|{str(seqM[1])}% aligned of {str(seqM[2])}bp|{primers[0]} {primers[1]}"
        merged_seq = SeqRecord(seq=seqM[0], name=sample+"_[merge]", id=sample+"_[merge]",
                               description=description)
        merged_seq.primer = "merge"
        merged_seq.type = ".fasta"
        merge_dict[sample].append(merged_seq)
        # 拼接完成后是否输出summary
        if verbose:
            print("-----------merge-----------")
            showRecordDict(merge_dict)
        # 拼接完后是否导出数据
        if out == "none":
            return merge_dict
        elif out == "single":
            SeqIO.write(merged_seq, os.path.join(folder,merged_seq.name+merged_seq.type), "fasta")
        elif out == "one":
            with open(os.path.join(folder,"output.fasta"), "a") as output_file:
                SeqIO.write(merged_seq, output_file, "fasta")
    if verbose:
        print("\n")
    return merge_dict

if __name__ == "__main__":
    import argparse
    # 创建参数解析器
    parser = argparse.ArgumentParser(description='Process some integers.')
    # 添加参数
    parser.add_argument('--folder', type=str, default="", help='the folder containing seq files')    
    parser.add_argument('--primers', nargs='+', default=['27F', '1492R'], help='list of primer names')
    parser.add_argument('--mod', type=str, default='auto', help='mode of trimming')
    parser.add_argument('--check_len', type=int, default=0, help='minimum length of trimmed sequence')
    parser.add_argument('--thredhold', type=int, default=30, help='phred threshold')
    parser.add_argument('--break_num', type=int, default=3, help='number of consecutive errors to stop trimming')
    parser.add_argument('--fastacut', nargs='+', default=[20, 750], help='range of fasta sequence to trim')
    parser.add_argument('--verbose', type=bool, default=True, help='whether to display processing information')
    parser.add_argument('--out', type=str, default='none', help='output data: none- no output, single- one sequence per file, one- output all data to one file')
    parser.add_argument('--regex_para', nargs='+', default=[2, 1,"\\w{4}_(.*)_\\[(\\w+)\\]" ], help='list of primer names')
    parser.add_argument('--set_reget', type=bool, default=False, help='whether set regex for file name parsing')


    # 解析参数
    args = parser.parse_args()
    # 运行main函数
    # 测试时使用
    #args.folder = "E:\\Desktop\\python-learn-2023\\4-mergeSeq\\mergeSeq\\examples\\sequences"
    args.out = "single"
    #args.set_regex=True


    main(folder = args.folder, out=args.out, regex_para=args.regex_para)
