import argparse
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

"""从ab1文件中，找到第二峰对应的序列
    前题：峰高满足要求
    重叠峰，但是有效识别。点突变，移码涂布，高相似度序列的污染。找到第二个模板对应的序列。
"""
def checkQuality(record, cutoff, check_len,verbose):
    """
    Richard Mott修建算法，找到符合正确率的范围，但是将cutoff设置为0.1，来自biopython，原来是0.05也就是95%的正确率
    http://www.phrap.org/phredphrap/phred.html
    http://www.clcbio.com/manual/genomics/Quality_trimming.html
    返回的index是实际的位置
    Returns:
        index: <tuple> 符合的质量标准的范围
    """
    if len(record.seq)<check_len:
        raise ValueError("Seqeunce is too short!")
    if any([num < 100 for num in record.annotations['abif_raw']['S/N%1']]):
        raise ValueError("S/N is too low!")
    else:
        # 设置修剪标志
        start = False
        trim_start = 0
        # 从用于计算phred qual值的公式中计算概率回来
        score_list = [cutoff - (10 ** (qual/-10.0)) for 
                        qual in record.letter_annotations["phred_quality"]]

        # 计算累积score_list
        # 如果累积值<0，则设置为0
        # 第一个值设置为0（假设：trim_start始终>0）
        running_sum = [0]
        for i in range(1, len(score_list)):
            num = running_sum[-1] + score_list[i]
            if num < 0:
                running_sum.append(0)
            else:
                running_sum.append(num)
                if not start:
                    # trim_start = value when cummulative starts to be > 0
                    trim_start = i
                    start = True

        # trim_finish = index of the highest cummulative value,
        # marking the segment with the highest cummulative score 
        trim_finish = running_sum.index(max(running_sum)) 
        percentage = 100*(1-cutoff) 
        index = (trim_start+1,trim_finish+1)
        if verbose:
            print(f"Correction rate at {percentage}% from {index[0]} to {index[1]}")
        return index
        # return record[trim_start:trim_finish]
def checkSecond(record, index, cutoff2, brk, check_len, verbose):
    """找到第二条序列的范围
    """
    h1=record.annotations['abif_raw']['P1AM1'][index[0]-1:index[1]-1]
    h2=record.annotations['abif_raw']['P2AM1'][index[0]-1:index[1]-1]
    # 设置修剪标志
    cutoff2=0.1
    brk = '0'* brk
    height_list=[]
    for m,n in zip(h1, h2):
        if (n/m)-cutoff2 >0 :
            height_list.append(int(((n/m)-cutoff2)*10))
        else:
            height_list.append(0)

    height_list = "".join(map(str,height_list))
    segments = [segment for segment in height_list.split(brk) if len(segment)>check_len]
    breaks={}
    # 找到每个片段的起始位置
    for i, segment in enumerate(segments):
        start_zeros = len(segment) - len(segment.lstrip('0'))
        end_zeros = len(segment) - len(segment.rstrip('0'))
        height_list.index(segment) 
        breaks[i]=(height_list.index(segment)+start_zeros+index[0], 
                    height_list.index(segment)+len(segment)-end_zeros+index[0])
    return breaks
def selectBreaks(record):
    """
    打开指定范围的序列图，通过点击设置截取识别的范围
    """
    print("in dev! set the range mannually")
    string = input("Input start and end (eg: 20,100;120,300):")
    breaks={}
    for i, s in enumerate(string.split(';')):
        if s == "":
            continue
        else:
            breaks[i] = tuple(map(int, s.split(',')))
    return breaks
def pickPeak(record, start_end,cutoff2):
    """根据图形选择解析的起始位置
    
    """
    seq3=Seq("")
    seq1=record.seq[start_end[0]-1:start_end[1]-1]
    seq2=Seq(record.annotations['abif_raw']['P2BA1'].decode('utf-8')[start_end[0]-1:start_end[1]-1])
    h1=record.annotations['abif_raw']['P1AM1'][start_end[0]-1:start_end[1]-1]
    h2=record.annotations['abif_raw']['P2AM1'][start_end[0]-1:start_end[1]-1]
    for i in range(len(h1)):
        if h2[i]/h1[i]>cutoff2: #第二条序列的高度大于第一条序列的比例cutoff2
            seq3+=seq2[i]
        else:
            seq3+=seq1[i]

    if len(record.seq) > start_end[1]:
        # 从record.annotations中获取abif_raw中的P2BA1序列
        second_seq = SeqRecord(seq3, id=record.id+f"_second[{start_end[0]}:{start_end[1]}]", 
        name=record.name, description=record.description)
        return second_seq

    else:
        raise ValueError("sequence is shorter than index!")
def main():
    """解析第二序列

    """
    # 创建参数解析器
    parser = argparse.ArgumentParser()
    # 添加参数
    parser.add_argument('--ab1_path', type=str, default="", help='the path used to store ab1 file')    
    parser.add_argument('--out', type=str, default="second_seqs.fasta", help='the output file name')
    parser.add_argument('--mod', type=str, default='default', choices=["auto","mannual","default"], 
                        help='mode of sequence cutting')
    parser.add_argument('--check_len', type=int, default=20, help='minimum length of trimmed sequence')
    parser.add_argument('--cutoff', type=float, default=0.05, help='error probability less than 5%')
    parser.add_argument('--cutoff2', type=float, default=0.25, help='error probability less than 5%')
    parser.add_argument('--brk', type=int, default=2, help='error probability less than 5%')
    parser.add_argument('--verbose', type=bool, default=True, help='whether to display processing information')
    parser.add_argument('--test', type=bool, default=False, help='do test')
    # 解析参数
    args = parser.parse_args()
    # 测试用 
    args.test = True
    if args.test:
        args.ab1_path="E:\\Desktop\\python-learn-2023\\4-mergeSeq\\mergeSeq\\data\\0030_31323021000779_(16)_[1492R].ab1"
        args.mod="mannual"
        args.check_len=20
        args.cutoff = 0.1
    # 检查文件是否存在
    if os.path.isfile(args.ab1_path):
        # 解析ABI文件
        try:
            record = next(SeqIO.parse(args.ab1_path, "abi"))
        except StopIteration:
            raise ValueError("No records found in ABI file")

        # 根据模式选择序列切割方式获得breaks字典
        if args.mod == "auto":
            index = checkQuality(record, args.cutoff, args.check_len, args.verbose)
            breaks= checkSecond(record, index, args.cutoff2, args.brk, args.check_len, args.verbose)
        elif args.mod == "mannual":
            # 手工输入
            breaks = selectBreaks(record)          
        # 检查序列切割是否成功

        for key in breaks.keys():
            second_seq = pickPeak(record, breaks[key], args.cutoff2)
            # 将第二条序列写入文件
            if args.out and second_seq:
                with open(os.path.join(os.path.dirname(args.ab1_path), args.out), "a") as output_file:
                    SeqIO.write(second_seq, output_file, "fasta")
    else:
        # 如果文件不存在，输出错误信息
        if args.verbose:
            print("Wrong Path!")
        sys.exit(1)
if __name__ == "__main__":  
    main()
