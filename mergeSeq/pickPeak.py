from Bio import Align
from mergeSeq.merge import phred2Index
import argparse
import os
from Bio import SeqIO
#// TODO:要保证能正确的import

"""
    目标：从ab1文件中，找到第二峰对应的序列
        找到污染序列（突变序列）：两条序列有很高的相似度，一条的原始模板含量较少，在高质量范围，本来一个位置只能检测到一个主导的信号，但是污染序列，也产生了一个清晰的信号，出现了重叠峰。
    问题：有的虽然解析了第二峰，但是第二峰位移严重，高度不足以达到认定为是一个明确的峰，在某些位置不存在第二峰，
    方案：能够按照一定标准识别清晰的第二峰，在峰图上可视化，然后再对区域进行选择，复制这个区域，当第二峰不足一定标准的时候，就用第一峰替代这些位置

"""
def checkQuality(record):
    """
    检测指定范围是否是正常测序结果，Phred>30，允许有断裂不超过x个，
    """
    # 根据阈值转变为boolean值
    Phred_list = record.letter_annotations["phred_quality"]
    
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

    return True
def selectIndex(record):
    """
    在指定范围检测是否存在第二条序列，第二条序列应该高度大于信号的10%
    """
    print("in dev! set --mod as default or auto")
    pass

def pickPeak(record, index, fragment):
    """
    根据图形选择解析的起始位置
    """
    # 判断是一条正常的测序结果
    second_seqs = []
    if second_seqs:
        return second_seqs
    else:
        return -1

def main():
    """

    """
    # 创建参数解析器
    parser = argparse.ArgumentParser()
    # 添加参数
    parser.add_argument('--ab1_path', type=str, default="", help='the path used to store ab1 file')    
    parser.add_argument('--out', type=str, default="second_seqs.fasta", help='the output file name')
    parser.add_argument('--mod', type=str, default='default', choices=["auto","mannual","default"], 
                        help='mode of sequence cutting')
    parser.add_argument('--check_len', type=int, default=0, help='minimum length of trimmed sequence')
    parser.add_argument('--thredhold', type=int, default=30, help='phred threshold')
    parser.add_argument('--break_num', type=int, default=3, help='number of consecutive errors to stop trimming')
    parser.add_argument('--default_index', nargs='+', default=[20, 700], help='default selected range')
    parser.add_argument('--fragment', type=bool, default=False, help='whether give mutilple fragments or one sequences')
    parser.add_argument('--verbos', type=bool, default=True, help='whether to display processing information')
    parser.add_argument('--test', type=bool, default=False, help='do test')
    # 解析参数
    args = parser.parse_args()
    # 测试用 
    args.test = True
    if args.test:
        args.ab1_path="E:\\Desktop\\python-learn-2023\\4-mergeSeq\\mergeSeq\\examples\\sequences\\0001_31319091700644_(GN15-15)_[27F].ab1"
        args.mod="auto"
        args.check_len=10
        args.thredhold=20
        args.break_num=3
        args.fragment=True
# 检查文件是否存在
    if os.path.isfile(args.ab1_path):
        # 解析ABI文件
        try:
            record = next(SeqIO.parse(args.ab1_path, "abi"))
        except StopIteration:
            print("Error: No records found in ABI file")
            sys.exit(1)
        # 检查序列质量
        if checkQuality(record):
            # 根据模式选择序列切割方式
            if args.mod == "auto":
                index = phred2Index(record, args.check_len, args.thredhold, args.break_num) # 函数在merge模块中
                if args.verbos:
                    print("Select the range based on phred value!")
            elif args.mod == "default":
                index = args.default_index
                if len(record.seq) < index[1]:
                    index[1] = len(record.seq)
                if args.verbos:
                    print("Select the range based on default setting")
            elif args.mod == "mannual":
                index = selectIndex(record)
            # 检查序列切割是否成功
            if index:
                # 按照指定范围检测是否存在第二条序列
                second_seqs = pickPeak(record, index, args.fragment)
                # 将第二条序列写入文件
                if args.out and second_seqs:
                    with open(os.path.join(os.path.dirname(args.ab1_path), args.out), "a") as output_file:
                        SeqIO.write(second_seqs, output_file, "fasta")
                        # 根据verbos参数输出信息
                        if args.verbos:
                            print("Second seqeucnes are exported!")
                return second_seqs  
            else:
                # 如果切割失败，根据verbos参数输出错误信息
                if args.verbos:
                    print("Please give an index!")
                sys.exit(1) 
        else:
            # 如果序列质量过低，输出错误信息
            if args.verbos:
                print("The quaility of sequence is too low to parse!")
            sys.exit(1)
    else:
        # 如果文件不存在，输出错误信息
        if args.verbos:
            print("Wrong Path!")
        sys.exit(1)

if __name__ == "__main__":  
    main()
