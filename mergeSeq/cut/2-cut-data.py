from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

# 清除不成对的数据
def showRecordDict(record_dict):
    for sample in record_dict.keys():
        print(sample, len(record_dict[sample]))
        for record in record_dict[sample]:
            print("  ", record.primer, record.type, len(record.seq))
    return record
def checkPair(record_dict):
    record_dict_check = {k: v for k, v in record_dict.items(
    ) if isinstance(v, list) and len(v) == 2}
    return record_dict_check

# 将Phred转化为index
def Phred2Index(record, check_len=0, thredhold=30, break_num=3):
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

# ab1转化为fasta
def Ab12Fas(record):
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

# 基于index剪切fasta数据

def trimFasta(record, index=[20, 750]):
    record_trim = record
    if len(record.seq) > index[1] > index[0]:
        record_trim.seq = record.seq[index[0]:index[1]]
        record_trim.description += f"trimmed {index[0]}:{index[1]}|"
        record_trim.index = index
    else:
        record.index = [-1,-1]
    return record_trim
    

######### main
# 生成相同数据结构
def main(record_dict, mod="auto", check_len=0, thredhold=30, break_num=3, fastacut=[20, 750]):
    """
    清理序列字典，选择自动或手动将每条序列剪切到需要的范围
    Args:
        mod: auto or mannual
        check_len: fasta检查长度
        fastacut：fasta剪切的范围
        thredhold：Q值的阈值
        break_num: 连续多少个比阈值低
        
    Returns:
        merge_dict: {sample:[record1,record2],}
    """
    record_dict_clean = defaultdict(list)
    record_dict_clean = checkPair(record_dict)
    # showRecordDict(record_dict_clean) #测试通过
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
        showRecordDict(merge_dict)  # 测试通过2个
        return merge_dict

############ 运行
### 数据初始化
import pickle

dict_path = 'example-python310-record_dict.pickle'
with open(dict_path, 'rb') as file:
    record_dict = pickle.load(file)


import argparse

def Phred2Index(record, check_len=0, thredhold=30, break_num=3):
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

"""
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--check_len', type=int, help='an integer for the check length')
    parser.add_argument('--fastacut_start', type=int, help='an integer for the start of the fasta cut')
    parser.add_argument('--fastacut_end', type=int, help='an integer for the end of the fasta cut')
    parser.add_argument('--break_num', type=int, help='an integer for the break number')
    parser.add_argument('--thredhold', type=int, help='an integer for the Q-value')
    args = parser.parse_args()
    if args.thredhold is None:
        args.thredhold = 30
    main(record_dict, check_len=args.check_len, fastacut=[args.fastacut_start, args.fastacut_end], break_num=args.break_num, thredhold=args.thredhold)
"""

main(record_dict, check_len=600,fastacut=[1,900],break_num=3, thredhold=20) 