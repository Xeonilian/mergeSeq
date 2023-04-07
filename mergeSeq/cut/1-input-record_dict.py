import os
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from Bio.Seq import Seq


###################### 函数定义 #########################

def ab1IO(path_var: list):
    # 导入ab1文件
    for record in SeqIO.parse(path_var['file_name'], "abi"):
        record.primer = path_var['primer']
        record.sample = path_var['sample']
        record.type = path_var['file_type']
        record.description += f"{record.type}|"
        record_dict[path_var['sample']].append(record)
        print("  import: ab1")

def fasIO(path_var: list):
    # 导入fasta文件
    try:  # 试了试文件对不对
        next(SeqIO.parse(path_var['file_name'], "fasta"))
        for record in SeqIO.parse(path_var['file_name'], "fasta"):  # 如果是fasta文件
            record.primer = path_var['primer']
            record.sample = path_var['sample']
            record.type = path_var['file_type']
            record.description += f"{record.type}|"
            print("  import fasta")
            record_dict[path_var['sample']].append(record)
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
        print("  import seq")
        record_dict[path_var['sample']].append(record)

def pathPaser(file_name: str, primers:list):
    # 解析文件名
    path_var = {"a":""}
    try:
        path_var['file_type'] = os.path.splitext(os.path.basename(file_name))[-1]
        path_var['name'] = os.path.splitext(os.path.basename(file_name))[0]
        match = re.search(r"(.*)_\[(\w+)\]", path_var['name'])
        if match.group(2) not in primers:
            path_var = {"a":""}
            path_var['check'] = False
        else:
            path_var['primer'] = match.group(2)
            path_var['sample'] = match.group(1)
            path_var['file_name'] = file_name
            path_var['check'] = True
    except:
        path_var = {"a":""}
        path_var['check'] = False
    return path_var

def checkRecordDict(record_dict, path_var: list):
    # 判断在record_dict中path_var中保持的样品和primer是不是存在
    if path_var['sample'] in record_dict:
        for record in record_dict[path_var['sample']]:
            if record.primer == path_var['primer']:
                return True
    return False

def checkFiles(file_names, primers):
    # 生成一个字典存放导入的结果
    record_dict = defaultdict(list)
    # print(type(record_dict))
    for file_name in file_names:
        # print(file_name)
        path_var = pathPaser(file_name, primers)
        if path_var['check'] == True:
            record_dict[path_var['sample']] = []
    return record_dict

def showRecordDict(record_dict):
    for sample in record_dict.keys():
        print(sample, len(record_dict[sample]))
        for record in record_dict[sample]:
            print(f"  {record.primer} {record.type} {len(record.seq)}bp")
        
######################## 运行
### 初始化阶段
# 文件路径输入
file_names = ['E:/Desktop/python-learn-2023/4-mergeSeq/mergeSeq/data-example/0005_31320101600206_(5)_[27F].ab1',
              'E:/Desktop/python-learn-2023/4-mergeSeq/mergeSeq/data-example/0005_31320101300621_(5)_[1492R].seq',
              'E:/Desktop/python-learn-2023/4-mergeSeq/mergeSeq/data-example/0005_31320101300621_(5)_[1492R].fasta',
              'E:/Desktop/python-learn-2023/4-mergeSeq/mergeSeq/data-example/0005_31320101300621_(5)_[1492R].ab1',
              'E:/Desktop/python-learn-2023/4-mergeSeq/mergeSeq/data-example/0005_31320101300621_(5)_[27F].ab1']
file_names.sort(key=lambda x: (not x.endswith('ab1'), x.endswith('ab1') and x))
# 设置primer对
primers = ['27F', '1492R']
# 解析文件名构建保存数据的字典
record_dict = checkFiles(file_names, primers)

###  运行阶段
# 遍历文件导入
for file_name in file_names:
    # 解析文件名 如果无法解析到primer信息，就停止导入
    path_var = pathPaser(file_name, primers)
    if path_var['check'] == True:
        # 判断该文件是否存在
        if not checkRecordDict(record_dict, path_var):
            if path_var['file_type'] == ".ab1":  
                ab1IO(path_var)
            else:
                fasIO(path_var)

### 输出阶段 
# 说明导入了2个样品，其中一个只有一个数据
for sample in record_dict.keys():
    print(sample, len(record_dict[sample]))
    for record in record_dict[sample]:
        print("  ", record.primer)

import pickle
with open('example-python310-record_dict.pickle', 'wb') as file:
    pickle.dump(record_dict, file)

