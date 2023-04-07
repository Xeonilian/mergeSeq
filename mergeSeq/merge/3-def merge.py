from Bio import Align
from Bio.Seq import Seq

def mergeSeq(seqF,seqR):
    #序列顺序交换
    if "R" in seqF.primer and "F" in seqR.primer:
        seq = seqF
        seqF= seqR
        seqR=seq
    aligner = Align.PairwiseAligner()

    # Set the aligner to perform local alignments
    aligner.mode = 'local'
    aligner.match_score = 2 
    aligner.mismatch_score = -2 
    aligner.gap_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5 
    alignments = aligner.align(seqF, seqR)  
    if len(alignments)==1:
        # 提取第一个结果
        res = alignments[0]
        
        if res.aligned[0]==2:
            merged=res.target[:res.aligned[0][0][1]]+res.query[res.aligned[1][0][1]:] #targe中第0组第二个数字，query中第0组的第二个数字
            print("Similarity:100%")

                
    else:
        print("multiple resutls!")