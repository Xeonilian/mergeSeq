from Bio import Align

# merge的主函数
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
        return [seqM, similarity, aligned_len,len(seqM)]