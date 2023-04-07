from Bio.Seq import Seq

seq1 = Seq("ATC")
seq2 = Seq("TTT")
def swapSeq(seq1, seq2):
    # 序列顺序交换
    seq = seq1
    seq1 = seq2
    seq2 = seq
    return seq1, seq2

seq1, seq2 = swapSeq(seq1,seq2)
print(seq1, seq2)