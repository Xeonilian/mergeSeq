# 基于Bio.Align.PairwiseAlignment对象的aligned，为了方便测试，传递一个做的tuple
# alignedSegs = res.aligned
aligned_segs=(((592,593),(594, 772),(777,789)),((0,1),(2, 180),(185,196)))
# 应该返回[772,180]
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

findLong(aligned_segs)