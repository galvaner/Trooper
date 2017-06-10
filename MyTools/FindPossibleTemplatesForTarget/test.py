import EmbossNeedle

from Bio.Emboss.Applications import NeedleCommandline

aseq = "acgtggcc"
bseq = "aaatttccggtt"
needle_cline = NeedleCommandline(asequence="asis:"+aseq, bsequence="asis:"+bseq, gapopen=10, gapextend=0.5, outfile="needle.txt")
stdout, stderr = needle_cline()
print stdout + stderr

with EmbossNeedle.MyEmboss('4L81_A', '4LVV_A') as MyEmbossInstance:
    print(MyEmbossInstance.GetSimilarity())
    print(MyEmbossInstance.GetGaps())
    print(MyEmbossInstance.GetAlignment())