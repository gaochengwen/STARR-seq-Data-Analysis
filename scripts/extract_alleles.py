import sys
import numpy as np
from Bio import pairwise2

seqcounts = sys.argv[1]
datafromfile = np.loadtxt(seqcounts, dtype="str")
seqn = datafromfile.shape[0]

for i in range(seqn):
   data = datafromfile[i]
   seq1 = data[1]
   seq2 = data[data.size-2]
   posinfo = data[data.size-3]
   pos = int(posinfo.split("_")[3])
   alleles = data[data.size-1]
   allele2 = alleles.split("/")[1]
   alignments = pairwise2.align.globalxx(seq1, seq2)
   A = alignments[0].seqA
   B = alignments[0].seqB
   #if len(allele2) > 1: 
   #print(alleles,seq2[pos-len(allele2):pos-200+len(seq2)])
   startp = B.find("-",0)
   while startp != -1:
      if startp == 0:
         A = A[1:]
         B = B[1:]
         startp = B.find("-",startp)
      else:
         A = A[:startp] + A[startp+1:]
         B = B[:startp] + B[startp+1:]
         startp = B.find("-",startp-1)
   print(data[0],alignments[0].score/200,alleles,seq2[pos-len(allele2):pos-200+len(seq2)],B[pos-len(allele2):pos-200+len(seq2)],A[pos-len(allele2):pos-200+len(seq2)])
