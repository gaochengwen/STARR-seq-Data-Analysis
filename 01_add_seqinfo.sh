#/bin/sh

awk {' if(FNR==NR) { seqs[$3"_"$7]=$5; alleles[$3"_"$7]=$4; } else { if(seqs[$(NF-1)"_"$NF]) print $0,seqs[$(NF-1)"_"$NF],alleles[$(NF-1)"_"$NF] } '} ./variants.txt $1 > $1"_seqinfo.txt"
