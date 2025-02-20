#/bin/sh

awk {' if(FNR==NR) { if($2>0.90 && $6==$4) { a[$1]=1; info[$1]=$0; } } else { if(a[$1]) { print $0,info[$1] } } '} $1"_alleles.txt" $1"_seqinfo.txt" > $1"_qc.txt"
