module load anaconda/Anaconda3_2021 && python extract_alleles.py $1"_seqinfo.txt" | awk {' print $0 '} > $1"_alleles.txt"
