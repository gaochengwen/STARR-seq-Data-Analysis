awk {' OFS="\t"; $2=$2"_"$1; $1=""; if(FNR==1) { $2="" } print $0 '} $1"_sum.txt" | sed 's/^[\t]*//g' > $1"_count.txt"
