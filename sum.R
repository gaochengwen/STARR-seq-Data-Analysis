library(dplyr)

argv <- commandArgs(T)
df1=read.table(paste0(argv[1],"_qc.txt"),sep=" ")
expname = strsplit(readLines(argv[1], n=1), split = " ")[[1]]
dn = length(expname)-1

sumdata = aggregate(df1[,3:dn], by=list(df1[,dn+2],df1[,dn+1]), FUN=sum)
colnames(sumdata) <- c("Alleles", "SNP", expname[3:dn-1])

sumdata[1:10,]

write.table(sumdata, file = paste0(argv[1],"_sum.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
