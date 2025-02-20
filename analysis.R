library(tidyverse)
library(dplyr)
library(stringr)

argv <- commandArgs(T)
vars <- read.table(argv[1],sep="\t",head=T)
head(vars)
rownames(vars)=paste(vars$SNP,vars$Type,sep="_")
exps <- read.table(argv[2],sep="\t")
head(exps)

data = array(dim = 5)
n=dim(exps)[1]

for (i in 1:n) {
   tpm <- read.table(exps[i,2],sep="\t",head=T)
   rownames(tpm) <- tpm$gene_ID
   tpm = tpm[,-1]
   if(exps[i,3]+exps[i,4] != dim(tpm)[2]) { print(paste0("Error found in ",exps[i,2])) }
   tpm[tpm<as.double(argv[3])]=NA
   for(j in 1:exps[i,3])
   { 
     DNA = colnames(tpm)[j]
     print(DNA)
     for(k in 1:exps[i,4])
     {
       CELL = colnames(tpm)[exps[i,3]+k]
       print(CELL)
       for(m in 1:dim(tpm)[1])
       {
         Ratio=tpm[m,exps[i,3]+k]/tpm[m,j];
         data=rbind(data,c(rownames(tpm)[m],exps[i,1],DNA,CELL,Ratio))
       }
     }
   }
}

data = data[-1,]
colnames(data) = c("Variant","Exp","DNA","CELL","Ratio")
write.table(data, file = "ratios.txt", sep = "\t", row.names = FALSE, quote = FALSE)

data <- read.table("ratios.txt",sep="\t",head=T)
data$Ratio1 =log(data$Ratio)

pdf("ALL.QQ.pdf")
data1<- data[data$DNA=="DNA",]
ggplot(data1, aes(sample = Ratio1, colour = factor(CELL), shape = factor(Exp)))  +  stat_qq() +  stat_qq_line()
dev.off()

pdf("EXP1.QQ.pdf")
data2<- data[data$Exp=="Exp",]
ggplot(data2, aes(sample = Ratio1, colour = factor(CELL), shape = factor(DNA)))  +  stat_qq() +  stat_qq_line()
dev.off()

data$Ratio2=NA
rownames(data)=paste(data$Variant,data$Exp,data$DNA,data$CELL,sep=":")

n=dim(data)[1]
for (i in 1:n) {
    varname = rownames(data)[i]
    if(str_detect(varname, "Allele_1"))
    {
          varname2 = str_replace(varname, "Allele_1", "Allele_2")
          data[i,"Ratio2"]=log(data[i,"Ratio"]/data[varname2,"Ratio"])
    }
}

pdf("ALL.RR.QQ.pdf")
data1<- data[data$DNA=="DNA",]
ggplot(data1, aes(sample = Ratio2, colour = factor(CELL), shape = factor(Exp)))  +  stat_qq() +  stat_qq_line()
dev.off()

pdf("EXP1.RR.QQ.pdf")
data2<- data[data$Exp=="Exp",]
ggplot(data2, aes(sample = Ratio2, colour = factor(CELL), shape = factor(DNA)))  +  stat_qq() +  stat_qq_line()
dev.off()

data3 <- data[str_detect(rownames(data), "Allele_1"),]

snps = unique(data3$Variant)
exps = unique(paste(data3$Exp,data3$DNA,data3$CELL,sep=":"))
sn = length(snps)
en = length(exps)
rrdata = array(dim=c(sn,en))
rownames(rrdata)=str_replace(snps, "_Allele_1_",":")
colnames(rrdata)=exps

for(i in 1:sn)
{
     for(j in 1:en)
     {
         rrdata[i,j]=data3[paste0(snps[i],":",exps[j]),"Ratio2"]
      }
}

write.table(rrdata, file = "rrdata.txt", sep = "\t", quote = FALSE)
