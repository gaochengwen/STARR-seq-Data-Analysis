library(dplyr)
library(readr)
library(tibble)
library(purrr)
library(RColorBrewer)
library(ape)
library(reshape2)

low_expr_TPM = 0.1
low_expr_TPM_percent = 0.2
RLEFilterPercent = 0.05
DSFilterPercent = 0.05
topk_genes = 50
cluster_percent = 0.6
pvalue_cutoff = 0.05
cluster_level = 5

pvalues.cut <- pvalue_cutoff
treesNum <- cluster_level

argv <- commandArgs(T)

count_data <- read.table(paste0(argv[1],"_count.txt"),sep="\t")
TPM_data <- apply(count_data,2,function(x){x/sum(x)})*5000

#### Filter out low exp genes
keep_genes_idx <- (rowMeans(TPM_data>low_expr_TPM)>low_expr_TPM_percent) 
TPM_data = TPM_data[keep_genes_idx,]
message(paste(sum(1 - keep_genes_idx), "genes are filtered by filter: >", low_expr_TPM_percent*100, "% samples have expression values <", low_expr_TPM))
output <- paste0(argv[1],"_tpm_qc1.txt")
message(paste(sum(keep_genes_idx), "genes left, saving output to", output))
TPM_data%>%as_tibble(rownames = "gene_ID")%>%write_delim(output,"\t")

    DSFilter <- DSFilterPercent*ncol(TPM_data)

    ### RLE
    # https://github.com/stormlovetao/eQTLQC/blob/86dcc388c8da7f1bd5b223f4b9b26f09c907eb15/Sample/src/report.Rmd#L71
    log_offset <- 0.0001
    logtpm = log10(TPM_data%>%as.matrix + log_offset)
    rle=logtpm-apply(logtpm, 1, median) 
    iqr = apply(rle,2,IQR)
    rle=melt( rle , variable.name = "Sample",value.name ="TPM", id="feature")
    names(rle) <- c("feature","Sample","TPM")
    rle_IQR <- rle %>% group_by(Sample) %>% summarise(IQR = IQR(TPM))
    rle_IQR_range <- rle_IQR$IQR %>% range %>% abs() %>% max()
    rle_IQR_range <- 2*rle_IQR_range %>% ceiling()
    bymedian <- with(rle, reorder(Sample, TPM, IQR))  # sort by IQR
    pdf(file = paste0(argv[1],".RLEplot.pdf"))
    boxplot(TPM ~ bymedian, data=rle, outline=F, ylim = c(-rle_IQR_range, rle_IQR_range), las=2, boxwex=1, col='gray', cex.axis=0.3, main="RLE plot before QC", xlab="", ylab="Residual expression levels", frame=F)
    dev.off()

    ### hcluster  
    sampleDists <- 1 - cor(logtpm, method='spearman')
    hc <- hclust(as.dist(sampleDists), method = "complete")
    hcphy <- as.phylo(hc)

    pdf(file = paste0(argv[1],".preQC_cluster.pdf"))
    plot(hcphy, type = "unrooted", cex=.2, lab4ut='axial',underscore = T, main="Sample clustering before QC (Spearman - Cor.)")
    dev.off()

    ### D-s
    D = apply(1-sampleDists, 1, median)
    pdf(file = paste0(argv[1],".D_stat_hist.pdf"))
    hist(D, breaks=100, ylab="Number of samples", xlab="D-statistic", main="Histogram of Sample D-statistics before data QC")
    dev.off()

    #DSFilter <- sort(D)[DSFilter]
    D<-as.data.frame(D)
    D<-data.frame(Sample = rownames(D),D = D$D)
    D_filterList = D%>%filter(D <= DSFilter)
    D_filterList <- D_filterList$Sample
    D_filterList<-as.character(D_filterList)
    message(paste0("The right most ", DSFilterPercent*100, "% samples (N=", length(D_filterList), ") are marked as candidate outliers in this step.") )
    
    outliersList <- D_filterList
    message("Outliers:")
    outliersList

    outliersIndex <- which(colnames(logtpm) %in% outliersList)
    if(!length(outliersIndex) == 0){
        TPM_data <- TPM_data[,-outliersIndex]
    }

    output <- paste0(argv[1],"_tpm_qc.txt")
    TPM_data%>%as_tibble(rownames = "gene_ID")%>%write_delim(output,delim = "\t",col_names = T, append = T)

