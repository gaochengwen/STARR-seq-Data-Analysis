# Load necessary library
library(mpra)

# Set the working directory to where the input files are located
# Please change this path according to your system
setwd('C:/path/to/your/data/')

# Function to safely read CSV files and handle errors
read_data <- function(file_path) {
  tryCatch({
    data <- read.csv(file_path, stringsAsFactors = FALSE)
    return(data)
  }, error = function(e) {
    stop(paste("Error reading file:", file_path, "\n", e))
  })
}

# Read DNA and RNA data
dna <- as.matrix(read_data("dna.csv"))
rna <- as.matrix(read_data("rna.csv"))

# Read eSeq data
eseq <- as.matrix(read_data("eseq.csv"))

# Read eid data (assuming the first column contains eid information)
eid_data <- read_data("eid.csv")
eid <- as.character(eid_data[, 1])

# Check if the dimensions of DNA and RNA matrices match
if (ncol(dna) != ncol(rna)) {
  stop("Error: The number of columns in the DNA matrix does not match the RNA matrix.")
}

# Create an MPRASet object with the data
mpraset_allelic_comparison <- tryCatch({
  MPRASet(DNA = dna, RNA = rna, eid = eid, eseq = eseq, barcode = NULL)
}, error = function(e) {
  stop("Error creating MPRASet object: ", e)
})

# Define the design matrix for the analysis
# 'intcpt' is an intercept term, 'allele2' checks for allele2 in the column names
design <- data.frame(intcpt = 1, allele2 = grepl("allele2", colnames(mpraset_allelic_comparison)))

# Create a block vector for grouping the data
# Here, we assume the data is grouped into three blocks
block_vector <- rep(1:3, length.out = ncol(mpraset_allelic_comparison))

# Perform the MPRA analysis
mpralm_allele_fit <- tryCatch({
  mpralm(object = mpraset_allelic_comparison, design = design, aggregate = "none", normalize = TRUE, 
         block = block_vector, model_type = "corr_groups", plot = TRUE)
}, error = function(e) {
  stop("Error performing MPRA analysis: ", e)
})

# Extract the results for the 'allele2' coefficient
# Display the top results
toptab_allele <- topTable(mpralm_allele_fit, coef = 2, number = Inf)
head(toptab_allele)

# Save the results to a CSV file
write.csv(toptab_allele, "toptab_allele.csv")

# Print a message indicating that the script has completed successfully
cat("MPRA analysis completed successfully. Results saved in 'toptab_allele.csv'.\n")
