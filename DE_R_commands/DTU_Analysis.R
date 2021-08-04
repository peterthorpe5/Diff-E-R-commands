getwd()

# Step 0: install the rnaseqDTU Bioconductor package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.13")

# Step 1: load the sample CSV file
# The CSV file columns are labeled with the sample IDs so they can be passed
# to DRIMSeq later
samples <- read.csv(file.path("/users/pa/bioinformatics/DTU_Analysis",
                              "samples.csv"))
samples

# Convert the conditions vector from the samples file to a factor data type
samples$condition <- factor(samples$condition)
table(samples$condition)

# Obtain the file path to each of the samples
files <- file.path("/users/pa/bioinformatics/DTU_Analysis/salmon_transcript_quant",
                   samples$sample_id, "quant.sf")
names(files) <- samples$sample_id
head (files)

# Step 2: construct a matrix of counts from the quantification files
BiocManager::install("tximport")
library(tximport)
txi <- tximport(files, type = "salmon", txOut = TRUE,
                countsFromAbundance = "scaledTPM")
cts <- txi$counts
cts <- cts[rowSums(cts) > 0,]
head(cts)

# Create a TxDb object (transcript database) which links transcripts to
# their parent genes
BiocManager::install("GenomicFeatures")
library(GenomicFeatures)
st.dir        <- "/users/pa/bioinformatics/DTU_Analysis/"
gtf           <- paste(st.dir, sep = "", "Mus_musculus.GRCm38.84.gtf.gz")
txdb.filename <- paste(st.dir, sep = "", "Mus_musculus.GRCm38.84.sqlite")
txdb <- makeTxDbFromGFF(gtf)
saveDb (txdb, txdb.filename)
# Note: this workflow creates the .sqlite file

# Use the select function on the TxDb to produce a data frame which matches
# transcript IDs to gene IDs
# TXNAME refers to transcript IDs & GENEID refers to gene IDs
txdb <- loadDb(txdb.filename)      # Reload the TxDb database
txdf <- select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID")
tab  <- table (txdf$GENEID)
txdf$ntx <- tab[match(txdf$GENEID, names(tab))]

# # # LEFT UP TO HERE # # #

# Create a data frame with the gene ID, transcript (feature) ID, and columns
# for (counts from) each of the samples
all(rownames(cts) %in% txdf$TXNAME)  # %in% indicates whether TXNAME matches
# a row name in counts
txdf   <- txdf[match(rownames(cts), txdf$TXNAME),]
all(rownames(cts) == txdf$TXNAME)
counts <- data.frame(gene_id = txdf$GENEID, feature_id = txdf$TXNAME, cts)
head(counts)

# Step 4: load the DRIMSeq package and create a dmDSdata object with the
# counts and samples data frames
BiocManager::install("DRIMSeq")
library(DRIMSeq)
dmDS <- dmDSdata(counts = counts, samples = samples)
dmDS # returns information about the number of genes

# Each row of the dmDSdata object contains all the transcripts corresponding
# to a particular gene
# Retrieving the first row therefore corresponds to all the transcripts from
# the first gene
methods(class = class(dmDS))
counts(dmDS[1,])[,1:4]

# Step 5: filter the dmDS object before estimating model parameters in order
# to:
#   A) Accelerate the fitting of proportions by DRIMSeq
#   B) Remove potentially troublesome transcripts (for instance, those which
#      have very low abundance)

n       <- 12 # the total number of samples
n.small <- 6  # sample size of the smallest group

# Filters used (inclusion criteria for transcripts):
#   1) Must have a count of at least 10 in at least n.small samples
#   2) Must have a relative abundance of at least 10% of a gene's expression
#   3) The total count of the corresponding gene must be at least 10 in each
#      of the n samples
dmDS <- dmFilter(dmDS,
                 min_samps_feature_expr = n.small, min_feature_expr = 10 ,
                 min_samps_feature_prop = n.small, min_feature_prop = 0.1,
                 min_samps_gene_expr    = n,       min_gene_expr    = 10)
dmDS

# Find out how many of the genes remaining after filtering have N isoforms
# by counting the number of unique gene IDs and tabulating the results
table(table(counts(dmDS)$gene_id))
# After running the line above, note that the dmDSdata object has no genes
# with only a single isoform, as it would be impossible to test for DTU in
# such a scenario

# Step 6: create a design matrix using a design formula and the sample
# information contained in the dmDS object (accessed via samples.csv)
design_full <- model.matrix(~condition, data = DRIMSeq::samples(dmDS))
colnames(design_full)

# Step 7: estimate the model parameters, test for differential transcript use
# Estimate the precision, a parameter inversely related to dispersion (counts
# varying more from their expected values) in the Dirichlet Multinomial model
# dispersion = 1 / (1 + precision)
# Subsequently, fit regression coefficients & perform null hypothesis testing
# on the coefficient of interest (that associated with the difference between
# condition 1 and condition 2 in this simple two-group model ("condition2"))
set.seed(1)
system.time({
  dmDS <- dmPrecision(dmDS, design = design_full )
  dmDS <- dmFit      (dmDS, design = design_full )
  dmDS <- dmTest     (dmDS, coef   = "condition2")
  # error with condition2: "'coef' does not match
  # the columns of the design matrix!'"
})

# Tabulate the results, including a p-value per gene or a p-value per transcript
#   p-value per gene: is there DTU within this gene?
#   p-value per transcript: has the proportion of this transcript changed within
#   its parent gene?
results <- DRIMSeq::results(dmDS)                               # per gene
head(results)
results.txp <- DRIMSeq::results(dmDS, level = "feature")        # per transcript
head(results.txp)

# Step 8: replace NA values in place of p-values with 1
# The NA values would otherwise cause problems for stage-wise analysis
# NA p-values occur when one condition group has all zero counts for a transcript,
# but the corresponding gene and other condition group both have sufficient counts
# DRIMSeq doesn't estimate a precision for these genes
no.na <- function(x) ifelse(is.na(x), 1, x)
results$pvalue     <- no.na(results$pvalue)
results.txp$pvalue <- no.na(results.txp$pvalue)

# Step 9: plot the estimated proportions for a significant gene showing evidence
# of switching (DTU)
index <- which(results$adj_pvalue < 0.05)[1]
results[index,]
plotProportions(dmDS, results$gene_id[index], "condition")     # ! NOT WORKING !
# Error in `levels<-`(`*tmp*`, value = as.character(levels)) :
#   factor level [2] is duplicated

# Step 10: load results tables (these would have been generated using the whole
# data set as opposed to just a subset, as was done in this workflow)
data(drim_tables)
nrow(results)     # returns 250 (too few; should be 7764 )
nrow(results.txp) # returns 648 (too few; should be 20711)

# Step 11: detection of DTU using stageR, which enables two-stage testing
#   Stage 1: which genes contain any evidence of DTU?
#   Stage 2: which transcripts from these genes (the subset from stage 1)
#            are participating in the DTU? 

# Construct a vector of p-values for the screening (first) stage
pScreen <- results$pvalue
# Remove gene & transcript version numbers from their Ensembl IDs so that
# stageR will be able to combine gene & transcript names
remove  <- function(x) substr(x, 1, 15)
names(pScreen) <- remove(results$gene_id)

# Construct a one-column matrix of the confirmation p-values
pConfirmation <- matrix(results.txp$pvalue, ncol = 1)
rownames(pConfirmation) <- remove(results.txp$feature_id)

# Arrange a two-column data frame with the gene & transcript IDs
tx2gene <- results.txp[,c("feature_id", "gene_id")]
for (i in 1:2) tx2gene[,i] <- remove(tx2gene[,i])

# Specify an alpha (overall FDR target for the analysis)
library(stageR)
stageR_obj <- stageRTx(pScreen = pScreen, pConfirmation = pConfirmation,
                       pScreenAdjusted = FALSE, tx2gene = tx2gene)
stageR_obj <- stageWiseAdjustment(stageR_obj, method = "dtu", alpha = 0.05)
suppressWarnings({
  drim.padj <- getAdjustedPValues(stageR_obj, order    = FALSE,
                                  onlySignificantGenes = TRUE)
})
head(drim.padj)

# Step 12: post-hoc filtering to improve the FDR & overall FDR control using the
# standard deviation (SD) of per-sample proportions as a filtering statistic
# This procedure is required as DRIMSeq can exceed its FDR bounds, especially on
# the transcript-level tests
# Set the p-values and adjusted p-values for transcripts with a small per-sample
# proportion SD to 1
results.txp.filter <- DRIMSeq::results(dmDS, level = "feature")
smallProportionSD  <- function(dmDS, filter = 0.1) {
  cts       <- as.matrix(subset(counts(dmDS), select = -c(gene_id, feature_id)))
  gene.cts  <- rowsum(cts, counts(dmDS)$gene_id)
  total.cts <- gene.cts[match(counts(dmDS)$gene_id, rownames(gene.cts)),]
  props     <- cts / total.cts
  propsSD   <- sqrt(rowVars(props))
  propsSD   < filter
}
filter <- smallProportionSD(dmDS)
results.txp.filter$pvalue[filter]     <- 1
results.txp.filter$adj_pvalue[filter] <- 1

# Step 13: load the DEXSeq package and create a DEXSeqDataSet object from the
# data in the dmDStest object (the dmDSdata object after results are added)
library(DEXSeq)
sample.data <- DRIMSeq::samples(dmDS)
count.data  <- round(as.matrix(counts(dmDS)[,-c(1:2)]))
dxd <- DEXSeqDataSet(countData = count.data, sampleData = sample.data,
                     design    = ~sample + exon + condition:exon,
                     featureID = counts(dmDS)$feature_id,
                     groupID   = counts(dmDS)$gene_id)

# Run the DEXSeq analysis
system.time({
  dxd <- estimateSizeFactors(dxd)
  dxd <- estimateDispersions(dxd, quiet = TRUE)
  dxd <- testForDEU(dxd, reducedModel   = ~sample + exon)
})

# Step 14: extract the results table and compute a per-gene adjusted p-value
# The adjusted p-value is obtained as follows:
#    1) Aggregate the results of multiple tests within a gene into a p-value
#       for that gene
#    2) Correct the p-values for testing across several genes
dxr   <- DEXSeqResults(dxd, independentFiltering = FALSE)
qval  <- perGeneQValue(dxr)
dxr.g <- data.frame(gene = names(qval), qval)

# Reduce the transcript-level results table to a simple data frame for size
# considerations
col <- c("featureID", "groupID", "pvalue")
dxr <- as.data.frame(dxr[,col])
head(dxr)

# Load results tables (these would have been generated using the whole data set
# as opposed to just a subset, as was done in this workflow)
data(dex_tables)

# Step 15: re-load stageR and run it with a target alpha of 0.05
library(stageR)
remove <- function(x) substr(x, 1, 15)
pConfirmation <- matrix(dxr$pvalue, ncol = 1)
dimnames(pConfirmation) <- list(remove(dxr$featureID), "transcript")
pScreen <- qval
names(pScreen) <- remove(names(pScreen))
tx2gene <- as.data.frame(dxr[,c("featureID", "groupID")])
for (i in 1:2) tx2gene[,i] <- remove(tx2gene[,i])

# Generate a table with the OFDR control
# This table will show genes which have passed screening, as well as the transcripts
# associated with them
# Less than 5% of these genes either will not contain any DTU or will harbor falsely
# confirmed transcripts
stageR_obj <- stageRTx(pScreen = pScreen, pConfirmation = pConfirmation,
                       pScreenAdjusted = TRUE, tx2gene  = tx2gene)
stageR_obj <- stageWiseAdjustment(stageR_obj, method = "dtu", alpha = 0.05)
suppressWarnings({
  dex.padj <- getAdjustedPValues(stageR_obj, order   = FALSE,
                                 onlySignificantGenes = TRUE)
})
head(dex.padj)