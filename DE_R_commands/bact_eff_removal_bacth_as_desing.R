library(edgeR)

rnaseqMatrix = read.table("XC.isoform.counts.matrix", header=T, row.names=1, com='', check.names=F)
batches = read.table("samples_described.txt", header=F, row.names=2, check.names=F)
batch_factors = as.factor(batches[colnames(rnaseqMatrix),])


exp_study = DGEList(counts=rnaseqMatrix)
exp_study = calcNormFactors(exp_study)
logCPM <- cpm(exp_study,log=TRUE)

seq_batches <- c('A', 'A', 'A', 'A', 'B', 'B', 'B', 'B', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'B', 'B', 'B', 'B', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'B', 'B', 'B', 'B', 'A', 'A', 'A', 'A', 'D', 'A', 'D', 'D', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'B', 'B', 'B', 'B', 'A', 'A', 'A', 'A', 'D', 'D', 'D', 'A', 'B', 'B', 'B', 'B', 'A', 'A', 'A', 'A', 'B', 'B', 'B', 'B', 'D', 'D', 'D', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A')

replicate_factors <- c('XCoutput_Ga_Br', 'XCoutput_Ga_Br', 'XCoutput_Ga_Br', 'XCoutput_Ga_Br', 'XCoutput_Ga_E', 'XCoutput_Ga_E', 'XCoutput_Ga_E', 'XCoutput_Ga_E', 'XCoutput_Ga_Fbr', 'XCoutput_Ga_Fbr', 'XCoutput_Ga_Fbr', 'XCoutput_Ga_Fbr', 'XCoutput_Ga_O', 'XCoutput_Ga_O', 'XCoutput_Ga_O', 'XCoutput_Ga_O', 'XCoutput_Ga_Te', 'XCoutput_Ga_Te', 'XCoutput_Ga_Te', 'XCoutput_Ga_Te', 'XCoutput_Ga_T', 'XCoutput_Ga_T', 'XCoutput_Ga_T', 'XCoutput_Ga_T', 'XCoutput_Gm_Br', 'XCoutput_Gm_Br', 'XCoutput_Gm_Br', 'XCoutput_Gm_Br', 'XCoutput_Gm_E', 'XCoutput_Gm_E', 'XCoutput_Gm_E', 'XCoutput_Gm_E', 'XCoutput_Gm_Fbr', 'XCoutput_Gm_Fbr', 'XCoutput_Gm_Fbr', 'XCoutput_Gm_Fbr', 'XCoutput_Gm_O', 'XCoutput_Gm_O', 'XCoutput_Gm_O', 'XCoutput_Gm_O', 'XCoutput_Gm_Te', 'XCoutput_Gm_Te', 'XCoutput_Gm_Te', 'XCoutput_Gm_Te', 'XCoutput_Xc_Br', 'XCoutput_Xc_Br', 'XCoutput_Xc_Br', 'XCoutput_Xc_Br', 'XCoutput_Xc_E', 'XCoutput_Xc_E', 'XCoutput_Xc_E', 'XCoutput_Xc_E', 'XCoutput_Xc_Fbr', 'XCoutput_Xc_Fbr', 'XCoutput_Xc_Fbr', 'XCoutput_Xc_Fbr', 'XCoutput_Xc_O', 'XCoutput_Xc_O', 'XCoutput_Xc_O', 'XCoutput_Xc_O', 'XCoutput_Xc_Te', 'XCoutput_Xc_Te', 'XCoutput_Xc_Te', 'XCoutput_Xc_Te', 'XCoutput_Xc_T', 'XCoutput_Xc_T', 'XCoutput_Xc_T', 'XCoutput_Xc_T', 'XCoutput_Xr_Br', 'XCoutput_Xr_Br', 'XCoutput_Xr_Br', 'XCoutput_Xr_Br', 'XCoutput_Xr_Fbr', 'XCoutput_Xr_Fbr', 'XCoutput_Xr_Fbr', 'XCoutput_Xr_Fbr', 'XCoutput_Xr_O', 'XCoutput_Xr_O', 'XCoutput_Xr_O', 'XCoutput_Xr_O', 'XCoutput_Xr_Te', 'XCoutput_Xr_Te', 'XCoutput_Xr_Te', 'XCoutput_Xr_Te')

Patient <- factor(c(8,8,33,33,51,51))
Tissue <- factor(c("N","T","N","T","N","T"))
data.frame(Sample=colnames(exp_study),seq_batches,replicate_factors)
design <- model.matrix(~seq_batches+replicate_factors)

logCPM <- removeBatchEffect(logCPM,batch=seq_batches, design=model.matrix(~replicate_factors))

millions = colSums(rnaseqMatrix)/1e6
myCPM = 2^logCPM
myRecounts = apply(myCPM, 1, function (x) x*millions)
write.table(t(myRecounts), file='XC.counts.matrix.batch_eff_removal.matrix', quote=F, sep='\t')
