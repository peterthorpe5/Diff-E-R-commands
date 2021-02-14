library(edgeR)

rnaseqMatrix = read.table("/storage/home/users/pjt6/newton/final_genome2/RNAseq_mapping/DE_analysis/Gp_genes.counts.matrix", header=T, row.names=1, com='')
rnaseqMatrix = round(rnaseqMatrix)
exp_study = DGEList(counts=rnaseqMatrix, group=factor(colnames(rnaseqMatrix)))
exp_study = calcNormFactors(exp_study)
exp_study$samples$eff.lib.size = exp_study$samples$lib.size * exp_study$samples$norm.factors
write.table(exp_study$samples, file="/storage/home/users/pjt6/newton/final_genome2/RNAseq_mapping/DE_analysis/Gp_genes.counts.matrix.TMM_info.txt", quote=F, sep="\t", row.names=F)
