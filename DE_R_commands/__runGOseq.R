library(goseq)
library(GO.db)
library(qvalue)
# capture list of genes for functional enrichment testing
factor_labeling = read.table("Gp_genes.counts.matrix.Gp_J2_vs_Gp_MALE.DESeq2.DE_results.P1e-3_C2.DE.subset", row.names=1)
factor_labeling[,1] = rep('custom_list', dim(factor_labeling)[1])
factor_labeling = factor_labeling[,1,drop=F]
colnames(factor_labeling) = c('type')
factor_list = unique(factor_labeling[,1])
DE_genes = rownames(factor_labeling)


# get gene lengths
gene_lengths = read.table("/storage/home/users/pjt6/newton/final_genome2/RNAseq_mapping/DE_analysis/DE_with_GO_enrichment//gene_len_GO.txt", header=T, row.names=1, com='')
gene_lengths = as.matrix(gene_lengths[,1,drop=F])


# get background gene list
background = read.table("Gp_genes.counts.matrix.Gp_J2_vs_Gp_MALE.DESeq2.count_matrix", header=T, row.names=1)
background.gene_ids = rownames(background)
background.gene_ids = unique(c(background.gene_ids, DE_genes))


# parse GO assignments
GO_info = read.table("/storage/home/users/pjt6/newton/final_genome2/RNAseq_mapping/DE_analysis/DE_with_GO_enrichment//Go_final.txt", header=F, row.names=1,stringsAsFactors=F)
GO_info_listed = apply(GO_info, 1, function(x) unlist(strsplit(x,',')))
names(GO_info_listed) = rownames(GO_info)
get_GO_term_descr =  function(x) {
    d = 'none';
    go_info = GOTERM[[x]];
    if (length(go_info) >0) { d = paste(Ontology(go_info), Term(go_info), sep=' ');}
    return(d);
}


#organize go_id -> list of genes
GO_to_gene_list = list()
for (gene_id in names(GO_info_listed)) {
    go_list = GO_info_listed[[gene_id]]
    for (go_id in go_list) {
        GO_to_gene_list[[go_id]] = c(GO_to_gene_list[[go_id]], gene_id)
    }
}


# GO-Seq protocol: build pwf based on ALL DE features
sample_set_gene_ids = background.gene_ids
sample_set_gene_lengths = gene_lengths[sample_set_gene_ids,]
GO_info_listed = GO_info_listed[ names(GO_info_listed) %in% sample_set_gene_ids ]
cat_genes_vec = as.integer(sample_set_gene_ids %in% rownames(factor_labeling))
pwf=nullp(cat_genes_vec, bias.data=sample_set_gene_lengths)
rownames(pwf) = sample_set_gene_ids


# perform functional enrichment testing for each category.
for (feature_cat in factor_list) {
   message('Processing category: ', feature_cat)
    gene_ids_in_feature_cat = rownames(factor_labeling)[factor_labeling$type == feature_cat]
    cat_genes_vec = as.integer(sample_set_gene_ids %in% gene_ids_in_feature_cat)
    pwf$DEgenes = cat_genes_vec
    res = goseq(pwf,gene2cat=GO_info_listed, use_genes_without_cat=TRUE)
    ## over-represented categories:
     pvals = res$over_represented_pvalue
     pvals[pvals > 1 - 1e-10] = 1 - 1e-10
     q = qvalue(pvals)
     res$over_represented_FDR = q$qvalues
go_enrich_filename = paste("Gp_genes.counts.matrix.Gp_J2_vs_Gp_MALE.DESeq2.DE_results.P1e-3_C2.DE.subset", '.GOseq.enriched', sep='')
    result_table = res[res$over_represented_pvalue<=0.05,]
    descr = unlist(lapply(result_table$category, get_GO_term_descr))
    result_table$go_term = descr;
    result_table$gene_ids = do.call(rbind, lapply(result_table$category, function(x) { 
            gene_list = GO_to_gene_list[[x]]
            gene_list = gene_list[gene_list %in% rownames(factor_labeling)]
            paste(gene_list, collapse=', ');
     }) )
    write.table(result_table[order(result_table$over_represented_pvalue),], file=go_enrich_filename, sep='	', quote=F, row.names=F)
    ## under-represented categories:
     pvals = res$under_represented_pvalue
     pvals[pvals>1-1e-10] = 1 - 1e-10
     q = qvalue(pvals)
     res$under_represented_FDR = q$qvalues
    go_depleted_filename = paste("Gp_genes.counts.matrix.Gp_J2_vs_Gp_MALE.DESeq2.DE_results.P1e-3_C2.DE.subset", '.GOseq.depleted', sep='')
    result_table = res[res$under_represented_pvalue<=0.05,]
    descr = unlist(lapply(result_table$category, get_GO_term_descr))
    result_table$go_term = descr;
    write.table(result_table[order(result_table$under_represented_pvalue),], file=go_depleted_filename, sep='	', quote=F, row.names=F)
}
