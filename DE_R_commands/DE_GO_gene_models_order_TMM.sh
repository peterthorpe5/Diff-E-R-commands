#!/bin/bash
#$ -cwd

#Abort on any error,
set -e

# conda activate trinity
# need to put this into here:
# /shelf/apps/pjt6/conda/envs/trinity/lib/site_perl/5.26.2/DelimParser.pm
# https://salsa.debian.org/med-team/trinityrnaseq/blob/d5a2e951f12e8812d4dfb3582eb2c8bbf79298e1/PerlLib/DelimParser.pm
cuwodi=$1

outdir=DE_analysis_EdgeRLOG2

cd ${cuwodi}/
# remove previous runs as this will cuase failure
rm -rf DE_analysis_EdgeRLOG2 DE_analysis_DESeq2LOG2 DE_analysis_EdgeRLOG1.6

trin_path=/mnt/shared/scratch/pthorpe/apps/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/
v9=/mnt/shared/scratch/pthorpe/apps/trinityrnaseq-v2.9.0/Analysis/DifferentialExpression/


# run DE
${trin_path}/run_DE_analysis_GLM.pl --matrix \
*genes.counts.matrix --samples_file samples_described.txt --method edgeR --output ${outdir}


# dont need to run this command if using TPM from slamon or kalisto. 
${trin_path}/run_TMM_normalization_write_FPKM_matrix.pl --matrix *genes.counts.matrix --lengths genes_feature_lengths.txt

# QC steps

${trin_path}/PtR  --order_columns_by_samples_file  --matrix \
*genes.counts.matrix --samples ${cuwodi}/samples_described.txt --CPM --log2 --min_rowSums 10 --compare_replicates

wait
${trin_path}/PtR  --order_columns_by_samples_file  --matrix \
*genes.counts.matrix --min_rowSums 2 --samples ${cuwodi}/samples_described.txt --log2 --CPM --sample_cor_matrix

wait
${trin_path}/PtR --matrix *genes.counts.matrix --pdf_width 14 --pdf_height 14 \
-s ${cuwodi}/samples_described.txt --min_rowSums 10 --log2 --CPM --center_rows --prin_comp 3

# does not need gene lenghts for this method 
### ${trin_path}/run_TMM_normalization_write_FPKM_matrix.pl --matrix *genes.counts.matrix --just_do_TMM_scaling


## .genes.counts.matrix.TMM_normalized.FPKM_normalized.FPKM - if you use gene lengths
## *.TMM_normalized.FPKM - if you use -just_do_TMM_scaling

##################################################################################################################################
# Edge R log 2

cd ${cuwodi}/${outdir}
pvalues="0.05 1e-2 1e-3"
ptree="20 30 40 50 60 70"


for pvalue in ${pvalues}
do
    analyse_de="${trin_path}/analyze_diff_expr.pl  --order_columns_by_samples_file --matrix ../*.TMM_normalized.FPKM 
    -P ${pvalue} -C 2  --samples ${cuwodi}/samples_described.txt"
    # may need to add this in if too many DE genes: --max_DE_genes_per_comparison 50000
    # but his above option changes the name of the outputed R script, so wilol need to 
    # chnage inthe following commands
    echo ${analyse_de}
    eval ${analyse_de}
    for pt in ${ptree}
    do
        ptcmd="${trin_path}/define_clusters_by_cutting_tree.pl   --Ptree ${pt} -R diffExpr.P${pvalue}_C2.matrix.RData"
        echo ${ptcmd}
        eval ${ptcmd}
        mvcmd="mv ${cuwodi}/${outdir}/*_P_${pt}.heatmap.heatmap.pdf 
        ${cuwodi}/${outdir}/diffExpr.P${pvalue}_C2.matrix.RData.clusters_fixed_P_${pt}"
        echo ${mvcmd}
        eval ${mvcmd}
        # change the kmeans number here to force into K clusters
        ${trin_path}/define_clusters_by_cutting_tree.pl   --Ktree 8 -R diffExpr.P${pvalue}_C2.matrix.RData
    done
     analyse_de="${v9}/analyze_diff_expr.pl  --order_columns_by_samples_file --examine_GO_enrichment --include_GOplot --gene_lengths ${cuwodi}/gene_len_GO.txt 
    --GO_annots    ${cuwodi}/Go_final.txt --order_columns_by_samples_file  --matrix ../*.TMM_normalized.FPKM 
    -P ${pvalue} -C 2 --samples ${cuwodi}/samples_described.txt"
    echo ${analyse_de}
    eval ${analyse_de}
done

################################################################################################################
# DESeq2 log2

outdir=DE_analysis_DESeq2LOG2

cd ${cuwodi}/

deseq2="${trin_path}/run_DE_analysis.pl --matrix *genes.counts.matrix --samples_file ${cuwodi}/samples_described.txt 
        --method DESeq2 --output ${outdir}"
echo ${deseq2}
eval ${deseq2}

cd ${cuwodi}/${outdir}
pvalues="0.05 1e-2 1e-3 "
ptree="20 30 40 50 60 70"


for pvalue in ${pvalues}
do
    
    analyse_de="${trin_path}/analyze_diff_expr.pl  --order_columns_by_samples_file --matrix ../*.TMM_normalized.FPKM 
    -P ${pvalue} -C 2 --samples ${cuwodi}/samples_described.txt"
    #may need this if too many DE genes --max_DE_genes_per_comparison 50000
    echo ${analyse_de}
    eval ${analyse_de}
    for pt in ${ptree}
    do
        ptcmd="${trin_path}/define_clusters_by_cutting_tree.pl   --Ptree ${pt} -R diffExpr.P${pvalue}_C2.matrix.RData"
        echo ${ptcmd}
        eval ${ptcmd}
        mvcmd="mv ${cuwodi}/${outdir}/*_P_${pt}.heatmap.heatmap.pdf 
        ${cuwodi}/${outdir}/diffExpr.P${pvalue}_C2.matrix.RData.clusters_fixed_P_${pt}"
        echo ${mvcmd}
        eval ${mvcmd}
        # change the kmeans number here to force into K clusters
        ${trin_path}/define_clusters_by_cutting_tree.pl   --Ktree 8 -R diffExpr.P${pvalue}_C2.matrix.RData
    done
    analyse_de="${trin_path}/analyze_diff_expr.pl  --order_columns_by_samples_file --examine_GO_enrichment --include_GOplot --gene_lengths ${cuwodi}/gene_len_GO.txt 
    --GO_annots    ${cuwodi}/Go_final.txt --order_columns_by_samples_file  --matrix ../*.TMM_normalized.FPKM 
    -P ${pvalue} -C 2 --samples ${cuwodi}/samples_described.txt"
    # may need this if too many DE: --max_DE_genes_per_comparison 50000 
    echo ${analyse_de}
    eval ${analyse_de}
done


################################################################################################

outdir=DE_analysis_EdgeRLOG1.6

cd ${cuwodi}/

${trin_path}/run_DE_analysis.pl --matrix \
*genes.counts.matrix --samples_file samples_described.txt --method edgeR --output ${outdir}


##################################################################################################################################
# Edge R log 1.6

cd ${cuwodi}/${outdir}
pvalues="0.05 1e-2 1e-3"
ptree="20 30 40 50 60 70"

for pvalue in ${pvalues}
do

    analyse_de="${trin_path}/analyze_diff_expr.pl  --order_columns_by_samples_file  --matrix ../*.TMM_normalized.FPKM 
    -P ${pvalue} -C 1.6  --samples ${cuwodi}/samples_described.txt"
    # may need thi: --max_DE_genes_per_comparison 50000
    echo ${analyse_de}
    eval ${analyse_de}
    for pt in ${ptree}
    do
        ptcmd="${trin_path}/define_clusters_by_cutting_tree.pl   --Ptree ${pt} -R diffExpr.P${pvalue}_C1.6.matrix.RData"
        echo ${ptcmd}
        eval ${ptcmd}
        mvcmd="mv ${cuwodi}/${outdir}/*_P_${pt}.heatmap.heatmap.pdf 
        ${cuwodi}/${outdir}/diffExpr.P${pvalue}_C1.6.matrix.RData.clusters_fixed_P_${pt}"
        echo ${mvcmd}
        eval ${mvcmd}
        # change the kmeans number here to force into K clusters
        ${trin_path}/define_clusters_by_cutting_tree.pl   --Ktree 8 -R diffExpr.P${pvalue}_C1.6.matrix.RData
        #mv ${cuwodi}/${outdir}/*_P_${pt}.heatmap.heatmap.pdf 
        #${cuwodi}/${outdir}/diffExpr.P${pvalue}_C1.6.matrix.RData.clusters_fixed_P_${pt}
    done
    analyse_de="${trin_path}/analyze_diff_expr.pl  --order_columns_by_samples_file --examine_GO_enrichment --include_GOplot --gene_lengths ${cuwodi}/gene_len_GO.txt 
    --GO_annots    ${cuwodi}/Go_final.txt --order_columns_by_samples_file  --matrix ../*.TMM_normalized.FPKM 
    -P ${pvalue} -C 1.6 --max_DE_genes_per_comparison 50000 --samples ${cuwodi}/samples_described.txt"
    echo ${analyse_de}
    eval ${analyse_de}
done
