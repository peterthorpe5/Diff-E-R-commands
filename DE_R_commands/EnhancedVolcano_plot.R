#Volcano plot with EnhancedVolcano (https://github.com/kevinblighe/EnhancedVolcano)
#
#with results<- a dataframe containing the gene ids, log2FoldChange, adjusted-pvalues (padj)
#label_volcano<-a dataframe with the gene ids to label.
#Sec<-a dataframe with the gene ids of the genes encoding secreted proteins
#keyval.shape<-a numeric to assign a specific point shape to each category of interest (Sec/DEGs for instance)


#Here I make a volcano plot where I (1) label only the genes encoding effectors and (2) modify the size and point shape of the gene encoding secreted proteins. 
#shape 9 is a diamond, shape 16 a normal plain dot 
keyvals.shape_sec <- ifelse(
rownames(results) %in% unique(Sec$Gene), 9,
16)
names(keyvals.shape_sec)[keyvals.shape_sec == 9] <- 'Sec'
names(keyvals.shape_sec)[keyvals.shape_sec == 16] <- 'DEGs'

#I already have loaded in R the df Sec, results (can use the DEseq2 contrast alternatively) and label_volcano.

EnhancedVolcano(results, lab=results$Gene, selectLab = unique((label_volcano$Gene)), x = 'log2FoldChange', y = 'padj',xlim = c(-2.5, 2.5),xlab = bquote(~Log[2]~ 'fold change'),ylab = bquote(~-Log[10]~adjusted~italic(P)),FCcutoff = 0.5,pCutoff = 0.001,colAlpha = 0.6,pointSize = c(ifelse(results$Gene%in%Sec$Gene, 3, 1)),gridlines.major = FALSE,gridlines.minor = FALSE,vline=0.5,boxedLabels = FALSE,labCol = "black",labFace = "bold",title = "Regulated genes",shapeCustom = keyvals.shape_sec)