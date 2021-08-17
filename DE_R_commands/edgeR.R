library("edgeR")

##2.1 Read in the data
targets <- read.delim(file = "groups.tsv", stringsAsFactors = FALSE)
dgList <- readDGE(targets, comment.char = "_")
dgList$samples$group <- relevel(dgList$samples$group, ref = "Non cancer")

#### Number of genes and samples:
dim(dgList)

##2.2 Filter the data
dgList <- dgList[rowSums(dgList$counts) >= 20, ]

#### Number of genes after filtering:
nrow(dgList)

dgList$samples$lib.size <- colSums(dgList$counts)

##2.3 Normalise the data
dgList <- calcNormFactors(dgList)

##2.4 Explore the data
plotMDS(dgList, col = as.integer(dgList$samples$group), labels = dgList$samples$patient)
dgList2 <- dgList[, which(dgList$samples$patient != "A13")]

##2.5 Perform a simple two group comparison

#2.5.1 Calculate the dispersion factors
dgList <- estimateCommonDisp(dgList)
dgList <- estimateTagwiseDisp(dgList)
plotBCV(dgList)

#2.5.2 Perform the differential expression calling
et <- exactTest(dgList)
num <- length(et$table[,1])
top <- topTags(et, n = num)
de_table <- top[top$table[,4] < 0.05,]

#### Number of differentially expressed genes:
nrow(de_table)

#2.5.3 Visualise the results
detags <- rownames(de_table)
plotSmear(et, de.tags = detags)
abline(h = c(-1, 1), col = "blue")

#2.5.4 Write the results to a file
write.table(de_table,"de.txt")

##2.6 Perform a matched pair comparison

#2.6.1 Specify the experiment design
patient <- as.factor(dgList$samples$patient)
group <- as.factor(dgList$samples$group)
design <- model.matrix(~ patient + group)

#2.6.2 Calculate the dispersion factors
dgGlm <- estimateGLMCommonDisp(dgList, design)
dgGlm <- estimateGLMTrendedDisp(dgGlm, design)
dgGlm <- estimateGLMTagwiseDisp(dgGlm, design)
plotBCV(dgGlm)

#2.6.3 Perform the differential expression calling
fit <- glmFit(dgGlm, design)
lrt <- glmLRT(fit)
num <- length(lrt$table[,1])
top <- topTags(lrt, n = num)
de_table_lrt <- top[top$table[,5] < 0.05,]

#### Number of differentially expressed genes:
nrow(de_table_lrt)

#2.6.4 Visualise the results
detags_glm <- rownames(de_table_lrt)
plotSmear(lrt, de.tags = detags_glm)
abline(h = c(-1, 1), col = "blue")

#2.6.5 Write the results to a file
write.table(de_table_lrt,"de_glm.txt")
