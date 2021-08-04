########################################### #
# interactive volcano plot
#############

# this works:
# import the data and packages
library(plotly)
library(ggplot2)
library(tidyverse)
cancer <- read.csv("data/cancer_gene_results.csv")

# custom thresholds
Pthreshold <- 0.01
logFC_threshold <- 2

# make a new column if the gene is sig yes or no
cancer$Significant<-ifelse(cancer$PValue <=Pthreshold, "Yes", "No")

# assign ggplot to var called volc
volc <- ggplot(data=cancer, aes(logFC, -log10(PValue), label=gene_id)) +
   geom_point(aes(col=Significant))  +
  geom_text_repel(aes(label=ifelse(abs(logFC >= logFC_threshold) 
                       & PValue < Pthreshold, as.character(gene_id), ""))) +
  theme_grey() +
  geom_vline(xintercept=c(logFC_threshold,-logFC_threshold),
             linetype = "dashed") +
  geom_hline(yintercept = -log10(Pthreshold), linetype = "dashed") +
  scale_color_manual(values = c("black","red"))
  
# can save voc as a standard image if required 
# this is the interactive pacakge
ggplotly(volc)

# save it to html
library(htmlwidgets)

g_plotly<-ggplotly(volc)
saveWidget(g_plotly, "volc_plotly.html")

######################################################
##################
##############

cancer <- read.csv("data/cancer_gene_results.csv")

Pthreshold <- 0.01
logFC_threshold <- 2

cancer$Significant<-ifelse(cancer$PValue <=Pthreshold, "Yes", "No")


#G < - cancer  %>% 

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)

volc <- ggplot(data=cancer, aes(logFC, -log10(PValue))) + 
  geom_point(aes(col=Significant))  +
  geom_text_repel(aes(label=ifelse(abs(logFC >= logFC_threshold) & 
                                     PValue < Pthreshold, 
                                   as.character(gene_id), ""), 
                      text=as.character(gene_id))) +
  theme_grey() +
  geom_vline(xintercept=c(logFC_threshold,-logFC_threshold),
             linetype = "dashed")+
  geom_hline(yintercept = -log10(Pthreshold), linetype = "dashed") +
  scale_color_manual(values = c("black","red")) 
  
  
ggplotly(volc, text = paste("geneid: ", cancer$gene_id,
               "<br>Pvalue: ", cancer$PValue,
               "<br>LogFC number: ", cancer$logFC), hoverinfo = 'text')


##################################################










genes<-read.csv("data/cancer_gene_results.csv")
genes$Significant<-ifelse(genes$PValue <=0.05, "Yes", "No")

ggplot(data=genes, aes(logFC, -log10(PValue))) + 
  geom_point(aes(col=Significant))+
  geom_text_repel(aes(label=ifelse(abs(logFC>=4) & PValue<0.05,as.character(gene_name),'')))+
  theme_classic()+
  geom_vline(xintercept=c(2,-2), linetype = "dashed")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")
  
  
or using tydyverse

cancer <-read.table("data/cancer_gene_results.csv", sep =",", header = TRUE)

cancer  %>% 


10 Tips to Customize Text Color, Font, Size in ggplot2 with element_text() - Python and R Tips

https://www.rdocumentation.org/packages/argparser/versions/0.7.1

cancer <- read.csv("data/cancer_gene_results.csv")

Pthreshold <- 0.01
logFC_threshold <- 2

cancer$Significant<-ifelse(cancer$PValue <=Pthreshold, "Yes", "No")

https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html

cancer  %>% 
  ggplot(aes(logFC, -log10(PValue))) + 
  geom_point(aes(col=Significant))  +
  geom_text_repel(aes(label=ifelse(abs(logFC >= logFC_threshold) & 
                                     PValue < Pthreshold, 
                                   as.character(gene_id), ''))) +
  theme_grey() +
  geom_vline(xintercept=c(logFC_threshold,-logFC_threshold),
             linetype = "dashed")+
  geom_hline(yintercept = -log10(Pthreshold), linetype = "dashed") +
  scale_color_manual(values = c("black","red"))
  
  
 #################
 interactive
 
 volc <- ggplot(data=cancer, aes(logFC, -log10(PValue))) + 
  geom_point(aes(col=Significant))  +
  geom_text_repel(aes(label=ifelse(abs(logFC >= logFC_threshold) & 
                                     PValue < Pthreshold, 
                                   as.character(gene_id), ""), 
                      text=as.character(gene_id))) +
  theme_grey() +
  geom_vline(xintercept=c(logFC_threshold,-logFC_threshold),
             linetype = "dashed")+
  geom_hline(yintercept = -log10(Pthreshold), linetype = "dashed") +
  scale_color_manual(values = c("black","red"))
  

ggplotly(volc)

 
library(htmlwidgets)

g_plotly<-ggplotly(volc)

saveWidget(g_plotly, "volc_plotly.html")


##################################################
# shiny   https://shiny.rstudio.com/tutorial/ 
 
library(shiny)














################################################################


#after the internet:

de <- cancer
  
  # The basic scatter plot: x is "log2FoldChange", y is "pvalue"
ggplot(data=de, aes(x=logFC, y=PValue)) + geom_point()
# Convert directly in the aes()
p <- ggplot(data=de, aes(x=logFC, y=-log10(PValue))) + geom_point()
# Add more simple "theme"
p <- ggplot(data=de, aes(x=logFC, y=-log10(PValue))) + geom_point() + theme_minimal()
# Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red")
    
# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)

# add a column of NAs
de$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
de$diffexpressed[de$logFC > 0.6 & de$PValue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
de$diffexpressed[de$logFC < -0.6 & de$PValue < 0.05] <- "DOWN"

# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(data=de, aes(x=logFC, y=-log10(PValue), col=diffexpressed)) + geom_point() + theme_minimal()

# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="black") +
        geom_hline(yintercept=-log10(0.05), col="black")
        
## Change point color 

# 1. by default, it is assigned to the categories in an alphabetical order):
p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))

# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)

# Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
de$delabel <- NA
de$delabel[de$diffexpressed != "NO"] <- de$gene_name[de$diffexpressed != "NO"]

ggplot(data=de, aes(x=de$logFC, y=-log10(de$PValue), col=diffexpressed, label=delabel)) + 
    geom_point() + 
    theme_minimal() +
    geom_text()


  