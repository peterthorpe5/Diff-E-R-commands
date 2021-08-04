#Introduction
    # This script is being written with the intention of performing PCA analysis to visualise underlying structure across RNAseq samples
#Install nerdy colour schemes
devtools::install_github("aljrico/gameofthrones")
# Load packages ------
library(tidyverse)
library(DT)
library(plotly)
library(gt)
library(gameofthrones)

# Identify variables of interest in study design file ----
#targets
#group <- targets$group
#group <- factor(group)


#Figure out what's different with Pete's data
Petes_data<-read.table('GM.isoform.counts.matrix.batch_eff_removal.matrix.minRow10.CPM.log2.dat', sep='\t')
Petes_cols<-colnames(Petes_data)

length(sampleIDs)
length(Petes_cols)
#filter Pete's data to match expectation of length
Petes_New_Targets<-filter(targets, targets$Sample %in% Petes_cols)
group <- Petes_New_Targets$group
group <- factor(group)


# Prepare your data -------
# use your normalized and filtered data in log2 cpm
#log2.GoodeidsCpmDGELIST.filtered.normalised

# Hierarchical clustering ---------------
#hierarchical clustering can only work on a data matrix, not a data frame
#try other distance methods (e.g. switch from 'maximum' to 'euclidean')...how does this change the result?
distance <- dist(t(Petes_data), method = "maximum") #other distance methods are "euclidean", maximum", "manhattan", "canberra", "binary" or "minkowski"
clusters <- hclust(distance, method = "average") #other agglomeration methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
plot(clusters, labels=Petes_cols, cex=0.5)

# Principal component analysis (PCA) -------------
pca.res <- prcomp(t(Petes_data), scale.=F, retx=T)
#look at the PCA result (pca.res) 
ls(pca.res)
summary(pca.res) # Prints variance summary for all principal components.
TranscriptScores_PCA<-pca.res$rotation #$rotation shows you how much each gene influenced each PC (called 'scores')
Loadings_Goodeids<-pca.res$x # 'x' shows you how much each sample influenced each PC (called 'loadings')
#note that these have a magnitude and a direction (this is the basis for making a PCA plot)
screeplot(pca.res) # A screeplot is a standard way to view eigenvalues for each PCA
pc.var<-pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per<-round(pc.var/sum(pc.var)*100, 1) # we can then use these eigenvalues to calculate the percentage variance explained by each PC
pc.per

pca.res.df <- pca.res$x[,1:8] %>% 
  as_tibble() %>%
  add_column(sample = Petes_cols,
             group = group)

# Categorize samples by tissue
pca.res.df <-  pca.res.df %>%
  mutate(Tissue=case_when(
    group == 'GMEmbryoXC' ~ "Embryo",
    group == 'GMEmbryoGA' ~ "Embryo",
    group == 'GMEmbryoGM' ~ "Embryo",
    group == 'GMFemale_BrainsGA' ~ "Female Brains",
    group == 'GMFemale_BrainsGM' ~ "Female Brains",
    group == 'GMFemale_BrainsXR' ~ "Female Brains",
    group == 'GMFemale_BrainsXC' ~ "Female Brains",
    group == 'GMMale_BrainsGA' ~ "Male Brains",
    group == 'GMMale_BrainXR' ~ "Male Brains",
    group == 'GMMale_BrainsGM' ~ "Male Brains",
    group == 'GMMale_BrainsXC' ~ "Male Brains",
    group == 'GMOvariesGA' ~ "Ovaries",
    group == 'GMOvariesGM' ~ "Ovaries",
    group == 'GMOvariesXR' ~ "Ovaries",
    group == 'GMOvariesXC' ~ "Ovaries",
    group == 'GMTestesGA' ~ "Testes",
    group == 'GMTestesGM' ~ "Testes",
    group == 'GMTestesXR' ~ "Testes",
    #group == 'GMTestesXC' ~ "Testes",
    group == 'GMTrophotaeniaGA' ~ "Trophotaenia",
    group == 'GMTrophotaeniaXC' ~ "Trophotaenia"))

pca.res.df <-  pca.res.df %>%
  mutate(Species=case_when(
    group == 'GMEmbryoXC' ~ "Xenoophorus captivus",
    group == 'GMEmbryoGA' ~ "Goodea atripinnis",
    group == 'GMEmbryoGM' ~ "Girardinichthys multiradiatus",
    group == 'GMFemale_BrainsGA' ~ "Goodea atripinnis",
    group == 'GMFemale_BrainsGM' ~ "Girardinichthys multiradiatus",
    group == 'GMFemale_BrainsXR' ~ "Xenotaenia resolanae",
    group == 'GMFemale_BrainsXC' ~ "Xenoophorus captivus",
    group == 'GMMale_BrainsGA' ~ "Goodea atripinnis",
    group == 'GMMale_BrainXR' ~ "Xenotaenia resolanae",
    group == 'GMMale_BrainsGM' ~ "Girardinichthys multiradiatus",
    group == 'GMMale_BrainsXC' ~ "Xenoophorus captivus",
    group == 'GMOvariesGA' ~ "Goodea atripinnis",
    group == 'GMOvariesGM' ~ "Girardinichthys multiradiatus",
    group == 'GMOvariesXR' ~ "Xenotaenia resolanae",
    group == 'GMOvariesXC' ~ "Xenoophorus captivus",
    group == 'GMTestesGA' ~ "Goodea atripinnis",
    group == 'GMTestesGM' ~ "Girardinichthys multiradiatus",
    group == 'GMTestesXR' ~ "Xenotaenia resolanae",
    #group == 'GMTestesXC' ~ "Xenoophorus captivus",
    group == 'GMTrophotaeniaGA' ~ "Goodea atripinnis",
    group == 'GMTrophotaeniaXC' ~ "Xenoophorus captivus"))

#Create traditional PCA plot seperated by PC axes that explain the most variation.
tiff("PCA_plot_GoodeidRNA_Species_Tissues_ForPete.tiff", units="in", width=9, height=6, res=300)
#pca.res.df <- as_tibble(pca.res$x)
ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=Petes_cols, colour=Tissue, shape=Species) +
  geom_point(size=4, alpha=0.7) +
  # geom_label() +
  #stat_ellipse() +
  scale_colour_got(option="Margaery", discrete = TRUE)+
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()
dev.off()


# Create a PCA 'small multiples' chart ----
# this is a way to view PCA laodings to understand impact of each sample on each pricipal component
pca.pivot <- pivot_longer(pca.res.df, # dataframe to be pivoted
                          cols = PC1:PC8, # column names to be stored as a SINGLE variable
                          names_to = "PC", # name of that new variable (column)
                          values_to = "loadings") # name of new variable (column) storing all the values (data)

ggplot(pca.pivot) +
  aes(x=Tissue, y=loadings, fill=Tissue) + 
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() +
  coord_flip()
