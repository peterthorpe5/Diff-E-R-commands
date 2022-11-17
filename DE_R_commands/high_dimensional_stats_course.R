
############## DAY 1 https://pad.carpentries.org/2022-09-27_ed-dash_high-dim-stats


# day 1 of the high dim stats course
library("lasso2")  
data(Prostate)  


view(data)

head(Prostate)

help(Prostate)

dim(Prostate)

names(Prostate)

pairs(Prostate)

# test 2 

cor(Prostate)


## correlation matrix for variables describing cancer/clinical variables
cor(Prostate[, c(1, 2, 4, 6, 9)])


## use linear regression to predict patient age from cancer progression variables
model <- lm(
  age ~ lcavol + lweight + lbph + lcp + lpsa + svi + gleason + pgg45,
  data = Prostate
)


summary(model)


## examine model residuals
plot(model)

# plot in in one wiondiw
par(mfrow=c(2,2))


plot(model)

par(mfrow=c(1,1))


#########################################
# DNA methylation


# methylation exmaple
# https://carpentries-incubator.github.io/high-dimensional-stats-r/01-introduction-to-high-dimensional-data/index.html

library("minfi")
library("here")
library("ComplexHeatmap")

methylation <- readRDS(here("data/methylation.rds"))
head(colData(methylation))

# assay extract the matrix
methyl_mat <- t(assay(methylation))
## calculate correlations between cells in matrix
cor_mat <- cor(methyl_mat)

hist(methyl_mat, breaks="FD", xlab="M-value")

head(cor_mat)

# make it a pretty output
knitr::kable(head(colData(methylation)), row.name=F)



age <- methylation$Age
lm_age_methyl1 <- lm(methyl_mat[1, ] ~ age)
lm_age_methyl1



plot(age, methyl_mat[1, ], xlab = "Age", ylab = "Methylation level", pch = 16)
abline(lm_age_methyl1)



library("broom")
tidy(lm_age_methyl1)


### we will do this with smoker

library("limma")

design_smoke <- model.matrix(~methylation$smoker)
methyl_mat <- assay(methylation)

fit_smoke <- lmFit(methyl_mat, design = design_smoke)
fit_smoke <- eBayes(fit_smoke)

toptab_smoke <- topTable(fit_smoke, coef = 2, number = nrow(fit_smoke))
plot(toptab_smoke$logFC, -log10(toptab_smoke$P.Value),
     xlab = "Effect size", ylab = bquote(-log[10](p)),
     pch = 19
)





p_raw <- toptab_age$P.Value
p_fwer <- p.adjust(p_raw, method = "bonferroni")
library("ggplot2")
ggplot() +
  aes(p_raw, p_fwer) +
  geom_point() +
  scale_x_log10() + scale_y_log10() +
  geom_abline(slope = 1, linetype = "dashed") +
  geom_hline(yintercept = 0.05, linetype = "dashed", col = "red") +
  geom_vline(xintercept = 0.05, linetype = "dashed", col = "red") +
  labs(x = "Raw p-value", y = "Bonferroni p-value")




fit_age <- lmFit(methyl_mat, design = design_age)
fit_age <- eBayes(fit_age)

toptab_age <- topTable(fit_age, coef = 2, number = nrow(fit_age))
orderEffSize <- rev(order(abs(toptab_age$logFC))) # order by effect size (absolute log-fold change)
head(toptab_age[orderEffSize, ])

p_raw <- toptab_age$P.Value
p_fwer <- p.adjust(p_raw, method = "hochberg")
library("ggplot2")
ggplot() +
  aes(p_raw, p_fwer) +
  geom_point() +
  scale_x_log10() + scale_y_log10() +
  geom_abline(slope = 1, linetype = "dashed") +
  geom_hline(yintercept = 0.05, linetype = "dashed", col = "red") +
  geom_vline(xintercept = 0.05, linetype = "dashed", col = "red") +
  labs(x = "Raw p-value", y = "hochberg p-value")


############################################################
# Day 3 

library("here")
library("minfi")
methylation <- readRDS(here("data/methylation.rds"))

## here, we transpose the matrix to have features as rows and samples as columns
# t to transpose
methyl_mat <- t(assay(methylation))
age <- methylation$Age

fit <- lm(age ~ methyl_mat)
# more feature that obser so, it doesnt work well .. 
# fewer obs that predictors
# Coefficients: (4964 not defined because of singularities)
summary(fit)



coef_horvath <- readRDS(here::here("data/coefHorvath.rds"))
methylation <- readRDS(here::here("data/methylation.rds"))

library("SummarizedExperiment")
age <- methylation$Age
methyl_mat <- t(assay(methylation))

# take the first 20 features, otherwise we cannot fit a model over ti. 
# now we have 20 bovs and 20 features
coef_horvath <- coef_horvath[1:20, ]
features <- coef_horvath$CpGmarker
# extract those which are in both the 
# known cpgmarkers. 
horvath_mat <- methyl_mat[, features]

## Generate an index to split the data
set.seed(42)
train_ind <- sample(nrow(methyl_mat), 25)


#Split the methylation data matrix and the age vector into training and test sets.

# we have already got a subset of the data, e.g. 20 row.
# from index vector (horvath_mat)

train_mat <- horvath_mat[train_ind, ]
train_age <- age[train_ind]
test_mat <- horvath_mat[-train_ind, ]
test_age <- age[-train_ind]



# Fit a model on the training data matrix and training age vector.
fit_horvath <- lm(train_age ~ train_mat)
fit_horvath <- lm(train_age ~ ., data = as.data.frame(train_mat))

mean(residuals(fit_horvath)^2)




mse <- function(true, prediction) {
  mean((true - prediction)^2)
}
pred_lm <- predict(fit_horvath, newdata = as.data.frame(test_mat))
err_lm <- mse(test_age, pred_lm)
err_lm



par(mfrow = c(1, 1))
plot(test_age, pred_lm, pch = 19)
abline(coef = 0:1, lty = "dashed")


# why would we want to restict our model?

library("glmnet")

## glmnet() performs scaling by default, supply un-scaled data:
horvath_mat <- methyl_mat[, features] # select the first 20 sites as before
train_mat <- horvath_mat[train_ind, ] # use the same individuals as selected before
test_mat <- horvath_mat[-train_ind, ]



ridge_fit <- glmnet(x = train_mat, y = train_age, alpha = 0)
plot(ridge_fit, xvar = "lambda")
abline(h = 0, lty = "dashed")


pred_ridge <- predict(ridge_fit, newx = test_mat)
err_ridge <- apply(pred_ridge, 2, function(col) mse(test_age, col))
min(err_ridge)



which_min_err <- which.min(err_ridge)
min_err_ridge <- min(err_ridge)
pred_min_ridge <- pred_ridge[, which_min_err]



chosen_lambda <- ridge_fit$lambda[which.min(err_ridge)]
plot(ridge_fit, xvar = "lambda")
abline(v = log(chosen_lambda), lty = "dashed")

# which performs beter. Ridge vs OLS: Ridged:

all <- c(pred_lm, test_age, pred_min_ridge)
lims <- range(all)
par(mfrow = 1:2)
plot(test_age, pred_lm,
     xlim = lims, ylim = lims,
     pch = 19
)
abline(coef = 0:1, lty = "dashed")
plot(test_age, pred_min_ridge,
     xlim = lims, ylim = lims,
     pch = 19
)
abline(coef = 0:1, lty = "dashed")


# laso model

fit_lasso <- glmnet(x = methyl_mat, y = age, alpha = 1)

par(mfrow = 1:2)
plot(fit_lasso, xvar = "lambda")
plot(fit_lasso)


# cv = cross validation, creates many tests sets, and thests them all. 
lasso <- cv.glmnet(methyl_mat[, -1], age, alpha = 1)
plot(lasso)



coefl <- coef(lasso, lasso$lambda.min)
selected_coefs <- as.matrix(coefl)[which(coefl != 0), 1]

## load the horvath signature to compare features
coef_horvath <- readRDS(here::here("data/coefHorvath.rds"))
## We select some of the same features! Hooray
intersect(names(selected_coefs), coef_horvath$CpGmarker)




# Fit an elastic net model (hint: alpha = 0.5) without cross-validation and plot the model object.
#  You can see that coefficients tend to go exactly to zero, but the paths are a bit less extreme than with pure LASSO; similar to ridge.

elastic <- glmnet(methyl_mat[, -1], age, alpha = 0.5)
par(mfrow = 1:1)

plot(elastic)


# Fit an elastic net model with cross-validation and plot the error. Compare with LASSO.
# The process of model selection is similar for elastic net models as for LASSO models.

elastic_cv <- cv.glmnet(methyl_mat[, -1], age, alpha = 0.5)
plot(elastic_cv)


# Select the lambda within one standard error of the minimum cross-validation error (hint: lambda.1se). Compare the coefficients with the LASSO model.
coefe <- coef(elastic_cv, elastic_cv$lambda.1se)
sum(coefe[, 1] == 0)

plot(
  coefl[, 1], coefe[, 1],
  pch = 16,
  xlab = "LASSO coefficients",
  ylab = "Elastic net coefficients"
)
abline(0:1, lty = "dashed")



############################################################################
# PCA



library("lasso2")
data("Prostate")

head(Prostate)



pros2 <- Prostate[, c("lcavol", "lweight", "lbph", "lcp", "lpsa")]
head(pros2)

# get the variance of the variable. 
apply(pros2, 2, var)

# get the mean ...
apply(pros2, 2, mean)


par(mfrow = c(1, 2))
# see if the variance is the same betwen them. They are not....
hist(pros2$lweight, breaks = "FD")
hist(pros2$lbph, breaks = "FD")


# need to scale the variable between them so they dont contribute more than th eothers

pca.pros <- prcomp(pros2, scale = TRUE, center = TRUE)
pca.pros


# how many pca comps do we needs
summary(pca.pros)


# calculate variance explained
varExp <- (pca.pros$sdev^2) / sum(pca.pros$sdev^2) * 100
# calculate percentage variance explained using output from the PCA
varDF <- data.frame(Dimensions = 1:length(varExp), varExp = varExp)
# create new dataframe with five rows, one for each principal component

pca.pros

par(mfrow = c(1, 1))
biplot(pca.pros, xlim = c(-0.3, 0.3))

# each num rps the patient in the dataset. 


### now use a PCA tool 
library("PCAtools")

library("SummarizedExperiment")
# here::here instead of loading the package ... 
cancer <- readRDS(here::here("data/cancer_expression.rds"))
mat <- assay(cancer)
metadata <- colData(cancer)

View(mat)


View(metadata)


all(colnames(mat) == rownames(metadata))


help(pca)



pca(
  mat,
  metadata = NULL,
  center = TRUE,
  scale = FALSE,
  rank = NULL,
  removeVar = NULL,
  transposed = FALSE,
  BSPARAM = ExactParam()
)

# need to save it as a variable

pc <- pca(
  mat,
  metadata = metadata,
  center = TRUE,
  scale = FALSE,
  rank = NULL,
  removeVar = NULL,
  transposed = FALSE)


biplot(pc)

# remove 20% variance
pc <- pca(mat, metadata = metadata, removeVar = 0.2)

# get the loading
pc$rotated[1:5, 1:5]

pc$loadings[1:5, 1:5]

# find the gene that affect PC1 the most

which.max(pc$loadings[, 2])
# 27
pc$loadings[27, ]

# cheating:
screeplot(pc, axisLabSize = 5, titleLabSize = 10)

# put the line at 80% - llook like PC24 explains 80% of the variantion
screeplot(pc, axisLabSize = 5, titleLabSize = 10, hline = 80)


# biplot

biplot(pc, lab = rownames(pc$metadata), pointSize = 1, labSize = 2)

# 
# grade is the cancer grade ...
biplot(pc, lab = NULL, colby = 'Grade', legendPosition = 'top')

# there are two groups

plotloadings(pc, labSize = 3)

# cluster explained by ER (oestrogen recetor + or -) ? 
biplot(pc, lab = NULL, shape = "ER", colby = 'Grade', legendPosition = 'top')

# do the clusters get explained by Age?
biplot(pc, lab = NULL, shape = "ER", colby = 'Age', legendPosition = 'top')


# proper solution

biplot(pc,
       lab = paste0(pc$metadata$Age,'years'),
       colby = 'ER',
       hline = 0, vline = 0,
       legendPosition = 'right')


pairsplot(pc)

########################################
# factor analysis






library(lasso2)
data(Prostate)

View(Prostate)

nrow(Prostate)

#select five log-transformed clinical variables for further analysis
pros2 <- Prostate[, c(1, 2, 4, 6, 9)]
head(pros2)


# Use the factanal() function to identify the minimum number of factors necessary to explain most of the variation in the data

#factanal(x, factors, data = NULL, covmat = NULL, n.obs = NA,
#         subset, na.action, start = NULL,
#         scores = c("none", "regression", "Bartlett"),
#         rotation = "varimax", control = NULL, ...)

factanal(pros2, 3)





# methylation exmaple
# https://carpentries-incubator.github.io/high-dimensional-stats-r/01-introduction-to-high-dimensional-data/index.html

library("minfi")
library("here")
library("ComplexHeatmap")

methylation <- readRDS(here("data/methylation.rds"))
head(colData(methylation))

# assay extract the matrix
methyl_mat <- t(assay(methylation))
## calculate correlations between cells in matrix
cor_mat <- cor(methyl_mat)

hist(methyl_mat, breaks="FD", xlab="M-value")

head(cor_mat)

# make it a pretty output
knitr::kable(head(colData(methylation)), row.name=F)



age <- methylation$Age
lm_age_methyl1 <- lm(methyl_mat[1, ] ~ age)
lm_age_methyl1



plot(age, methyl_mat[1, ], xlab = "Age", ylab = "Methylation level", pch = 16)
abline(lm_age_methyl1)



library("broom")
tidy(lm_age_methyl1)


### we will do this with smoker

library("limma")

design_smoke <- model.matrix(~methylation$smoker)
methyl_mat <- assay(methylation)

fit_smoke <- lmFit(methyl_mat, design = design_smoke)
fit_smoke <- eBayes(fit_smoke)

toptab_smoke <- topTable(fit_smoke, coef = 2, number = nrow(fit_smoke))
plot(toptab_smoke$logFC, -log10(toptab_smoke$P.Value),
     xlab = "Effect size", ylab = bquote(-log[10](p)),
     pch = 19
)





p_raw <- toptab_age$P.Value
p_fwer <- p.adjust(p_raw, method = "bonferroni")
library("ggplot2")
ggplot() +
  aes(p_raw, p_fwer) +
  geom_point() +
  scale_x_log10() + scale_y_log10() +
  geom_abline(slope = 1, linetype = "dashed") +
  geom_hline(yintercept = 0.05, linetype = "dashed", col = "red") +
  geom_vline(xintercept = 0.05, linetype = "dashed", col = "red") +
  labs(x = "Raw p-value", y = "Bonferroni p-value")




fit_age <- lmFit(methyl_mat, design = design_age)
fit_age <- eBayes(fit_age)

toptab_age <- topTable(fit_age, coef = 2, number = nrow(fit_age))
orderEffSize <- rev(order(abs(toptab_age$logFC))) # order by effect size (absolute log-fold change)
head(toptab_age[orderEffSize, ])

p_raw <- toptab_age$P.Value
p_fwer <- p.adjust(p_raw, method = "hochberg")
library("ggplot2")
ggplot() +
  aes(p_raw, p_fwer) +
  geom_point() +
  scale_x_log10() + scale_y_log10() +
  geom_abline(slope = 1, linetype = "dashed") +
  geom_hline(yintercept = 0.05, linetype = "dashed", col = "red") +
  geom_vline(xintercept = 0.05, linetype = "dashed", col = "red") +
  labs(x = "Raw p-value", y = "hochberg p-value")


#######
# kmers clustering



# example to be careful - some clustering does not actually work as you think

library("ggplot2")

# get a radom pop
simdata_2d <- data.frame(x = rnorm(n = 100, sd = 1), y = rnorm(n =100, sd = 2))

qplot(data = simdata_2d, x = x, y = y)


# fuck it up - assign them to a cluster
simdata_2d$cluster <- ifelse(
  simdata_2d$y < (simdata_2d$x * 0.5 + 0.9), "a", "b"
)


# see it. 
simdata_2d

qplot(data = simdata_2d, x = x, y = y, colour = cluster)

# add elipse to it ...
qplot(data = simdata_2d, x = x, y = y, colour = cluster) + stat_ellipse()


# get a radom pop
simdata_2d_part2 <- data.frame(x = rnorm(n = 100, sd = 1, mean = ), y = rnorm(n =100, sd = 2))

###########################################################################################

# example from pages
# scater single cell - kmers exmaple
library("SingleCellExperiment")
library("scater")


# do a PCA then kmers on the pca data
scrnaseq <- readRDS(here::here("data/scrnaseq.rds"))
scrnaseq <- runPCA(scrnaseq, ncomponents = 15)
# pull out the first two components
pcs <- reducedDim(scrnaseq)[, 1:2]

View(pcs)

# kmers starts with a random start point, so set this for reproducibility
set.seed(42)
cluster <- kmeans(pcs, centers = 4)
View(cluster)
help(kmeans)
# we have an obj called cluster made above, which have someting in it called cluster, which is 1 - 4 as th ekemans 4 clusters
# as we asked for. 
scrnaseq$kmeans <- as.character(cluster$cluster)
# factor would be better
scrnaseq$kmeans <- as.factor(cluster$cluster)


plotReducedDim(scrnaseq, "PCA", colour_by = "kmeans")

# another way to view, example
qplot(data = as.data.frame(pcs), x = PC1, y= PC2)

# alternate plot approach

pcs_clustered <- cbind(as.data.frame(pcs), data.frame(kmeans_4 = as.factor(cluster$cluster)))

# ggplot of the clustering stuff .... 
qplot(data = pcs_clustered, x = PC1, y = PC2, colour = kmeans_4) + stat_ellipse()


# kmean = 5

set.seed(42)
cluster_5 <- kmeans(pcs, centers = 5)
View(cluster_5)
help(kmeans)
# we have an obj called cluster made above, which have someting in it called cluster, which is 1 - 4 as th ekemans 4 clusters
# as we asked for. 
scrnaseq$kmeans_5 <- as.character(cluster_5$cluster)
# factor would be better
scrnaseq$kmeans_5 <- as.factor(cluster_5$cluster)

plotReducedDim(scrnaseq, "PCA", colour_by = "kmeans_5")

#########


x <- rnorm(20)
y <- rnorm(20)
x[10] <- x[10] + 10
plot(x, y, pch = 16)
points(mean(x), mean(y), pch = 16, col = "firebrick")
points(median(x), median(y), pch = 16, col = "dodgerblue")


library("cluster")
dist_mat <- dist(pcs)
sil <- silhouette(cluster$cluster, dist = dist_mat)
plot(sil, border = NA)



pc <- as.data.frame(pcs)
colnames(pc) <- c("x", "y")
pc$sil <- sil[, "sil_width"]
pc$clust <- factor(cluster$cluster)
mean(sil[, "sil_width"])


ggplot(pc) +
  aes(x, y, shape = clust, colour = sil) +
  geom_point() +
  scale_colour_gradient2(
    low = "dodgerblue", high = "firebrick"
  ) +
  scale_shape_manual(
    values = setNames(1:4, 1:4)
  )


# cluster robustness and th boostrap


data <- 1:5
sample(data, 5)

replicate(10, sample(data, 5))

# replace means you can get the same value more than once


replicate(10, sample(data, 5, replace = TRUE))

boots <- replicate(1000, mean(sample(data, 5, replace = TRUE)))
hist(boots,
     breaks = "FD",
     main = "1,000 bootstrap samples",
     xlab = "Mean of sample"
)


# apply it to clustering 
library("pheatmap")
library("bluster")
library("viridis") #  allows us to seclect colours

km_fun <- function(x) {
  kmeans(x, centers = 4)$cluster
}

# test the fucn works. 
kmeans_4_test <- km_fun(pcs)

# boot starp stb is a fucntion froone of the packages
ratios <- bootstrapStability(pcs, FUN = km_fun, clusters = cluster$cluster)

#result: telling us that point are in the same cluster a large proportion of th etime
#1         2         3         4
#1 0.9998021 1.0000000 1.0000000 1.0000653
#2        NA 0.9973625 0.9925105 1.0000653
#3        NA        NA 0.9857277 0.9970179
#4        NA        NA        NA 0.9998021


pheatmap(ratios,
         cluster_rows = FALSE, cluster_cols = FALSE,
         col = viridis(10),
         breaks = seq(0, 1, length.out = 10)
)


# repeat with k = 5

km_fun_k <- function(x, centres = 4) {
  kmeans(x, centres)$cluster
}


# make it generics!!!!
# boot starp stb is a fucntion froone of the packages
ratios <- bootstrapStability(pcs, FUN = km_fun_k, centres = 9
, clusters = cluster$cluster)



pheatmap(ratios,
         cluster_rows = FALSE, cluster_cols = FALSE,
         col = viridis(10),
         breaks = seq(0, 1, length.out = 10)
)





 #################################################################
 # higherchiachal clustering
 

library("minfi")
library("here")
library("ComplexHeatmap") # the dude uses pheatmap instead

methyl <- readRDS(here("data/methylation.rds"))

methyl_mat <- t(assay(methyl))


Heatmap(methyl_mat,
        name = "Methylation level",
        cluster_rows = TRUE, cluster_columns = TRUE,
        row_dend_width = unit(0.2, "npc"),
        column_dend_height = unit(0.2, "npc"),
        show_row_names = FALSE, show_column_names = FALSE
)



# creating a dendrogram

#First, create some example data with two variables x1 and x2
set.seed(450)
example_data <- data.frame(
  x1 = rnorm(20, 8, 4.5),
  x2 = rnorm(20, 6, 3.4)
)

#plot the data and name data points by row numbers
plot(example_data$x1, example_data$x2, type = "n")
text(
  example_data$x1,
  example_data$x2,
  labels = rownames(example_data),
  cex = 0.7
)


dist_m <- dist(example_data, method = "euclidean")

# how to vis a dist matirx
# https://stackoverflow.com/questions/3081066/what-techniques-exists-in-r-to-visualize-a-distance-matrix


dist_mi <- 1/dist_m # one over, as qgraph takes similarity matrices as input
library(qgraph)
#jpeg('example_forcedraw.jpg', width=1000, height=1000, unit='px')
qgraph(dist_mi, layout='spring', vsize=3)
#dev.off()


clust <- hclust(dist_m, method = "complete")
plot(clust)



## k is a user defined parameter determining
## the desired number of clusters at which to cut the treee
cutree(clust, k = 3)


## both give same results 

four_cut <- cutree(clust, h = 4)

## we can produce the cluster each observation belongs to
## using the mutate and count functions
library(dplyr)
example_cl <- mutate(example_data, cluster = four_cut)
count(example_cl, cluster)


#plot cluster each point belongs to on original scatterplot
library(ggplot2)
ggplot(example_cl, aes(x = x2, y = x1, color = factor(cluster))) + geom_point()


different linkage models


## create a distance matrix using euclidean method
distmat <- dist(methyl_mat)
## hierarchical clustering using complete method
clust <- hclust(distmat, method = "complete")
## plot resulting dendrogram
plot(clust)

## draw border around three clusters
rect.hclust(clust, k = 3, border = 2:6)
## draw border around two clusters
rect.hclust(clust, k = 2, border = 2:6)


###

## cut tree at height = 4
cut <- cutree(clust, h = 50)

library("dendextend")
avg_dend_obj <- as.dendrogram(clust)      
## colour branches of dendrogram depending on clusters
plot(color_branches(avg_dend_obj, h = 50))



clust <- hclust(distmat, method = "ward.D")
plot(clust)

# more here: 

https://carpentries-incubator.github.io/high-dimensional-stats-r/07-hierarchical/index.html



