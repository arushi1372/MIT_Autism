setwd("~/Desktop/MITAutismWork")
ChildBehavior <- read.csv("Copy of EEG and NEU Participant Tracking - Child Behavioral Questionnaires.csv", skip = 1, header = TRUE)
library(corrplot)
ChildBehavior <- ChildBehavior[ , -c(3:5)]
ChildBehavior <- ChildBehavior[,colSums(is.na(ChildBehavior))<nrow(ChildBehavior)]
ChildBehavior[,3:27] <- as.data.frame(lapply(ChildBehavior[,3:27],function(x) as.numeric(as.character(x))))
ChildBehavior <- ChildBehavior[c(1:27), -c(28:54)]
row.names(ChildBehavior) <- ChildBehavior[,1]
ChildBehavior <- ChildBehavior[,-c(1:2)]
ChildBehavior[is.na(ChildBehavior)] <- 0

ChildBehaviorsubset <- read.csv("Copy of EEG and NEU Participant Tracking - Without ADOS.csv", skip = 1, header = TRUE)
ChildBehaviorsubset <- ChildBehaviorsubset[ , -c(3:5)]
ChildBehaviorsubset <- ChildBehaviorsubset[,colSums(is.na(ChildBehaviorsubset))<nrow(ChildBehaviorsubset)]
ChildBehaviorsubset[,3:31] <- as.data.frame(lapply(ChildBehaviorsubset[,3:31],function(x) as.numeric(as.character(x))))
ChildBehaviorsubset <- ChildBehaviorsubset[c(1:27), ]
row.names(ChildBehaviorsubset) <- ChildBehaviorsubset[,1]
ChildBehaviorsubset <- ChildBehaviorsubset[,-c(1:2)]
ChildBehaviorsubset[is.na(ChildBehaviorsubset)] <- 0

#rework the spreadsheet to create a dataframe for the correlational plot
ChildBehaviorrenamed <- read.csv("Copy of EEG and NEU Participant Tracking - Renamed Child Behaviors.csv", skip = 1, header = TRUE)
ChildBehaviorrenamed <- ChildBehaviorrenamed[ , -c(3:5)]
ChildBehaviorrenamed <- ChildBehaviorrenamed[,colSums(is.na(ChildBehaviorrenamed))<nrow(ChildBehaviorrenamed)]
ChildBehaviorrenamed[,3:35] <- as.data.frame(lapply(ChildBehaviorrenamed[,3:35],function(x) as.numeric(as.character(x))))
ChildBehaviorrenamed <- ChildBehaviorrenamed[c(1:27), -c(36:ncol(ChildBehaviorrenamed))]
row.names(ChildBehaviorrenamed) <- ChildBehaviorrenamed[,1]
ChildBehaviorrenamed <- ChildBehaviorrenamed[,-c(1:2)]
ChildBehaviorrenamed[is.na(ChildBehaviorrenamed)] <- 0


#THIS IS REALLY HELPFUL - CREATE A CORRELATIONAL PLOT
#correlate the behavior dataframe
corr_mat=cor(ChildBehaviorrenamed,method="s")
#name the file that you want to download it to
pdf("foo2.pdf", width = 20, height = 20)
#create the correlation plot
corrplot(corr_mat, order = "hclust", tl.col='black', tl.cex=.75) 
dev.off()



#################random methods to try to use factor analysis (principal component analysis and factanal eigenvalues)
#standardize the scores
ChildBehavior_stan = as.data.frame(scale(ChildBehavior))
res = factanal(ChildBehavior_stan, factors = 10, rotation = "none")
res$loadings

#without rotation of the factors
loadings_factor1 = res$loadings[,1]
eigenv_factor1 = sum(loadings_factor1^2) 
eigenv_factor1
#proportion variance (eigenvalue/# of variables)
eigenv_factor1/23
#uniqueness
res$uniquenesses
loadings_seeker = res$loadings[1,]
communality_seeker = sum(loadings_seeker^2)
uniqueness_seeker = 1-communality_seeker
uniqueness_seeker

#with rotation of the factors
res1 = factanal(ChildBehavior_stan, factors = 10, rotation = "varimax", na.action = na.omit, lower = 0.01)
res1$loadings
load = res1$loadings[,1:2]
plot(load, type="n") # set up plot 
text(load, labels=names(ChildBehavior_stan), cex=.7)

fit.3 <- factanal(ChildBehavior_stan,factors=7,rotation="varimax", lower = 0.01)
print(fit.3)

#feature selection

#Principal Component Analysis
library(stats)
pca_existing <- prcomp(ChildBehavior, scale. = TRUE)
pca_existing$x
plot(pca_existing)
biplot(pca_existing, scale = 0)

scores_existing_df <- as.data.frame(pca_existing$x)
plot(PC1~PC2, data=scores_existing_df, cex = .1, lty = "solid")
text(PC1~PC2, data=scores_existing_df, labels=rownames(ChildBehavior),cex=.8)
set.seed(1234)
existing_clustering <- kmeans(ChildBehavior, centers = 3)
existing_cluster_groups <- existing_clustering$cluster
plot(PC1~PC2, data=scores_existing_df, cex = .1, lty = "solid", col=existing_cluster_groups)
text(PC1~PC2, data=scores_existing_df, labels=rownames(ChildBehavior),cex=.8, col=existing_cluster_groups)

plot(cumsum(pca_existing$sdev^2/sum(pca_existing$sdev^2)))
pc.use <- 6 # explains 90% of variance
trunc <- pca_existing$x[,1:pc.use] %*% t(pca_existing$rotation[,1:pc.use])

#and add the center (and re-scale) back to data
if(pca_existing$scale != FALSE){
  trunc <- scale(trunc, center = FALSE , scale=1/pca_existing$scale)
}
if(pca_existing$center != FALSE){
  trunc <- scale(trunc, center = -1 * pca_existing$center, scale=FALSE)
}
dim(trunc)
dim(ChildBehavior)
trunc <- as.data.frame(trunc)

# ChildBehaviort <- as.data.frame(t(ChildBehavior))
# ChildBehaviort$`116` <- NULL
# ChildBehaviort$`235 pending elligibility`<- NULL
# pca_existingt <- prcomp(ChildBehaviort, scale. = TRUE)
# pca_existingt$x
# plot(pca_existingt)
# scores_existing_dft <- as.data.frame(pca_existingt$x)
# plot(PC1~PC2, data=scores_existing_dft, cex = .1, lty = "solid")
# text(PC1~PC2, data=scores_existing_dft, labels=rownames(ChildBehaviort),cex=.8)
# set.seed(1234)
# existing_clusteringt <- kmeans(ChildBehaviort, centers = 3)
# existing_cluster_groupst <- existing_clusteringt$cluster
# plot(PC1~PC2, data=scores_existing_dft, cex = .1, lty = "solid", col=existing_cluster_groupst)
# text(PC1~PC2, data=scores_existing_dft, labels=rownames(ChildBehaviort),cex=.8, col=existing_cluster_groupst)

# plot(cumsum(res$sdev^2/sum(res$sdev^2)))
# pc.use <- 4 # explains 93% of variance
# trunc <- res$x[,1:pc.use] %*% t(res$rotation[,1:pc.use])
# 
# #and add the center (and re-scale) back to data
# if(res$scale != FALSE){
#   trunc <- scale(trunc, center = FALSE , scale=1/res$scale)
# }
# if(res$center != FALSE){
#   trunc <- scale(trunc, center = -1 * res$center, scale=FALSE)
# }
# dim(trunc); dim(Xt)
