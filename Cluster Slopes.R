#grouping the electrodes together and then correlating


#you can run this entire code for auditory data
#just change where the data is accessed from (for the slopes) and the working directories for where you want the plots to go

#set working directory
setwd("~/Desktop/MITAutismWork")
library("R.matlab", lib.loc="/Library/Frameworks/R.framework/Versions/3.3/Resources/library")
library(repmis)
n
library(BrailleR)
2
library(dplyr)
library(data.table)
library(gtools)

#import matlab files into r formatted data
#set path to where you have the slopes data
#CHANGE THIS TO AUD OR VISCHECK BASED ON STIMULI
text.files <- list.files(path="~/Dropbox/Arushi_MIT/Stats_Regression (Arushi)/Slopes4Stats/VisCheck", recursive=T, pattern="*.mat", full.names=T) 
readFile <- function(f) { dat.fl <- readLines(f) } 
text.data <- sapply(text.files, readFile) 
cls <- names(text.data)
list.filenames <- list(cls)
#create a list of the participant IDs (aka the folder names)
datadf <- as.data.frame(list.filenames)
colnames(datadf) <- "ID"
#clean up the IDs so that only the numbers remain
datadf$ID <- gsub(".*ERP","",datadf$ID)
#CHANGE THIS TO AUD OR VISCHECK BASED ON STIMULI
datadf$ID <- gsub("VisCheck.*", "", datadf$ID)
#read in the slopes data from each file in each folder
slopesdata <- do.call(smartbind,lapply(cls,as.data.frame(readMat)))
slopesdata$value.warning <- NULL
#combine the IDs with the slopes data
slopesdata <- cbind(datadf, slopesdata)

#only include the slopes data
slopesdata <- slopesdata[ , -grep("peaks", colnames(slopesdata))]
#separate out min slopes from max slopes
maxslopes <- slopesdata
maxslopes <- maxslopes[ , -grep("min", colnames(maxslopes))]
minslopes <- slopesdata
minslopes <- minslopes[ , -grep("max", colnames(minslopes))]

#match ID names
maxslopes$ID <- sub("_", "", maxslopes$ID)
maxslopes$ID <- sub("_", "", maxslopes$ID)
minslopes$ID <- sub("_", "", minslopes$ID)
minslopes$ID <- sub("_", "", minslopes$ID)

#######################
#set the working directory for where the eeg channel names are
setwd("~/Desktop/MITAutismWork")
channelnames <- as.data.frame(readMat("EEGChannels.mat"))
channelnames <- channelnames[c(3,12), ]
channelnames <- as.data.frame(t(channelnames))
#for each channel name in the dataframe, replace the number with the actual name
for (b in 1:nrow(channelnames)) {
  channel <- channelnames[b, 1]
  enumber <- channelnames[b, 2]
  current <- paste("value.max.slope.",enumber, sep = "")
  channels <- paste("value.max.slope.", channel, sep = "")
  current2 <- paste("value.min.slope.",enumber, sep = "")
  channel2 <- paste("value.min.slope.", channel, sep = "")
  colnames(maxslopes) <- gsub(paste("\\b", current, "\\b", sep = ""), channels, colnames(maxslopes))
  colnames(minslopes) <- gsub(paste("\\b", current2, "\\b", sep = ""), channel2, colnames(minslopes))
}

clustered_vis_max <- data.frame(0,0)
#for each slope group of electrodes, take the means and then include in the clustered_vis_max dataframe
#this is for grouping the maximum slopes
for(i in 1:nrow(maxslopes)) {
mean <- c(maxslopes[i, 'value.max.slope.Fp1'], maxslopes[i, 'value.max.slope.Fp2'], maxslopes[i,'value.max.slope.AF3'], maxslopes[i,'value.max.slope.AF4'])
x <- mean(mean, na.rm = TRUE)
clustered_vis_max[1, 1] <- "Fp1.Fp2.AF3.AF4"
clustered_vis_max[i+1, 1] <- x
}
for(i in 1:nrow(maxslopes)) {
  mean <- c(maxslopes[i, 'value.max.slope.O1'], maxslopes[i, 'value.max.slope.O2'], maxslopes[i,'value.max.slope.Oz'])
  x <- mean(mean, na.rm = TRUE)
  clustered_vis_max[1, 2] <- "O1.O2.Oz"
  clustered_vis_max[i+1, 2] <- x
}
for(i in 1:nrow(maxslopes)) {
  mean <- c(maxslopes[i, 'value.max.slope.O1'], maxslopes[i, 'value.max.slope.O2'], maxslopes[i,'value.max.slope.Oz'], maxslopes[i,'value.max.slope.PO3'], maxslopes[i,'value.max.slope.PO4'])
  x <- mean(mean, na.rm = TRUE)
  clustered_vis_max[1, 3] <- "O1.O2.Oz.PO3.PO4"
  clustered_vis_max[i+1, 3] <- x
}
for(i in 1:nrow(maxslopes)) {
  mean <- c(maxslopes[i, 'value.max.slope.C1'], maxslopes[i, 'value.max.slope.C2'], maxslopes[i,'value.max.slope.Cz'])
  x <- mean(mean, na.rm = TRUE)
  clustered_vis_max[1, 4] <- "C1.C2.Cz"
  clustered_vis_max[i+1, 4] <- x
}
for(i in 1:nrow(maxslopes)) {
  mean <- c(maxslopes[i, 'value.max.slope.P3'], maxslopes[i, 'value.max.slope.P4'], maxslopes[i,'value.max.slope.Pz'])
  x <- mean(mean, na.rm = TRUE)
  clustered_vis_max[1, 5] <- "P3.P4.Pz"
  clustered_vis_max[i+1, 5] <- x
}
for(i in 1:nrow(maxslopes)) {
  mean <- c(maxslopes[i, 'value.max.slope.T7'], maxslopes[i, 'value.max.slope.T8'])
  x <- mean(mean, na.rm = TRUE)
  clustered_vis_max[1, 6] <- "T7.T8"
  clustered_vis_max[i+1, 6] <- x
}
for(i in 1:nrow(maxslopes)) {
  mean <- c(maxslopes[i, 'value.max.slope.P3'], maxslopes[i, 'value.max.slope.P4'], maxslopes[i,'value.max.slope.P7'], maxslopes[i,'value.max.slope.P8'])
  x <- mean(mean, na.rm = TRUE)
  clustered_vis_max[1, 7] <- "P3.P4.P7.P8"
  clustered_vis_max[i+1, 7] <- x
}
for(i in 1:nrow(maxslopes)) {
  mean <- c(maxslopes[i, 'value.max.slope.FC5'], maxslopes[i, 'value.max.slope.FC1'], maxslopes[i,'value.max.slope.C3'])
  x <- mean(mean, na.rm = TRUE)
  clustered_vis_max[1, 8] <- "FC5.FC1.C3"
  clustered_vis_max[i+1, 8] <- x
}
for(i in 1:nrow(maxslopes)) {
  mean <- c(maxslopes[i, 'value.max.slope.C4'], maxslopes[i, 'value.max.slope.FC2'], maxslopes[i,'value.max.slope.FC6'])
  x <- mean(mean, na.rm = TRUE)
  clustered_vis_max[1, 9] <- "C4.FC2.FC6"
  clustered_vis_max[i+1, 9] <- x
}
for(i in 1:nrow(maxslopes)) {
  mean <- c(maxslopes[i, 'value.max.slope.Fz'], maxslopes[i, 'value.max.slope.Cz'], maxslopes[i,'value.max.slope.FC1'], maxslopes[i,'value.max.slope.FC2'])
  x <- mean(mean, na.rm = TRUE)
  clustered_vis_max[1, 10] <- "Fz.Cz.FC1.FC2"
  clustered_vis_max[i+1, 10] <- x
}
for(i in 1:nrow(maxslopes)) {
  mean <- c(maxslopes[i, 'value.max.slope.F3'], maxslopes[i, 'value.max.slope.Fz'], maxslopes[i,'value.max.slope.F4'], maxslopes[i,'value.max.slope.FC1'], maxslopes[i,'value.max.slope.FC2'])
  x <- mean(mean, na.rm = TRUE)
  clustered_vis_max[1, 11] <- "F3.Fz.F4.FC1.FC2"
  clustered_vis_max[i+1, 11] <- x
}

string <- as.character(clustered_vis_max[1,])
#clustered slopes for the minimum peak slopes
clustered_vis_min <- data.frame(0,0)
for(a in 1:length(string)) {
  for(i in 1:nrow(minslopes)) {
    c <- unlist(strsplit(string[1], ".", fixed=TRUE))
    mean <- c(minslopes[i, paste("value.min.slope.", c[1], sep = "")], minslopes[i, paste("value.min.slope.", c[2], sep = "")],  minslopes[i, paste("value.min.slope.", c[3], sep = "")], minslopes[i, paste("value.min.slope.", c[4], sep = "")], minslopes[i, paste("value.min.slope.", c[5], sep = "")]  )
    x <- mean(mean, na.rm = TRUE)
    clustered_vis_min[1, a] <- string[a]
    clustered_vis_min[i+1, a] <- x
  }
}

#merge data - max

#format your downloaded csv sheet as according to what it looks like now
ChildBehavior <- read.csv("Copy of EEG and NEU Participant Tracking - Correlated Totals.csv", skip = 1, header = TRUE)
ChildBehavior <- ChildBehavior[ , -c(3:5)]
ChildBehavior <- ChildBehavior[,colSums(is.na(ChildBehavior))<nrow(ChildBehavior)]
ChildBehavior[,3:ncol(ChildBehavior)] <- as.data.frame(lapply(ChildBehavior[,3:ncol(ChildBehavior)],function(x) as.numeric(as.character(x))))
ChildBehavior <- ChildBehavior[c(1:27),c(1:2, 9, 15, 39, 49, 52)]
#ChildBehavior <- ChildBehavior[, -c(25:54)]
#ChildBehavior <- ChildBehavior[c(1:27),-c(3:8)]
colnames(clustered_vis_max) <- clustered_vis_max[1,]
clustered_vis_max <- clustered_vis_max[-1, ]
clustered_vis_max <- cbind(datadf, clustered_vis_max)
clustered_vis_max$ID <- sub("_", "", clustered_vis_max$ID)
clustered_vis_max$ID <- sub("_", "", clustered_vis_max$ID)
vis_maxmerged <- merge(clustered_vis_max, ChildBehavior, by.x = "ID", by.y = "ID")

#correlations - max
library(Hmisc)
vis_maxmerged[sapply(vis_maxmerged, is.nan)] <- NA
res <- rcorr(as.matrix(vis_maxmerged[,-c(1,13)]))
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
resvis2 <- flattenCorrMatrix(res$r, res$P)
x <- c("Average", "Group")
#x <- c("SSP2", "score", "total", "X.4", "X.6", "Social", "Restricted", "Overall", "ADOS")
resvis2 <- resvis2[grepl(paste(x, collapse = "|"), resvis2$column),]
resvis2 <- resvis2[!grepl(paste(x, collapse = "|"), resvis2$row), ]
significantresvis2 <- subset(resvis2, ((resvis2$cor > 0.5) & (resvis2$p < 0.1)) | ((resvis2$cor < -0.5) & (resvis2$p < 0.1)), )

#maximums + vis correlational plots
vis_maxmerged[,-c(1, 13)] <- as.data.frame(lapply(vis_maxmerged[,-c(1,13)],function(x) as.numeric(as.character(x))))
resvis2[,1:2] <- lapply(resvis2[,1:2], as.character)
vis_maxmerged[ ,"Diagnosis"] <- as.character(vis_maxmerged[ ,"Diagnosis"])
vis_maxmerged$plotcolor <- "black"
vis_maxmerged$plotcolor[vis_maxmerged$Diagnosis=="ASD"]="red"
setwd("~/Desktop/SinhaLab/Aud/Factor/max")
for (a in 1:nrow(resvis2)) {
  tryCatch({
  x <- resvis2[a, 1]
  y <- resvis2[a, 2]
  jpeg(paste(x, y, ".jpeg", sep = ""))
  plot(vis_maxmerged[,x], vis_maxmerged[,y], xlab = x, ylab = y, pch = 16, cex = 1.5, col = vis_maxmerged$plotcolor)
  points <- WhereXY(vis_maxmerged[,x], vis_maxmerged[,y])
  title(paste("n=", points[5,5], "   ", "r=", round(resvis2[a,3], digits = 3), sep = ""))
  legend("topright", c("ASD","NT"), lty=c(1,1),lwd=c(2.5,2.5),col=c("red","black"))
  abline(lm(vis_maxmerged[,y]~vis_maxmerged[,x]), col = "red")
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  dev.off()
}

#minimums + vis correlation plots
colnames(clustered_vis_min) <- clustered_vis_min[1,]
clustered_vis_min <- clustered_vis_min[-1, ]
clustered_vis_min <- cbind(datadf, clustered_vis_min)
clustered_vis_min$ID <- sub("_", "", clustered_vis_min$ID)
clustered_vis_min$ID <- sub("_", "", clustered_vis_min$ID)
vis_minmerged <- merge(clustered_vis_min, ChildBehavior, by.x = "ID", by.y = "ID")

res <- rcorr(as.matrix(vis_minmerged[,-c(1,13)]))
resvismin2 <- flattenCorrMatrix(res$r, res$P)
#change the vector of subsetting just the behaviors in the column based on what you named in the spreadsheet
#this is for my own groups based on the correlation plot
x <- c("Average", "Group")
#this is for all total score correlations
#x <- c("SSP2", "score", "total", "X.4", "X.6", "Social", "Restricted", "Overall", "ADOS")
resvismin2 <- resvismin2[grepl(paste(x, collapse = "|"), resvismin2$column),]
resvismin2 <- resvismin2[!grepl(paste(x, collapse = "|"), resvismin2$row), ]
significantresvismin2 <- subset(resvismin2, ((resvismin2$cor > 0.5) & (resvismin2$p < 0.1)) | ((resvismin2$cor < -0.5) & (resvismin2$p < 0.1)), )

vis_minmerged[,-c(1, 13)] <- as.data.frame(lapply(vis_minmerged[,-c(1,13)],function(x) as.numeric(as.character(x))))
resvismin2[,1:2] <- lapply(resvismin2[,1:2], as.character)
vis_minmerged[ ,"Diagnosis"] <- as.character(vis_minmerged[ ,"Diagnosis"])
vis_minmerged$plotcolor <- "black"
vis_minmerged$plotcolor[vis_minmerged$Diagnosis=="ASD"]="red"
setwd("~/Desktop/SinhaLab/Aud/Factor/min")
for (a in 1:nrow(resvismin2)) {
  tryCatch({
    x <- resvismin2[a, 1]
    y <- resvismin2[a, 2]
    jpeg(paste(x, y, ".jpeg", sep = ""))
    plot(vis_minmerged[,x], vis_minmerged[,y], xlab = x, ylab = y, pch = 16, cex = 1.5, col = vis_minmerged$plotcolor)
    points <- WhereXY(vis_minmerged[,x], vis_minmerged[,y])
    title(paste("n=", points[5,5], "   ", "r=", round(resvismin2[a,3], digits = 3), sep = ""))
    legend("topright", c("ASD","NT"), lty=c(1,1),lwd=c(2.5,2.5),col=c("red","black"))
    abline(lm(vis_minmerged[,y]~vis_minmerged[,x]), col = "red")
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  dev.off()
}


#correlations with clusters of slopes and all behaviors
setwd("~/Desktop/MITAutismWork")
ChildBehaviorrenamed <- read.csv("Copy of EEG and NEU Participant Tracking - Renamed Child Behaviors.csv", skip = 1, header = TRUE)
ChildBehaviorrenamed <- ChildBehaviorrenamed[ , -c(3:5)]
ChildBehaviorrenamed <- ChildBehaviorrenamed[,colSums(is.na(ChildBehaviorrenamed))<nrow(ChildBehaviorrenamed)]
ChildBehaviorrenamed[,3:35] <- as.data.frame(lapply(ChildBehaviorrenamed[,3:35],function(x) as.numeric(as.character(x))))
ChildBehaviorrenamed <- ChildBehaviorrenamed[c(1:27), -c(36:ncol(ChildBehaviorrenamed))]
maxmerged2 <- merge(clustered_vis_max, ChildBehaviorrenamed, by.x = "ID", by.y = "ID")
res <- rcorr(as.matrix(maxmerged2[,-c(1,13)]))
resslopes <- flattenCorrMatrix(res$r, res$P)
x <- c("FP1", "O1", "C1", "P3", "T7", "FC5", "C4", "Fz", "F3")
resslopes <- resslopes[!grepl(paste(x, collapse = "|"), resslopes$column),]
resslopes <- resslopes[grepl(paste(x, collapse = "|"), resslopes$row), ]
significantresslopes <- subset(resslopes, ((resslopes$cor > 0.5) & (resslopes$p < 0.1)) | ((resslopes$cor < -0.5) & (resslopes$p < 0.1)), )

vis_maxmerged2[,-c(1, 13)] <- as.data.frame(lapply(vis_maxmerged2[,-c(1,13)],function(x) as.numeric(as.character(x))))
resslopes[,1:2] <- lapply(resslopes[,1:2], as.character)
vis_maxmerged2[ ,"Diagnosis"] <- as.character(vis_maxmerged2[ ,"Diagnosis"])
vis_maxmerged2$plotcolor <- "black"
vis_maxmerged2$plotcolor[vis_maxmerged2$Diagnosis=="ASD"]="red"
setwd("~/Desktop/SinhaLab/Aud/Factor/allbehaviors/max")
for (a in 1:nrow(significantresslopes)) {
  tryCatch({
    x <- significantresslopes[a, 1]
    y <- significantresslopes[a, 2]
    jpeg(paste(x, y, ".jpeg", sep = ""))
    plot(vis_maxmerged2[,x], vis_maxmerged2[,y], xlab = x, ylab = y, pch = 16, cex = 1.5, col = vis_maxmerged2$plotcolor)
    points <- WhereXY(vis_maxmerged2[,x], vis_maxmerged2[,y])
    title(paste("n=", points[5,5], "   ", "r=", round(significantresslopes[a,3], digits = 3), sep = ""))
    legend("topright", c("ASD","NT"), lty=c(1,1),lwd=c(2.5,2.5),col=c("red","black"))
    abline(lm(vis_maxmerged2[,y]~ vis_maxmerged2[,x]), col = "red")
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  dev.off()
}

#can do the above analysis with minimum slopes too
