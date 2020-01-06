
#calling EEG spreadsheet, dropbox slope files, channel names
#run it line by line
#in the future, change the color of the dots in the plots depending on whether ASD low or high to see whether the spread of ASD is due to a spectrum
#set the working directory for where the CSV downloaded spreadsheet is located
setwd("~/Desktop/MITAutismWork")

#read in the csv spreadsheet
ChildBehavior <- read.csv("Copy of EEG and NEU Participant Tracking - Child Behavioral Questionnaires.csv", skip = 1, header = TRUE)

#load packages necessary for future commands
library("R.matlab", lib.loc="/Library/Frameworks/R.framework/Versions/3.3/Resources/library")
library(repmis)
n
library(BrailleR)
2
library(dplyr)
library(data.table)
library(gtools)

#import matlab files into r formatted data

#set path to dropbox on your computer
text.files <- list.files(path="~/Dropbox/Arushi_MIT/Stats_Regression (Arushi)/Slopes4Stats/Aud", recursive=T, pattern="*.mat", full.names=T) 
readFile <- function(f) { dat.fl <- readLines(f) } 
text.data <- sapply(text.files, readFile) 
cls <- names(text.data)
list.filenames <- list(cls)
#reads in the names of the folders (aka the participant IDs)
datadf <- as.data.frame(list.filenames)
colnames(datadf) <- "ID"
#removes any extra text so that only the numbered ID remains
datadf$ID <- gsub(".*ERP","",datadf$ID)
datadf$ID <- gsub("Aud.*", "", datadf$ID)
#read in each file from each folder and merge it into one dataframe
slopesdata <- do.call(smartbind,lapply(cls,as.data.frame(readMat)))
slopesdata$value.warning <- NULL
#combine the participant IDs and the file data values
slopesdata <- cbind(datadf, slopesdata)

#only include the slopes data
slopesdata <- slopesdata[ , -grep("peaks", colnames(slopesdata))]
#separate out min slopes from max slopes
maxslopes <- slopesdata
maxslopes <- maxslopes[ , -grep("min", colnames(maxslopes))]
minslopes <- slopesdata
minslopes <- minslopes[ , -grep("max", colnames(minslopes))]

#child behavior t tests in NT vs ASD

#remove columns 3 to 5 (unneeded)
ChildBehavior <- ChildBehavior[ , -c(3:5)]
#remove any columns with all NA values
ChildBehavior <- ChildBehavior[,colSums(is.na(ChildBehavior))<nrow(ChildBehavior)]
#change the class of the data to numeric
ChildBehavior[,3:35] <- as.data.frame(lapply(ChildBehavior[,3:35],function(x) as.numeric(as.character(x))))

#run two sided ttests for each behavior
pvals <- as.data.frame(0,0)
for (i in 3:35) {
  x <- t.test(ChildBehavior[,i][ChildBehavior$Diagnosis == "ASD"], ChildBehavior[ ,i][ChildBehavior$Diagnosis == "NT"], na.rm = TRUE)
  print(x$p.value)
  pvals[i,] <- x$p.value
}

pvals <- as.data.frame(pvals[3:nrow(pvals),])
colnames(pvals) <- "pval"
#make the pvalues numeric
pvals$pval <- (lapply(pvals$pval, function(x) as.numeric(as.character(x))))
#use the bonferroni adjustment and add pvalues to the dataset
pvals$bonferroni <- p.adjust(pvals$pval, method = "bonferroni")
#alternative hypothesis true or false with significance value of 0.05
pvals$ah <- ifelse(pvals$bonferroni < 0.05,1,0)

### so pvals contains the pvalues for the ttests between NT and ASD behavioral scores

#match ID names across data
maxslopes$ID <- sub("_", "", maxslopes$ID)
maxslopes$ID <- sub("_", "", maxslopes$ID)
minslopes$ID <- sub("_", "", minslopes$ID)
minslopes$ID <- sub("_", "", minslopes$ID)

#merge data - max
#remove all behaviors after SCQ
ChildBehavior <- ChildBehavior[, -c(36:54)]
#merge the maximum peak slopes with the behavior scores
maxmerged <- merge(maxslopes, ChildBehavior, by.x = "ID", by.y = "ID" , all = FALSE)

#correlations - max
library(Hmisc)
#create one large matrix where each column is column with every other column in maxmerged
res <- rcorr(as.matrix(maxmerged[,-c(1,34)]))
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
res2 <- flattenCorrMatrix(res$r, res$P)
#remove any correlations that are slopes in the second variable
res2 <- res2[!grepl("value.max", res2$column),]
#keep only correlations that are slopes in the first variable
res2 <- res2[grepl("value.max", res2$row), ]
#result is a dataset (res2) of correlations between slopes (first variable) and behaviors (second variable)

#merge minimum slopes data with behavioral scores datad
minmerged <- merge(minslopes, ChildBehavior, by.x = "ID", by.y = "ID" , all = FALSE)

#correlations for the minimum slopes
#create the correlation matrix with each variable correlated to each variable
resmin <- rcorr(as.matrix(minmerged[,-c(1,34)]))
resmin2 <- flattenCorrMatrix(resmin$r, resmin$P)
#remove any correlations that are slopes in the second variable
resmin2 <- resmin2[!grepl("value.min", resmin2$column),]
#keep only correlations that are slopes in the first variable
resmin2 <- resmin2[grepl("value.min", resmin2$row), ]
#result is a dataset (resmin2) of correlations between slopes (first variable) and behaviors (second variable)

#test that above is correct -- crosscheck
correlation.test <- cor.test(maxmerged$value.max.slope.1, maxmerged$Seeking.Seeker.Quadrant..Raw.Score.out.of.35.)

#plot each significant correlation

#subset res2 to only include the significant correlations greater than 0.5 and with pvalues less than 0.1
significantres2 <- subset(res2, ((res2$cor > 0.5) & (res2$p < 0.1)) | ((res2$cor < -0.5) & (res2$p < 0.1)), )
#subset resmin2 to only include the significant correlations greater than 0.5 and with pvalues less than 0.1
significantresmin2 <- subset(resmin2, ((resmin2$cor > 0.5) & (resmin2$p < 0.1)) | ((resmin2$cor < -0.5) & (resmin2$p < 0.1)), )

#changing labels of EEG channels in significant correlation dataframes
#set the working directory to where the channel names for each electrode data is located
setwd("~/Desktop/MITAutismWork")
#read in the channel names data
channelnames <- as.data.frame(readMat("EEGChannels.mat"))
#include only the rows with the labels and the channel number
channelnames <- channelnames[c(3,12), ]
#transpose the dataset
channelnames <- as.data.frame(t(channelnames))
#for loop to replace the number with the channel name in significantres2 and significantresmin2
for (b in 1:nrow(channelnames)) {
  channel <- channelnames[b, 1]
  enumber <- channelnames[b, 2]
  current <- paste("value.max.slope.",enumber, sep = "")
  channels <- paste("value.max.slope.", channel, sep = "")
  current2 <- paste("value.min.slope.",enumber, sep = "")
  channel2 <- paste("value.min.slope.", channel, sep = "")
  significantres2$row <- gsub(paste("\\b", current, "\\b", sep = ""), channels, significantres2$row)
  significantresmin2$row <- gsub(paste("\\b", current2, "\\b", sep = ""), channel2, significantresmin2$row)
}

#changing labels of EEG channels in maxmerged and minmerged
#for loop to replace the number with the channel name in maxmerged and minmerged
for (b in 1:nrow(channelnames)) {
  channel <- channelnames[b, 1]
  enumber <- channelnames[b, 2]
  current <- paste("value.max.slope.",enumber, sep = "")
  channels <- paste("value.max.slope.", channel, sep = "")
  current2 <- paste("value.min.slope.",enumber, sep = "")
  channel2 <- paste("value.min.slope.", channel, sep = "")
  colnames(maxmerged)[[b+1]] <- gsub(paste("\\b", current, "\\b", sep = ""), channels, colnames(maxmerged)[[b+1]])
  colnames(minmerged)[[b+1]] <- gsub(paste("\\b", current2, "\\b", sep = ""), channel2, colnames(minmerged)[[b+1]])
}

#maximum peak slopes plots
###################################
#turn all the slopes values into numeric values
maxmerged[,-c(1, 34)] <- as.data.frame(lapply(maxmerged[,-c(1,34)],function(x) as.numeric(as.character(x))))
#make columns 1 and 2 character values
significantres2[,1:2] <- lapply(significantres2[,1:2], as.character)
#create a new column to assign plot color based on diagnosis
maxmerged[ ,"Diagnosis"] <- as.character(maxmerged[ ,"Diagnosis"])
maxmerged$plotcolor <- "black"
maxmerged$plotcolor[maxmerged$Diagnosis=="ASD"]="red"
#set the working directory to where you want to store the plots
setwd("~/Desktop/SinhaLab/Aud/max")
#for loop to generate plots of all the significantres2 correlations with regression lines
for (a in 1:nrow(significantres2)) {
  x <- significantres2[a, 1]
  y <- significantres2[a, 2]
  jpeg(paste(x, y, ".jpeg", sep = ""))
  plot(maxmerged[,x], maxmerged[,y], xlab = x, ylab = y, pch = 16, cex = 1.5, col = maxmerged$plotcolor)
  points <- WhereXY(maxmerged[,x], maxmerged[,y])
  #add a title showing n value and r value
  title(paste("n=", points[5,5], "   ", "r=", round(significantres2[a,3], digits = 3), sep = ""))
  legend("topright", c("ASD","NT"), lty=c(1,1),lwd=c(2.5,2.5),col=c("red","black"))
  abline(lm(maxmerged[,y]~maxmerged[,x]), col = "red")
  dev.off()
}

#minimum peak slopes plots
#################################################
#set the working directory to where you want to store the plots
setwd("~/Desktop/SinhaLab/Aud/min")
#turn columns 1 and 2 into character values
significantresmin2[,1:2] <- lapply(significantresmin2[,1:2], as.character)
#create a new column to assign plot color based on diagnosis
minmerged[ ,"Diagnosis"] <- as.character(minmerged[ ,"Diagnosis"])
minmerged$plotcolor <- "black"
minmerged$plotcolor[minmerged$Diagnosis=="ASD"]="red"
#for loop to generate plots of all the significantres2 correlations with regression lines
for (a in 1:nrow(significantresmin2)) {
  x <- significantresmin2[a, 1]
  y <- significantresmin2[a, 2]
  jpeg(paste(x, y, ".jpeg", sep = ""))
  plot(minmerged[,x], minmerged[,y], xlab = x, ylab = y, pch = 16, cex = 1.5, col = minmerged$plotcolor)
  points <- WhereXY(minmerged[,x], minmerged[,y])
  #add a title showing n value and r value
  title(paste("n=", points[5,5], "   ", "r=", round(significantresmin2[a,3], digits = 3), sep = ""))
  legend("topright", c("ASD","NT"), lty=c(1,1),lwd=c(2.5,2.5),col=c("red","black"))
  abline(lm(minmerged[,y]~minmerged[,x]), col = "red")
  dev.off()
}


#######################################################################
#anova
# One test could be to look at the whether the habituation (maximum slopes)
# differ greatly based on whether the ADOS score is 4/5, the ADOS score is
# 7/8/9, or the subject is NT. Or a test to look at whether the behavioral
# scores differ between these three groups. The purpose would be to add another
# level of analysis within the ASD subjects.

#create a new ChildBehavior dataset to modify for anova analyses
ChildBehaviorAnova <- ChildBehavior
#make the ADOS score numeric
ChildBehaviorAnova$ADOS.2.Comparison.Score <- as.numeric(as.character(ChildBehaviorAnova$ADOS.2.Comparison.Score))
#make the diagnosis a character value
ChildBehaviorAnova$Diagnosis <- as.character(ChildBehaviorAnova$Diagnosis)
#change the diagnosis to ASDHigh or ASD based on whether the scores in 4/5 or 7/8/9 !!!!!!!!!!!!! this changes based on which rows have ASD as the diagnosis
#c(1:5, 16:17, 21:22, 26) states the rows that say ASD
for(i in c(1:5, 16:17, 21:22, 26)) {
  if(ChildBehaviorAnova[i, "ADOS.2.Comparison.Score"] == 7 | ChildBehaviorAnova[i, "ADOS.2.Comparison.Score"] == 8 |ChildBehaviorAnova[i, "ADOS.2.Comparison.Score"] == 9) {
    ChildBehaviorAnova[i, "Diagnosis"] <- "ASDHigh"
  }
}
#use only the first 26 rows !!!!!!!!!!!!! this is viable to change based on the number of participants
ChildBehaviorAnova <- ChildBehaviorAnova[1:26, ]
ChildBehaviorAnova$Diagnosis = factor(ChildBehaviorAnova$Diagnosis)
anova1 <- data.frame(0,0)

#run anova on each behavior based on the diagnosis and us tukeyhsd to find what
#the p values are for each comparison (ASDHigh-ASD, ASD-NT, and ASDHigh-NT)
for(i in 3:ncol(ChildBehaviorAnova)) {
  tryCatch({
  x <- ChildBehaviorAnova[[i]]
  fitanova <- aov(x ~ Diagnosis, data = ChildBehaviorAnova)
  tuk <- TukeyHSD(fitanova)
  anova1[i, 2] <- tuk[[1]][1,4]
  if(is.na(anova1[i,3] <- tuk[[1]][2,4])) {next}
  anova1[i,3] <- tuk[[1]][2,4]
  anova1[i,4] <- tuk[[1]][3,4]
  anova1[i, 1] <- colnames(ChildBehaviorAnova)[[i]]
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

#make the dataframe pretty
anova1 <- anova1[c(3:29, 34:35),]
colnames(anova1)[[2]] <- "ASDHigh-ASD"
colnames(anova1)[[3]] <- "NT-ASD"
colnames(anova1)[[4]] <- "NT-ASDHigh"

###MAXIMUM Slopes
#anova to see whether there is a difference in the slopes across the three categories (ASDHigh, ASD, and NT)
maxdiagnose <- maxmerged[,c(1:34,57)]
maxdiagnose$ADOS.2.Comparison.Score <- as.numeric(as.character(maxdiagnose$ADOS.2.Comparison.Score))
maxdiagnose$Diagnosis <- as.character(maxdiagnose$Diagnosis)
#classify diagnosis based on ADOS score !!!!!!!!!! this changes based on which rows have ASD as the diagnosis
for(i in c(1:2, 5, 9:12, 16, 21)) {
  if(maxdiagnose[i, "ADOS.2.Comparison.Score"] == 7 | maxdiagnose[i, "ADOS.2.Comparison.Score"] == 8 |maxdiagnose[i, "ADOS.2.Comparison.Score"] == 9) {
    maxdiagnose[i, "Diagnosis"] <- "ASDHigh"
  }
}
maxdiagnose$Diagnosis = factor(maxdiagnose$Diagnosis)
anovamax <- data.frame(0,0)
#run anova on each behavior based on the diagnosis and us tukeyhsd to find what
#the p values are for each comparison (ASDHigh-ASD, ASD-NT, and ASDHigh-NT)
for(i in 2:33) {
  tryCatch({
    x <- maxdiagnose[[i]]
    anovamax[i, 1] <- colnames(maxdiagnose)[[i]]
    fitanova <- aov(x ~ Diagnosis, data = maxdiagnose)
    tuk <- TukeyHSD(fitanova)
    anovamax[i, 2] <- tuk[[1]]["ASDHigh-ASD",4]
    anovamax[i,3] <- tuk[[1]]["NT-ASD",4]
    anovamax[i,4] <- tuk[[1]]["NT-ASDHigh",4]
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
colnames(anovamax)[[2]] <- "ASDHigh-ASD"
colnames(anovamax)[[3]] <- "NT-ASD"
colnames(anovamax)[[4]] <- "NT-ASDHigh"
anovamax <- anovamax[2:33,]

#slopes t tests- maximum peak slopes
maxslopettest <- maxmerged[,c(1:34)]
pvalslopes <- data.frame(0,0)
for (i in 2:33) {
  tryCatch({
  pvalslopes[i,1] <- colnames(maxslopettest)[[i]]
  x <- t.test(maxslopettest[,i][maxslopettest$Diagnosis == "ASD"], maxslopettest[ ,i][maxslopettest$Diagnosis == "NT"], na.rm = TRUE)
  print(x$p.value)
  pvalslopes[i,2] <- x$p.value
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

#make the pvalslopes dataframe pretty/readable and understandable
pvalslopes <- pvalslopes[2:nrow(pvalslopes),]
colnames(pvalslopes)[[1]] <- "channel"
colnames(pvalslopes)[[2]] <- "pval"
pvalslopes$pval <- (lapply(pvalslopes$pval, function(x) as.numeric(as.character(x))))
pvalslopes$ah <- ifelse(pvalslopes$pval < 0.05,1,0)

##MAXIMUMS
#adding scores and correlations
library(gtools)
#create dataframe
additivescorecorr <- data.frame(0,0)
colnames(additivescorecorr) <- c("correlation.p.value", "correlation.estimate")
#for each of the slopes, for each of the behaviors, for every other behavior, create a correlation
for (i in 2:33) {
  for (b in 35:67) {
    for(a in b:66) {
      tryCatch({
      x <- maxmerged[,a+1]
      correlation <- cor.test(maxmerged[,i], maxmerged[,b] + x)
      #rowname <- paste(colnames(maxmerged)[[i]], ":", colnames(maxmerged)[[b]], "+", colnames(maxmerged)[[a+1]])
      additivescorecorr = smartbind(additivescorecorr, data.frame(colnames(maxmerged)[[i]], colnames(maxmerged)[[b]], colnames(maxmerged)[[a+1]], correlation$p.value, correlation$estimate))
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
  }
}

additivescorecorr <- additivescorecorr[2:nrow(additivescorecorr), ]
#subset the significant correlations into a separate dataset
sigaddscorecorr <- subset(additivescorecorr, additivescorecorr$correlation.p.value <= 0.05, )
sigaddscorecorr <- subset(sigaddscorecorr, sigaddscorecorr$correlation.estimate > 0.5 | sigaddscorecorr$correlation.estimate < -0.5, )

#table of r squared and p values from each significant additive correlation
tworegress <- data.frame(0)
for (a in 1:nrow(sigaddscorecorr)) {
  tryCatch({
  x <- as.character(sigaddscorecorr[a, 3])
  y <- as.character(sigaddscorecorr[a, 4])
  z <- as.character(sigaddscorecorr[a, 5])
  A <- maxmerged[,y]
  B <- maxmerged[,z]
  line <- lm((A + B)~maxmerged[,x])
  test <- summary(line)
  #from the regression line store the rsquared and pvalues
  tworegress = smartbind(tworegress, data.frame(x, y, z, test$adj.r.squared, coef(test)[2,4]))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
tworegress <- tworegress[2:379, 2:6]
colnames(tworegress)[[5]] <- "p.val"
#order the rsquared and pvalues in order of descending and ascending respectively
tworegress <- tworegress[with(tworegress, order(-test.adj.r.squared, p.val)), ]

#plots of the highest r squareds
#set working directory to where you want the plots to go
setwd("~/Desktop/SinhaLab/Aud/maxmultreg")
for (a in 1:100) { #arbitrary number of plots since many high r squared
  x <- as.character(tworegress[a, 1])
  y <- as.character(tworegress[a, 2])
  z <- as.character(tworegress[a, 3])
  A <- maxmerged[,y]
  B <- maxmerged[,z]
  jpeg(paste(x, y, z, ".jpeg", sep = ""))
  ylab <- paste(y, "+", z)
  plot(maxmerged[,x], A+B, xlab = x, ylab = ylab, pch = 16, cex = 1.5, col = maxmerged$plotcolor)
  points <- WhereXY(maxmerged[,x], A+B)
  title(paste("n=", points[5,5], "   ", "R^2=", round(tworegress[a,4], digits = 3), sep = ""))
  abline(lm(A+B~maxmerged[,x]), col = "red")
  dev.off()
}

#MINIMUMS
#anova of the minimum slopes based on ASDHigh, ASD, and NT
mindiagnose <- minmerged[,c(1:34,65)]
mindiagnose$ADOS.2.Comparison.Score <- as.numeric(as.character(mindiagnose$ADOS.2.Comparison.Score))
mindiagnose$Diagnosis <- as.character(mindiagnose$Diagnosis)
###this is subject to change depending on which rows in the dataframe are ASD
#c(1:2, 5, 9:12, 16, 21) are all rows diagnosed as ASD
for(i in c(1:2, 5, 9:12, 16, 21)) {
  if(mindiagnose[i, "ADOS.2.Comparison.Score"] == 7 | mindiagnose[i, "ADOS.2.Comparison.Score"] == 8 |mindiagnose[i, "ADOS.2.Comparison.Score"] == 9) {
    mindiagnose[i, "Diagnosis"] <- "ASDHigh"
  }
}
mindiagnose$Diagnosis = factor(mindiagnose$Diagnosis)
anovamin <- data.frame(0,0)
#run anova on each behavior based on the diagnosis and us tukeyhsd to find what
#the p values are for each comparison (ASDHigh-ASD, ASD-NT, and ASDHigh-NT)
for(i in 2:33) {
  tryCatch({
    x <- mindiagnose[[i]]
    anovamin[i, 1] <- colnames(mindiagnose)[[i]]
    fitanova <- aov(x ~ Diagnosis, data = mindiagnose)
    tuk <- TukeyHSD(fitanova)
    anovamin[i, 2] <- tuk[[1]]["ASDHigh-ASD",4]
    anovamin[i,3] <- tuk[[1]]["NT-ASD",4]
    anovamin[i,4] <- tuk[[1]]["NT-ASDHigh",4]
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
colnames(anovamin)[[2]] <- "ASDHigh-ASD"
colnames(anovamin)[[3]] <- "NT-ASD"
colnames(anovamin)[[4]] <- "NT-ASDHigh"
anovamin <- anovamin[2:33,]

#slopes t tests for minimums
minslopettest <- minmerged[,c(1:34)]
pvalslopesmin <- data.frame(0,0)
for (i in 2:33) {
  tryCatch({
    pvalslopesmin[i,1] <- colnames(minslopettest)[[i]]
    x <- t.test(minslopettest[,i][minslopettest$Diagnosis == "ASD"], minslopettest[ ,i][minslopettest$Diagnosis == "NT"], na.rm = TRUE)
    print(x$p.value)
    pvalslopesmin[i,2] <- x$p.value
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

#make the pvalslopesmin dataframe pretty/readable and understandable
pvalslopesmin <- pvalslopesmin[2:nrow(pvalslopesmin),]
colnames(pvalslopesmin)[[1]] <- "channel"
colnames(pvalslopesmin)[[2]] <- "pval"
pvalslopesmin$pval <- (lapply(pvalslopesmin$pval, function(x) as.numeric(as.character(x))))
pvalslopesmin$ah <- ifelse(pvalslopesmin$pval < 0.05,1,0)

#adding scores and correlations for minimum slopes
library(gtools)
#create new dataframe
additive_score_min <- data.frame(0,0)
colnames(additive_score_min) <- c("correlation.p.value", "correlation.estimate")
#for each slope, for each behavior, for every other behavior, correlation the two
for (i in 2:33) {
  for (b in 35:67) {
    for(a in b:66) {
      tryCatch({
        x <- minmerged[,a+1]
        correlation <- cor.test(minmerged[,i], minmerged[,b] + x)
        #rowname <- paste(colnames(maxmerged)[[i]], ":", colnames(maxmerged)[[b]], "+", colnames(maxmerged)[[a+1]])
        additive_score_min = smartbind(additive_score_min, data.frame(colnames(minmerged)[[i]], colnames(minmerged)[[b]], colnames(minmerged)[[a+1]], correlation$p.value, correlation$estimate))
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
  }
}

additive_score_min <- additive_score_min[2:nrow(additive_score_min), ]
#subset the significant correlations with good pvalues into separate dataframe
sigaddscore_min <- subset(additive_score_min, additive_score_min$correlation.p.value <= 0.05, )
sigaddscore_min <- subset(sigaddscore_min, sigaddscore_min$correlation.estimate > 0.5 | sigaddscore_min$correlation.estimate < -0.5, )

#table of r squared and p values
tworegress_min <- data.frame(0)
for (a in 1:nrow(sigaddscore_min)) {
  tryCatch({
    x <- as.character(sigaddscore_min[a, 3])
    y <- as.character(sigaddscore_min[a, 4])
    z <- as.character(sigaddscore_min[a, 5])
    A <- minmerged[,y]
    B <- minmerged[,z]
    line <- lm((A + B)~minmerged[,x])
    test <- summary(line)
    #from the regression line store the rsquared and pvalues
    tworegress_min = smartbind(tworegress_min, data.frame(x, y, z, test$adj.r.squared, coef(test)[2,4]))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
tworegress_min <- tworegress_min[2:nrow(tworegress_min), 2:6]
colnames(tworegress_min)[[5]] <- "p.val"
#order the rsquared and pvalues in order of descending and ascending respectively
tworegress_min <- tworegress_min[with(tworegress_min, order(-test.adj.r.squared, p.val)), ]

#plots of the highest r squareds
#set the working directory to where you want to store the plots generated
setwd("~/Desktop/SinhaLab/Aud/minmultreg")
for (a in 1:100) { #arbitrary number of plots for now
  x <- as.character(tworegress_min[a, 1])
  y <- as.character(tworegress_min[a, 2])
  z <- as.character(tworegress_min[a, 3])
  A <- minmerged[,y]
  B <- minmerged[,z]
  jpeg(paste(x, y, z, ".jpeg", sep = ""))
  ylab <- paste(y, "+", z)
  plot(minmerged[,x], A+B, xlab = x, ylab = ylab, pch = 16, cex = 1.5, col = minmerged$plotcolor)
  points <- WhereXY(minmerged[,x], A+B)
  title(paste("n=", points[5,5], "   ", "R^2=", round(tworegress_min[a,4], digits = 3), sep = ""))
  abline(lm(A+B~minmerged[,x]), col = "red")
  dev.off()
}

#modes of the most common electrodes
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#for the maximumslopes correlations with two behaviors added together, this was
#the most commonly well correlated electrode
getmode(sigaddscorecorr$colnames.maxmerged...i..)
#value.max.slope.F8

#for the minimumslopes correlations with two behaviors added together, this was
#the most commonly well correlated electrode
getmode(sigaddscore_min$colnames.minmerged...i..)
#value.min.slope.F3

#for the maximumslopes correlations with one behavior in auditory, this was the
#most commonly well correlated electrode
getmode(significantres2$row)
#value.max.slope.F8

#for the minimumslopes correlations with one behavior in auditory, this was the
#most commonly well correlated electrode
getmode(significantresmin2$row)
#value.min.slope.F3

#exporting this significant correlation tables
#set the working directory to where you want to export the files
setwd("~/Desktop/MITAutismWork")
write.csv(significantres2, file = "aud_sig_max_corr.csv" )
write.csv(significantresmin2, file = "aud_sig_min_corr.csv")
write.csv(sigaddscorecorr, file = "aud_sig_max_addcorr.csv")
write.csv(sigaddscore_min, file = "aud_sig_min_addcorr.csv")

#Compare WISC scores
...

#discount behaviors with not good t tests
# library(gtools)
# pvalssub <- subset(pvals, pvals$ah == 1, )
# pvalssub <- as.data.frame(t(pvalssub))
# ChildBehaviorSub <- merge(ChildBehavior, pvalssub, by.x = colnames(ChildBehavior), by.y = colnames(pvalssub), all.y = FALSE)
# ChildBehaviorSub <- smartbind(ChildBehavior, pvalssub)


#dropbox links
#https://www.dropbox.com/sh/b0l4bczaoit3h7o/AABSEoecSl9OROnekP_zuI2Ia?dl=0
#https://www.dropbox.com/sh/qx8b50qwalqkclj/AAAowwFjdFJOdH-roFf6hbDFa?dl=0

# 
# #nonlinear models
# # install.packages("mosaic")
# # library(mosaic)
# x <- significantres2[22, 1]
# y <- significantres2[22, 2]
# plot(maxmerged[,x], maxmerged[,y], xlab = x, ylab = y, pch = 16, cex = 1.5)
# abline(lm(maxmerged[,y]~maxmerged[,x]), col = "red")
# m2<-nls(Thought.Problems..unnormed.scaled.score ~ a/value.max.slope.26,data=maxmerged,start= list(a=1))
# v = maxmerged[(!(is.na(maxmerged$value.max.slope.26))) & (!(is.na(maxmerged$Thought.Problems..unnormed.scaled.score))), 15]
# abline(lm(maxmerged[,y]~log(maxmerged[,x])), col = "green")
# 
# new = data.frame(x = seq(min(maxmerged[,x]),max(maxmerged[,x]),len=200))
# lines(new$x,predict(m2,newdata=new), col = "orange")
