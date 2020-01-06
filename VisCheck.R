#same analyses for vischeck data as with aud corr
#don't run in same window since the names are the same so it would overwrite

#set the working directory to where you can access the downloaded spreadsheets csv file
setwd("~/Desktop/MITAutismWork")
#read in the childbehavior spreadsheet
ChildBehavior <- read.csv("EEG and NEU Participant Tracking - Child Behavioral Questionnaires.csv", skip = 1, header = TRUE)
#add all the packages necessary
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
text.files <- list.files(path="~/Dropbox/Arushi_MIT/Stats_Regression (Arushi)/Slopes4Stats/VisCheck", recursive=T, pattern="*.mat", full.names=T) 
readFile <- function(f) { dat.fl <- readLines(f) } 
text.data <- sapply(text.files, readFile) 
cls <- names(text.data)
list.filenames <- list(cls)
#reads in the names of the folders (aka the participant IDs)
datadf <- as.data.frame(list.filenames)
colnames(datadf) <- "ID"
#remove any extra text so that only the numbered ID remains
datadf$ID <- gsub(".*ERP","",datadf$ID)
datadf$ID <- gsub("VisCheck.*", "", datadf$ID)
#read in each file from each folder and combine into one dataframe
slopesdata <- do.call(smartbind,lapply(cls,as.data.frame(readMat)))
slopesdata$value.warning <- NULL
#combine the participant IDs with the slopes data
slopesdata <- cbind(datadf, slopesdata)

#only include the slopes data
slopesdata <- slopesdata[ , -grep("peaks", colnames(slopesdata))]
#separate out min slopes from max slopes
maxslopes <- slopesdata
maxslopes <- maxslopes[ , -grep("min", colnames(maxslopes))]
minslopes <- slopesdata
minslopes <- minslopes[ , -grep("max", colnames(minslopes))]

#child behavior t tests in NT vs ASD

#remove columns 3 to 5
ChildBehavior <- ChildBehavior[ , -c(3:5)]
#remove any columns with all NA values
ChildBehavior <- ChildBehavior[,colSums(is.na(ChildBehavior))<nrow(ChildBehavior)]
#make columns 3 to 35 numeric values
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
pvals$pval <- (lapply(pvals$pval, function(x) as.numeric(as.character(x))))
pvals$bonferroni <- p.adjust(pvals$pval, method = "bonferroni")
pvals$ah <- ifelse(pvals$bonferroni < 0.05,1,0)

#match ID names
maxslopes$ID <- sub("_", "", maxslopes$ID)
maxslopes$ID <- sub("_", "", maxslopes$ID)
minslopes$ID <- sub("_", "", minslopes$ID)
minslopes$ID <- sub("_", "", minslopes$ID)

#merge data - max
ChildBehavior <- ChildBehavior[, -c(36:54)]
maxmerged <- merge(maxslopes, ChildBehavior, by.x = "ID", by.y = "ID" , all = FALSE)

#correlations - max
library(Hmisc)
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
res2 <- res2[!grepl("value.max", res2$column),]
res2 <- res2[grepl("value.max", res2$row), ]

#merge data - min
minmerged <- merge(minslopes, ChildBehavior, by.x = "ID", by.y = "ID" , all = FALSE)

#correlations - min
resmin <- rcorr(as.matrix(minmerged[,-c(1,34)]))
resmin2 <- flattenCorrMatrix(resmin$r, resmin$P)
resmin2 <- resmin2[!grepl("value.min", resmin2$column),]
resmin2 <- resmin2[grepl("value.min", resmin2$row), ]

#test that above is correct -- check
correlation.test <- cor.test(maxmerged$value.max.slope.1, maxmerged$Seeking.Seeker.Quadrant..Raw.Score.out.of.35.)

#plot each significant correlation
significantres2 <- subset(res2, ((res2$cor > 0.5) & (res2$p < 0.1)) | ((res2$cor < -0.5) & (res2$p < 0.1)), )
significantresmin2 <- subset(resmin2, ((resmin2$cor > 0.5) & (resmin2$p < 0.1)) | ((resmin2$cor < -0.5) & (resmin2$p < 0.1)), )

#changing labels of EEG channels in significant correlation dataframes
setwd("~/Desktop/MITAutismWork")
channelnames <- as.data.frame(readMat("EEGChannels.mat"))
channelnames <- channelnames[c(3,12), ]
channelnames <- as.data.frame(t(channelnames))
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

#maximums
maxmerged[,-c(1, 34)] <- as.data.frame(lapply(maxmerged[,-c(1,34)],function(x) as.numeric(as.character(x))))
significantres2[,1:2] <- lapply(significantres2[,1:2], as.character)
maxmerged[ ,"Diagnosis"] <- as.character(maxmerged[ ,"Diagnosis"])
maxmerged$plotcolor <- "black"
maxmerged$plotcolor[maxmerged$Diagnosis=="ASD"]="red"
setwd("~/Desktop/SinhaLab/Vis/max")
for (a in 1:nrow(significantres2)) {
  x <- significantres2[a, 1]
  y <- significantres2[a, 2]
  jpeg(paste(x, y, ".jpeg", sep = ""))
  plot(maxmerged[,x], maxmerged[,y], xlab = x, ylab = y, pch = 16, cex = 1.5, col = maxmerged$plotcolor)
  points <- WhereXY(maxmerged[,x], maxmerged[,y])
  title(paste("n=", points[5,5], "   ", "r=", round(significantres2[a,3], digits = 3), sep = ""))
  abline(lm(maxmerged[,y]~maxmerged[,x]), col = "red")
  dev.off()
}

#minimums
minmerged[ ,"Diagnosis"] <- as.character(minmerged[ ,"Diagnosis"])
minmerged$plotcolor <- "black"
minmerged$plotcolor[minmerged$Diagnosis=="ASD"]="red"
setwd("~/Desktop/SinhaLab/Vis/min")
significantresmin2[,1:2] <- lapply(significantresmin2[,1:2], as.character)
for (a in 1:nrow(significantresmin2)) {
  x <- significantresmin2[a, 1]
  y <- significantresmin2[a, 2]
  jpeg(paste(x, y, ".jpeg", sep = ""))
  plot(minmerged[,x], minmerged[,y], xlab = x, ylab = y, pch = 16, cex = 1.5, col = minmerged$plotcolor)
  points <- WhereXY(minmerged[,x], minmerged[,y])
  title(paste("n=", points[5,5], "   ", "r=", round(significantresmin2[a,3], digits = 3), sep = ""))
  abline(lm(minmerged[,y]~minmerged[,x]), col = "red")
  dev.off()
}

#anova
# One test could be to look at the whether the habituation (maximum slopes)
# differ greatly based on whether the ADOS score is 4/5, the ADOS score is
# 7/8/9, or the subject is NT. Or a test to look at whether the behavioral
# scores differ between these three groups. The purpose would be to add another
# level of analysis within the ASD subjects.

ChildBehaviorAnova <- ChildBehavior
ChildBehaviorAnova$ADOS.2.Comparison.Score <- as.numeric(as.character(ChildBehaviorAnova$ADOS.2.Comparison.Score))
ChildBehaviorAnova$Diagnosis <- as.character(ChildBehaviorAnova$Diagnosis)
#for the rows which have a diagnosis of ASD, change the diagnosis to ASDHigh based on what the score is
#!!!!!!! this will change with the different spreadsheets
for(i in c(1:5, 16:17, 21:22, 26)) {
  if(ChildBehaviorAnova[i, "ADOS.2.Comparison.Score"] == 7 | ChildBehaviorAnova[i, "ADOS.2.Comparison.Score"] == 8 |ChildBehaviorAnova[i, "ADOS.2.Comparison.Score"] == 9) {
    ChildBehaviorAnova[i, "Diagnosis"] <- "ASDHigh"
  }
}
#trying to make the above easier but haven't figured out a way that works
#asd <- "ASD"
#ChildBehaviorAnova[ChildBehaviorAnova$Diagnosis %in% asd,][6,1]
#ChildBehaviorAnova[is.element(ChildBehaviorAnova$Diagnosis, asd),][2]
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

anova1 <- anova1[c(3:29, 34:35),]
colnames(anova1)[[2]] <- "ASDHigh-ASD"
colnames(anova1)[[3]] <- "NT-ASD"
colnames(anova1)[[4]] <- "NT-ASDHigh"

###MAXIMUMS
#anova - habituation - max
maxdiagnose <- maxmerged[,c(1:34,65)]
maxdiagnose$ADOS.2.Comparison.Score <- as.numeric(as.character(maxdiagnose$ADOS.2.Comparison.Score))
maxdiagnose$Diagnosis <- as.character(maxdiagnose$Diagnosis)
#this changes based on what rows are ASD
for(i in c(1, 4, 8:11, 14, 18)) {
  if(maxdiagnose[i, "ADOS.2.Comparison.Score"] == 7 | maxdiagnose[i, "ADOS.2.Comparison.Score"] == 8 |maxdiagnose[i, "ADOS.2.Comparison.Score"] == 9) {
    maxdiagnose[i, "Diagnosis"] <- "ASDHigh"
  }
}
maxdiagnose$Diagnosis = factor(maxdiagnose$Diagnosis)
anovamax <- data.frame(0,0)
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

pvalslopes <- pvalslopes[2:nrow(pvalslopes),]
colnames(pvalslopes)[[1]] <- "channel"
colnames(pvalslopes)[[2]] <- "pval"
pvalslopes$pval <- (lapply(pvalslopes$pval, function(x) as.numeric(as.character(x))))
pvalslopes$ah <- ifelse(pvalslopes$pval < 0.05,1,0)

##MAXIMUMS
#adding scores and correlations
additivescorecorr <- data.frame(0,0)
colnames(additivescorecorr) <- c("correlation.p.value", "correlation.estimate")
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
sigaddscorecorr <- subset(additivescorecorr, additivescorecorr$correlation.p.value <= 0.05, )
sigaddscorecorr <- subset(sigaddscorecorr, sigaddscorecorr$correlation.estimate > 0.5 | sigaddscorecorr$correlation.estimate < -0.5, )

#table of r squared and p values
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
    tworegress = smartbind(tworegress, data.frame(x, y, z, test$adj.r.squared, coef(test)[2,4]))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
tworegress <- tworegress[2:379, 2:6]
colnames(tworegress)[[5]] <- "p.val"
tworegress <- tworegress[with(tworegress, order(-test.adj.r.squared, p.val)), ]

#plots of the highest r squareds
setwd("~/Desktop/SinhaLab/Vis/maxmultreg")
#these hundred plots are arbitrary, you can print more since there are a lot of significant r squared values
for (a in 1:100) {
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
#anova - habituation - min
mindiagnose <- minmerged[,c(1:34,65)]
mindiagnose$ADOS.2.Comparison.Score <- as.numeric(as.character(mindiagnose$ADOS.2.Comparison.Score))
mindiagnose$Diagnosis <- as.character(mindiagnose$Diagnosis)
for(i in c(1, 4, 8:11, 14, 18)) {
  if(mindiagnose[i, "ADOS.2.Comparison.Score"] == 7 | mindiagnose[i, "ADOS.2.Comparison.Score"] == 8 |mindiagnose[i, "ADOS.2.Comparison.Score"] == 9) {
    mindiagnose[i, "Diagnosis"] <- "ASDHigh"
  }
}
mindiagnose$Diagnosis = factor(mindiagnose$Diagnosis)
anovamin <- data.frame(0,0)
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

#slopes t tests- min
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

pvalslopesmin <- pvalslopesmin[2:nrow(pvalslopesmin),]
colnames(pvalslopesmin)[[1]] <- "channel"
colnames(pvalslopesmin)[[2]] <- "pval"
pvalslopesmin$pval <- (lapply(pvalslopesmin$pval, function(x) as.numeric(as.character(x))))
pvalslopesmin$ah <- ifelse(pvalslopesmin$pval < 0.05,1,0)

#adding scores and correlations
additive_score_min <- data.frame(0,0)
colnames(additive_score_min) <- c("correlation.p.value", "correlation.estimate")
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
    tworegress_min = smartbind(tworegress_min, data.frame(x, y, z, test$adj.r.squared, coef(test)[2,4]))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
tworegress_min <- tworegress_min[2:nrow(tworegress_min), 2:6]
colnames(tworegress_min)[[5]] <- "p.val"
tworegress_min <- tworegress_min[with(tworegress_min, order(-test.adj.r.squared, p.val)), ]

#plots of the highest r squareds
setwd("~/Desktop/SinhaLab/Vis/minmultreg")
#100 plots is arbitrary considering that there are a lot of significant r squareds
for (a in 1:100) {
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

#single correlations and factor analysis
#in the future, one option could be to make the dots different colors based on
#place on the ASD spectrum (to see whether the spread of ASD dots is due to
#different ADOS scores) 

#for the multiple regression, one option could be to
#restrict correlations to ADOS scores and subscores to find the correlations
#that aren't just with the ASD participants

#mode of the significant correlations to see which electrode channel shows up most often
#can do the same in the future with behaviors
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
getmode(sigaddscorecorr$colnames.maxmerged...i..)
#value.max.slope.T7
getmode(sigaddscore_min$colnames.minmerged...i..)
#value.min.slope.O2
getmode(significantres2$row)
#value.max.slope.O2
getmode(significantresmin2$row)
#value.min.slope.O2

#exporting
setwd("~/Desktop/MITAutismWork")
write.csv(significantres2, file = "vis_sig_max_corr.csv" )
write.csv(significantresmin2, file = "vis_sig_min_corr.csv")
write.csv(sigaddscorecorr, file = "vis_sig_max_addcorr.csv")
write.csv(sigaddscore_min, file = "vis_sig_min_addcorr.csv")
