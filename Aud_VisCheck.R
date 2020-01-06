setwd("~/Desktop/MITAutismWork")
#collect all the previously downloaded csv files from AudCorr and VisCheck
temp = list.files(pattern="*.csv")
#load those significant files into R environment
list2env(lapply(setNames(temp, make.names(gsub("*.csv$", "", temp))), read.csv), envir = .GlobalEnv)


#see which combinations of correlations are common across aud and visual
max_corr <- merge(aud_sig_max_corr,vis_sig_max_corr,by=c("row", "column")) 
min_corr <- merge(aud_sig_min_corr, vis_sig_min_corr, by=c("row", "column"))

max_corr$X.x <- NULL
max_corr$X.y <- NULL
min_corr$X.x <- NULL
min_corr$X.y <- NULL

max_addcorr <- merge(aud_sig_max_addcorr, vis_sig_max_addcorr, by=c("colnames.maxmerged...i..", "colnames.maxmerged...b..", "colnames.maxmerged...a...1.."))
min_addcorr <- merge(aud_sig_min_addcorr, vis_sig_min_addcorr, by=c("colnames.minmerged...i..", "colnames.minmerged...b..", "colnames.minmerged...a...1.."))

max_addcorr$X.x <- NULL
max_addcorr$X.y <- NULL
min_addcorr$X.x <- NULL
min_addcorr$X.y <- NULL

colnames(max_corr) <- gsub("x", "aud", colnames(max_corr))
colnames(max_corr) <- gsub("y", "vis", colnames(max_corr))

colnames(min_corr) <- gsub("x", "aud", colnames(min_corr))
colnames(min_corr) <- gsub("y", "vis", colnames(min_corr))

colnames(max_addcorr) <- gsub('[.]x', ".aud", colnames(max_addcorr))
colnames(max_addcorr) <- gsub('[.]y', ".vis", colnames(max_addcorr))

colnames(min_addcorr) <- gsub('[.]x', ".aud", colnames(min_addcorr))
colnames(min_addcorr) <- gsub('[.]y', ".vis", colnames(min_addcorr))

