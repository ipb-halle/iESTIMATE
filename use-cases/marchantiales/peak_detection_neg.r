


# ############################## MS1 ##############################



# ---------- MS1 Preparations ----------
# MS1 variables
polarity <- "negative"
pol <- substr(x=polarity, start=1, stop=3)

# General variables
mzml_files_neg <- NULL
mzml_names_neg <- NULL
mzml_times_neg <- NULL

# Load files
mzml_files_neg <- list.files(mzml_dir, pattern="*.mzML", recursive=T, full.names=T)
mzml_files_neg <- mzml_files_neg[grep(pol, mzml_files_neg, invert=FALSE)]

# Basenames of files without path and without extension
mzml_names_neg <- gsub('(.*)\\..*', '\\1', gsub('( |-|,)', '.', basename(mzml_files_neg)))

mzml_names_neg <- gsub(x=mzml_names_neg, pattern="^1\\.", replacement="R.cavernosa.SWE.", perl=TRUE)
mzml_names_neg <- gsub(x=mzml_names_neg, pattern="^11\\.", replacement="R.canaliculata.GOT.", perl=TRUE)
mzml_names_neg <- gsub(x=mzml_names_neg, pattern="^12\\.", replacement="R.bifurca.GOT.", perl=TRUE)
mzml_names_neg <- gsub(x=mzml_names_neg, pattern="^13\\.", replacement="R.ciliifera.GOT.", perl=TRUE)
mzml_names_neg <- gsub(x=mzml_names_neg, pattern="^17\\.", replacement="R.subbifurca.OE1.", perl=TRUE)
mzml_names_neg <- gsub(x=mzml_names_neg, pattern="^2\\.", replacement="R.sorocarpa.SWE.", perl=TRUE)
mzml_names_neg <- gsub(x=mzml_names_neg, pattern="^20\\.", replacement="R.beyrichiana.OEL.", perl=TRUE)
mzml_names_neg <- gsub(x=mzml_names_neg, pattern="^21\\.", replacement="R.subbifurca.OE2.", perl=TRUE)
mzml_names_neg <- gsub(x=mzml_names_neg, pattern="^25\\.", replacement="R.ciliifera.HAL.", perl=TRUE)
mzml_names_neg <- gsub(x=mzml_names_neg, pattern="^26\\.", replacement="R.gougetiana.HAL.", perl=TRUE)
mzml_names_neg <- gsub(x=mzml_names_neg, pattern="^45\\.", replacement="R.huebeneriana.GOT.", perl=TRUE)
mzml_names_neg <- gsub(x=mzml_names_neg, pattern="^5\\.", replacement="A.gracilis.SWE.", perl=TRUE)
mzml_names_neg <- gsub(x=mzml_names_neg, pattern="^7\\.", replacement="R.gothica.GOT.", perl=TRUE)
mzml_names_neg <- gsub(x=mzml_names_neg, pattern="^8\\.", replacement="M.fragrans.GOT.", perl=TRUE)
mzml_names_neg <- gsub(x=mzml_names_neg, pattern="^8c\\.", replacement="R.hemisphaerica.GOT.", perl=TRUE)
mzml_names_neg <- gsub(x=mzml_names_neg, pattern="^9\\.", replacement="A.hyalina.GOT.", perl=TRUE)

# Create phenodata based on species
mzml_pheno_neg <- data.frame(sample_name=mzml_names_neg, sample_group=as.factor(gsub(pattern="\\.\\d.*", replacement="", x=mzml_names_neg, perl=TRUE)))
mzml_pheno_neg$location <- as.factor(gsub(x=gsub(x=mzml_pheno_neg$sample_group, pattern=".*\\.", replacement=""), pattern="\\d", replacement="L"))
mzml_pheno_neg$species <- as.factor(gsub(x=mzml_pheno_neg$sample_group, pattern="\\.([^\\.]*)$", replacement=""))
mzml_pheno_samples_neg <- as.factor(mzml_pheno_neg$sample_group)
#mzml_pheno_colors_neg <- c("lemonchiffon4","orange4","mediumaquamarine","palegreen4","limegreen","royalblue","turquoise3","limegreen","palegreen4","seagreen4","tomato4","orchid3","green3","magenta3","violetred","mediumpurple3")
mzml_pheno_colors_neg <- c("orchid3","violetred","magenta3","turquoise3","mediumaquamarine","orange4","lemonchiffon4","palegreen4","steelblue3","seagreen4","mediumpurple3","tomato4","royalblue","limegreen")
#mzml_pheno_colors_neg <- RColorBrewer::brewer.pal(n=nlevels(mzml_pheno_samples_neg), name="Set1")
mzml_pheno_colors_samples_neg <- sapply(mzml_pheno_samples_neg, function(x) { x <- mzml_pheno_colors_neg[which(x==levels(mzml_pheno_samples_neg))] } )

# Save timestamps of samples
for (i in 1:length(mzml_files_neg)) {
	fl <- mzR::openMSfile(mzml_files_neg[i])
	run_info <- mzR::runInfo(fl)
	mzR::close(fl)
	mzml_times_neg <- c(mzml_times_neg, run_info$startTimeStamp)
}

# Display MSn levels
mzml_msn_neg <- NULL
for (i in 1:length(mzml_files_neg)) {
	mzml_data_neg <- readMSData(mzml_files_neg[i], mode="onDisk")
	mzml_msn_neg <- rbind(mzml_msn_neg, t(as.matrix(table(msLevel(mzml_data_neg)))))
}
colnames(mzml_msn_neg) <- c("MS1", "MS2")
rownames(mzml_msn_neg) <- mzml_names_neg

# Plot MSn levels
pdf(file="plots/neg_msn_levels.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=16, family="Helvetica")
par(mfrow=c(2,1), mar=c(16,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.8)
boxplot(mzml_msn_neg, main="Number of spectra")

model_boxplot <- boxplot(t(mzml_msn_neg[,2]), main="Number of MS2 spectra per sample", xaxt="n")
tick <- seq_along(model_boxplot$names)
axis(1, at=tick, labels=FALSE)
text(tick, par("usr")[3]-par("usr")[3]/10, model_boxplot$names, adj=0, srt=270, xpd=T)
dev.off()



# ---------- Peak detection ----------
# Import raw data as MSnbase object
raw_data_neg <- readMSData(files=mzml_files_neg, pdata=new("NAnnotatedDataFrame", mzml_pheno_neg), mode="onDisk", centroided=TRUE)
table(msLevel(raw_data_neg))
head(fData(raw_data_neg)[, c("isolationWindowTargetMZ", "isolationWindowLowerOffset",
							 "isolationWindowUpperOffset", "msLevel", "retentionTime")])
write.csv(fData(raw_data_neg), file="data/neg_raw_data.csv", row.names=FALSE)

# Restrict data to 1020 seconds (17 minutes)
raw_data_neg <- filterRt(raw_data_neg, c(min.rt, max.rt))

# Inspect mz values per file
raw_mz_neg <- mz(raw_data_neg)
raw_mz_neg <- split(raw_mz_neg, f=fromFile(raw_data_neg))
print(length(raw_mz_neg))

# Plot base peak chromatograms
chromas_bpc_neg <- chromatogram(raw_data_neg, aggregationFun="max")

pdf(file="plots/neg_chromas_bpc.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
plot(chromas_bpc_neg, main="Raw base peak chromatograms", xlab="retention time [s]", ylab="intensity", col=mzml_pheno_colors_samples_neg)
legend("topleft", bty="n", pt.cex=0.5, cex=0.7, y.intersp=0.7, text.width=0.5, pch=20, col=mzml_pheno_colors_neg, legend=levels(mzml_pheno_samples_neg))
dev.off()

# Plot TICs
chromas_tic_neg <- chromatogram(raw_data_neg, aggregationFun="sum")

pdf(file="plots/neg_chromas_tic.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
plot(chromas_tic_neg, main="Raw TIC chromatograms", xlab="retention time [s]", ylab="intensity", col=mzml_pheno_colors_samples_neg)
legend("topleft", bty="n", pt.cex=0.5, cex=0.7, y.intersp=0.7, text.width=0.5, pch=20, col=mzml_pheno_colors_neg, legend=levels(mzml_pheno_samples_neg))
dev.off()

# Get TICs
pdf(file="plots/neg_tics.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
par(mfrow=c(1,1), mar=c(19,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
tics_neg <- split(tic(raw_data_neg), f=fromFile(raw_data_neg))
boxplot(tics_neg, names=mzml_names_neg, las=2, col=mzml_pheno_colors_samples_neg, ylab="intensity", main="Total ion current")
dev.off()

# Grouping/binning the samples based on similarity of their base peak chromatogram to spot potentially problematic samples
chromas_bin_bpc_neg <- MSnbase::bin(chromas_bpc_neg, binSize=2)
chromas_bin_bpc_cor_neg <- cor(do.call(cbind, lapply(chromas_bin_bpc_neg, intensity)))
colnames(chromas_bin_bpc_cor_neg) <- rownames(chromas_bin_bpc_cor_neg) <- raw_data_neg$sample_name
pdf(file="plots/neg_chromas_bpc_bin_cor.pdf", encoding="ISOLatin1", pointsize=10, width=8, height=8, family="Helvetica")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
heatmap(chromas_bin_bpc_cor_neg, margins=c(15,15))
dev.off()

# Grouping/binning the samples based on similarity of their TIC to spot potentially problematic samples
chromas_bin_tic_neg <- MSnbase::bin(chromas_tic_neg, binSize=2)
chromas_bin_tic_cor_neg <- cor(do.call(cbind, lapply(chromas_bin_tic_neg, intensity)))
colnames(chromas_bin_tic_cor_neg) <- rownames(chromas_bin_tic_cor_neg) <- raw_data_neg$sample_name
pdf(file="plots/neg_chromas_tic_bin_cor.pdf", encoding="ISOLatin1", pointsize=10, width=8, height=8, family="Helvetica")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
heatmap(chromas_bin_tic_cor_neg, margins=c(15,15))
dev.off()

# Assess retention times and intensities of first file
head(rtime(chromas_neg[1, 1]))
head(intensity(chromas_neg[1, 1]))

# Peak detection in MS1 data
if (polarity=="negative") {
	ms1_params_neg <- CentWaveParam(ppm=25, mzCenterFun="mean", peakwidth=c(7.2, 36), prefilter=c(2, 50), mzdiff=0.0012, snthresh=7, noise=0, integrate=1,
									firstBaselineCheck=TRUE, verboseColumns=TRUE, fitgauss=FALSE, roiList=list(), roiScales=numeric())
} else {
	ms1_params_pos <- CentWaveParam(ppm=25, mzCenterFun="mean", peakwidth=c(9.4, 32), prefilter=c(6, 51), mzdiff=-0.0043, snthresh=2, noise=0, integrate=1,
									firstBaselineCheck=TRUE, verboseColumns=FALSE, fitgauss=FALSE, roiList=list(), roiScales=numeric())
}
ms1_data_neg <- findChromPeaks(raw_data_neg, param=ms1_params_neg)

# Per file summary
ms1_summary_neg <- lapply(split.data.frame(chromPeaks(ms1_data_neg), f=chromPeaks(ms1_data_neg)[, "sample"]), FUN=function(z) { c(peak_count=nrow(z), rt=quantile(z[, "rtmax"] - z[, "rtmin"])) } )
ms1_summary_neg <- do.call(rbind, ms1_summary_neg)
rownames(ms1_summary_neg) <- basename(fileNames(ms1_data_neg))
print(ms1_summary_neg)
table(msLevel(ms1_data_neg))
write.csv(as.data.frame(table(msLevel(ms1_data_neg))), file="data/neg_ms1_data.csv", row.names=FALSE)

# To get a global overview of the peak detection we can plot the frequency of identified peaks per file along the retention time axis. This allows to identify time periods along the MS run in which a higher number of peaks was identified and evaluate whether this is consistent across files.
pdf(file="plots/neg_ms1_data.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
par(mfrow=c(1,1), mar=c(4,18,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
plotChromPeakImage(ms1_data_neg, main="Frequency of identified peaks per RT")
dev.off()

# Group peaks
if (polarity=="negative") {
	ms1_data_neg <- groupChromPeaks(ms1_data_neg, param=PeakDensityParam(sampleGroups=ms1_data_neg$sample_group, minFraction=0.7, bw=0.25, minSamples=1, binSize=0.5))
} else {
	ms1_data_neg <- groupChromPeaks(ms1_data_neg, param=PeakDensityParam(sampleGroups=ms1_data_neg$sample_group, minFraction=0.7, bw=0.25, minSamples=1, binSize=0.5))
}

# RT correction
if (polarity=="negative") {
	ms1_data_neg <- adjustRtime(ms1_data_neg, param=PeakGroupsParam(minFraction=0.7, smooth="loess", span=0.5, extra=1, family="gaussian"))
} else {
	ms1_data_neg <- adjustRtime(ms1_data_neg, param=PeakGroupsParam(minFraction=0.7, smooth="loess", span=0.5, extra=1, family="gaussian"))
}

# Plot the difference of raw and adjusted retention times
pdf(file="plots/neg_ms1_raw_adjusted.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=8, family="Helvetica")
par(mfrow=c(2,1), mar=c(4.5,4.2,4,1), cex=0.8)
plot(chromas_neg, main="Raw chromatograms", col=mzml_pheno_colors_samples_neg)
legend("topleft", bty="n", pt.cex=0.5, cex=0.7, y.intersp=0.7, text.width=0.5, pch=20, col=mzml_pheno_colors_neg, legend=levels(mzml_pheno_samples_neg))
plotAdjustedRtime(ms1_data_neg, lwd=2, main="Retention Time correction", col=mzml_pheno_colors_samples_neg)
legend("topleft", bty="n", pt.cex=0.5, cex=0.7, y.intersp=0.7, text.width=0.5, pch=20, col=mzml_pheno_colors_neg, legend=levels(mzml_pheno_samples_neg))
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.8)
dev.off()

# Group peaks
if (polarity=="negative") {
	ms1_data_neg <- groupChromPeaks(ms1_data_neg, param=PeakDensityParam(sampleGroups=ms1_data_neg$sample_group, minFraction=0.7, bw=0.25, minSamples=1, binSize=0.5))
} else {
	ms1_data_neg <- groupChromPeaks(ms1_data_neg, param=PeakDensityParam(sampleGroups=ms1_data_neg$sample_group, minFraction=0.7, bw=0.25, minSamples=1, binSize=0.5))
}

# Get integrated peak intensity per feature/sample
print(head(featureValues(ms1_data_neg, value="into")))

# Fill peaks
#ms1_data_neg <- fillChromPeaks(ms1_data_neg, param=FillChromPeaksParam(ppm=ppm, fixedRt=0, expandRt=5))
#head(featureValues(ms1_data_neg))
#head(featureSummary(ms1_data_neg, group=ms1_data_neg$sample_group))

# Evaluate grouping
pdf(file="plots/neg_ms1_grouping.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
ms1_pca_neg <- prcomp(t(na.omit(log2(featureValues(ms1_data_neg, value="into")))), center=TRUE)
plot(ms1_pca_neg$x[, 1], ms1_pca_neg$x[, 2], pch=19, main="PCA: Grouping of samples",
	 xlab=paste0("PC1: ", format(summary(ms1_pca_neg)$importance[2, 1] * 100, digits=3), " % variance"),
	 ylab=paste0("PC2: ", format(summary(ms1_pca_neg)$importance[2, 2] * 100, digits=3), " % variance"),
	 col=mzml_pheno_colors_samples_neg, cex=2)
grid()
text(ms1_pca_neg$x[, 1], ms1_pca_neg$x[, 2], labels=ms1_data_neg$sample_name, col=mzml_pheno_colors_samples_neg, pos=3, cex=0.5)
dev.off()

# Evaluate QC Deviations
pdf(file="plots/neg_ms1_mzrt_deviations.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=8, family="Helvetica")
par(mfrow=c(2,1), mar=c(14.5,4,4,1), cex=0.8, cex.axis=0.6)
plotQC(ms1_data_neg, what="mzdevsample", sampColors=mzml_pheno_colors_samples_neg, sampNames=mzml_pheno_neg$sample_name)
title(main="Median mz deviation for each sample")
plotQC(ms1_data_neg, what="rtdevsample", sampColors=mzml_pheno_colors_samples_neg, sampNames=mzml_pheno_neg$sample_name)
title(main="Median RT deviation for each sample")
dev.off()

# Show peaks
tail(chromPeaks(ms1_data_neg))
tail(chromPeakData(ms1_data_neg))

# Show process history
processHistory(ms1_data_neg)



# ---------- Build MS1 feature tables ----------
# Build feature matrix
ms1_matrix_neg <- featureValues(ms1_data_neg, method="medret", value="into")
colnames(ms1_matrix_neg) <- mzml_names_neg
dim(ms1_matrix_neg)
feat_list_neg <- t(ms1_matrix_neg)

# Build feature summary
ms1_summary_neg <- featureSummary(ms1_data_neg)
ms1_def_neg <- featureDefinitions(ms1_data_neg)

# Transform data
feat_list_neg <- log2(feat_list_neg)

# Missing value imputation
#feat_list_neg[which(is.na(feat_list_neg))] <- median(na.omit(as.numeric(unlist(feat_list_neg))))
feat_list_neg <- scale(feat_list_neg, scale=FALSE, center=FALSE)
feat_list_neg[is.na(feat_list_neg)] <- 0
feat_list_neg[which(feat_list_neg < 0)] <- 0
feat_list_neg[is.infinite(feat_list_neg)] <- 0
feat_list_neg <- feat_list_neg[!apply(feat_list_neg, MARGIN=1, function(x) max(x,na.rm=TRUE) == min(x,na.rm=TRUE)),]

# Plot histogram
pdf(file="plots/neg_feat_list_hist.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
hist(as.numeric(feat_list_neg), main="Histogram of feature table")
dev.off()

# PCA
pdf(file="plots/neg_ms1_feature_table_pca.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
ms1_pca_neg <- prcomp(feat_list_neg, center=TRUE)
plot(ms1_pca_neg$x[, 1], ms1_pca_neg$x[, 2], pch=19, main="PCA of feature table",
	 xlab=paste0("PC1: ", format(summary(ms1_pca_neg)$importance[2, 1] * 100, digits=3), " % variance"),
	 ylab=paste0("PC2: ", format(summary(ms1_pca_neg)$importance[2, 2] * 100, digits=3), " % variance"),
	 col=mzml_pheno_colors_samples_neg, cex=2)
grid()
text(ms1_pca_neg$x[, 1], ms1_pca_neg$x[, 2], labels=ms1_data_neg$sample_name, col=mzml_pheno_colors_samples_neg, pos=3, cex=0.5)

plot(ms1_pca_neg$x[, 3], ms1_pca_neg$x[, 4], pch=19, main="PCA of feature table",
	 xlab=paste0("PC3: ", format(summary(ms1_pca_neg)$importance[2, 3] * 100, digits=3), " % variance"),
	 ylab=paste0("PC4: ", format(summary(ms1_pca_neg)$importance[2, 4] * 100, digits=3), " % variance"),
	 col=mzml_pheno_colors_samples_neg, cex=2)
grid()
text(ms1_pca_neg$x[, 3], ms1_pca_neg$x[, 4], labels=ms1_data_neg$sample_name, col=mzml_pheno_colors_samples_neg, pos=3, cex=0.5)

plot(ms1_pca_neg$x[, 5], ms1_pca_neg$x[, 6], pch=19, main="PCA of feature table",
	 xlab=paste0("PC5: ", format(summary(ms1_pca_neg)$importance[2, 5] * 100, digits=3), " % variance"),
	 ylab=paste0("PC6: ", format(summary(ms1_pca_neg)$importance[2, 6] * 100, digits=3), " % variance"),
	 col=mzml_pheno_colors_samples_neg, cex=2)
grid()
text(ms1_pca_neg$x[, 5], ms1_pca_neg$x[, 6], labels=ms1_data_neg$sample_name, col=mzml_pheno_colors_samples_neg, pos=3, cex=0.5)

plot(ms1_pca_neg$x[, 7], ms1_pca_neg$x[, 8], pch=19, main="PCA of feature table",
	 xlab=paste0("PC7: ", format(summary(ms1_pca_neg)$importance[2, 7] * 100, digits=3), " % variance"),
	 ylab=paste0("PC8: ", format(summary(ms1_pca_neg)$importance[2, 8] * 100, digits=3), " % variance"),
	 col=mzml_pheno_colors_samples_neg, cex=2)
grid()
text(ms1_pca_neg$x[, 7], ms1_pca_neg$x[, 8], labels=ms1_data_neg$sample_name, col=mzml_pheno_colors_samples_neg, pos=3, cex=0.5)

plot(ms1_pca_neg$x[, 9], ms1_pca_neg$x[, 10], pch=19, main="PCA of feature table",
	 xlab=paste0("PC9: ", format(summary(ms1_pca_neg)$importance[2, 9] * 100, digits=3), " % variance"),
	 ylab=paste0("PC10: ", format(summary(ms1_pca_neg)$importance[2, 10] * 100, digits=3), " % variance"),
	 col=mzml_pheno_colors_samples_neg, cex=2)
grid()
text(ms1_pca_neg$x[, 9], ms1_pca_neg$x[, 10], labels=ms1_data_neg$sample_name, col=mzml_pheno_colors_samples_neg, pos=3, cex=0.5)

screeplot(ms1_pca_neg, bstick=TRUE, type=c("barplot", "lines"), npcs=min(20, length(ms1_pca_neg$sdev)), ptype="o", bst.col="red", bst.lty="solid", xlab="Component", ylab="Inertia", main="Broken stick test of PCA", legend=FALSE)
dev.off()

# Create single 0/1 matrix
bina_list_neg <- t(ms1_matrix_neg)
bina_list_neg[is.na(bina_list_neg)] <- 1
bina_list_neg <- log2(bina_list_neg)
bina_list_neg[bina_list_neg < log2(ms1_intensity_cutoff)] <- 0
bina_list_neg[bina_list_neg != 0] <- 1

# Only unique compounds in group mzml_pheno$ and not the others
uniq_list_neg <- apply(X=bina_list_neg, MARGIN=2, FUN=function(x) { if (length(unique(mzml_pheno_neg$sample_group[grepl("1", x)])) == 1) x else rep(0, length(x)) } )
colnames(uniq_list_neg) <- colnames(bina_list_neg)
rownames(uniq_list_neg) <- rownames(bina_list_neg)

# Create data frame
model_div_neg             <- data.frame(features=apply(X=bina_list_neg, MARGIN=1, FUN=function(x) { sum(x) } ))
model_div_neg$richness    <- apply(X=bina_list_neg, MARGIN=1, FUN=function(x) { sum(x) } )
#model_div_neg$menhinick   <- apply(X=bina_list_neg, MARGIN=1, FUN=function(x) { menhinick.diversity(x) } )
model_div_neg$shannon     <- apply(X=feat_list_neg, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") })
model_div_neg$pielou      <- apply(X=scale(feat_list_neg, center=F), MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") / log(vegan::specnumber(x)) })
#model_div_neg$chao        <- vegan::specpool2vect(X=vegan::specpool(feat_list_neg, species), index="chao")
model_div_neg$simpson     <- apply(X=feat_list_neg, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="simpson") })
model_div_neg$inverse     <- apply(X=feat_list_neg, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="inv") })
model_div_neg$fisher      <- apply(X=feat_list_neg, MARGIN=1, FUN=function(x) { fisher.alpha(round(x,0)) })
model_div_neg$unique      <- apply(X=uniq_list_neg, MARGIN=1, FUN=function(x) { sum(x) })
model_div_neg$hillfunc    <- as.numeric(unlist(calcDiv(feat_list_neg, compDisMat=scales::rescale(as.matrix(dist(t(feat_list_neg)), diag=TRUE, upper=TRUE)), q=1, type="FuncHillDiv")))

# Remove NAs if present
model_div_neg[is.na(model_div_neg)] <- 0



# ---------- QC compounds ----------
# Check whether QC compounds have been picked by XCMS
QC_neg <- NULL
QC_neg <- rbind(QC_neg, data.frame(compound="N-(3-Indolylacetyl)-L-Valine", rt=462, mz=273.1218))
QC_neg <- rbind(QC_neg, data.frame(compound="Biochanin A", rt=585, mz=283.0587))
#QC_neg <- rbind(QC_neg, data.frame(compound="Kinetin", rt=279, mz=216.086))
QC_neg$FT <- ""
QC_neg$median_into <- 0
QC_neg$sd_into <- 0

for (i in 1:nrow(QC_neg)) {
	QC_neg$FT[i] <- rownames(ms1_def_neg[ms1_def_neg$rtmed>=QC_neg$rt[i]-rtabs & ms1_def_neg$rtmed<=QC_neg$rt[i]+rtabs & ms1_def_neg$mzmed>=QC_neg$mz[i]-mzabs & ms1_def_neg$mzmed<=QC_neg$mz[i]+mzabs, ])[1]
	QC_neg$median_into[i] <- median(ms1_matrix_neg[which(rownames(ms1_matrix_neg)==QC_neg$FT[i]), ], na.rm=TRUE)
	QC_neg$sd_into[i] <- sd(ms1_matrix_neg[which(rownames(ms1_matrix_neg)==QC_neg$FT[i]), ], na.rm=TRUE)
}

# Plot variation of MM8 compounds in the samples
pdf(file="plots/neg_ms1_qc_compounds_into.pdf", encoding="ISOLatin1", pointsize=12, width=5, height=5, family="Helvetica")
par(mar=c(7,4,1,1))
boxplot(t(ms1_matrix_neg[rownames(ms1_matrix_neg) %in% QC_neg$FT,]), main="", xlab="", ylab="log(into)", names=NA)
text(1:nrow(QC_neg), par("usr")[3]-(par("usr")[4]-par("usr")[3])/24, srt=-40, adj=c(0,1), labels=QC_neg$compound, xpd=TRUE, cex=0.8)
dev.off()

# Plot extracted ion chromatograms
chromas_eic_neg <- chromatogram(ms1_data_neg, aggregationFun="sum", adjustedRtime=TRUE, include="any", msLevel=1L)

# Plot chromatograms of QC compounds
QC_neg[QC_neg$compound=="Biochanin A", "rt"] <- 588
pdf(file="plots/neg_ms1_qc_compounds_chromatogram.pdf", encoding="ISOLatin1", pointsize=10, width=3, height=9, family="Helvetica")
par(mfrow=c(3,1), mar=c(4,4,3,1), cex=0.8)
for (j in 1:nrow(QC_neg)) {
	ymax = NULL
	for (i in 1:ncol(ms1_matrix_neg)) ymax = c(ymax, chromas_eic_neg@.Data[[i]]@intensity[chromas_eic_neg@.Data[[i]]@rtime>=QC_neg$rt[j]-rtabs & chromas_eic_neg@.Data[[i]]@rtime<=QC_neg$rt[j]+rtabs])
	#for (i in 1:ncol(ms1_matrix_neg)) ymax = c(ymax, chromas_tic_neg@.Data[[i]]@intensity[chromas_tic_neg@.Data[[i]]@rtime>=QC_neg$rt[j]-rtabs & chromas_tic_neg@.Data[[i]]@rtime<=QC_neg$rt[j]+rtabs])
	#for (i in 1:ncol(ms1_matrix_neg)) ymax = c(ymax, chromas_bpc_neg@.Data[[i]]@intensity[chromas_bpc_neg@.Data[[i]]@rtime>=QC_neg$rt[j]-rtabs & chromas_bpc_neg@.Data[[i]]@rtime<=QC_neg$rt[j]+rtabs])
	ymin = floor(min(ymax))
	ymax = ceiling(max(ymax))
	plot(0, 0, xlim=c(QC_neg[j,"rt"]-rtabs,QC_neg[j,"rt"]+rtabs), ylim=c(ymin,ymax), main=QC_neg[j,"compound"], xlab="RT [s]", ylab="log(into)")
	for (i in 1:ncol(ms1_matrix_neg)) {
		lines(chromas_eic_neg@.Data[[i]]@rtime, chromas_eic_neg@.Data[[i]]@intensity, lwd=2, col=mzml_pheno_colors_samples_neg[i])
		#lines(chromas_tic_neg@.Data[[i]]@rtime, chromas_tic_neg@.Data[[i]]@intensity, lwd=2, col=mzml_pheno_colors_samples_neg[i])
		#lines(chromas_bpc_neg@.Data[[i]]@rtime, chromas_bpc_neg@.Data[[i]]@intensity, lwd=2, col=mzml_pheno_colors_samples_neg[i])
	}
}
dev.off()



# ############################## MS2 ##############################



# ---------- MS2 spectra detection ----------
# Estimate precursor intensity
precursor_intensity_neg <- xcms::estimatePrecursorIntensity(ms1_data_neg)
print(head(na.omit(precursor_intensity_neg)))

# Reconstruct MS2 spectra from MS1 data
ms2_data_neg <- chromPeakSpectra(ms1_data_neg, msLevel=2L, return.type="Spectra")
print(ms2_data_neg)
print(length(ms2_data_neg$peak_id))

# Extract all usable MS2 spectra
ms2_spectra_neg <- list()
for (i in 1:nrow(ms1_def_neg)) {
#ms2_spectra_neg <- foreach(i=1:nrow(ms1_def_neg)) %dopar% {
	# Extract existing MS2 spectra for feature
	feature_of_interest <- ms1_def_neg[i, "mzmed"]
	peaks_of_interest <- chromPeaks(ms1_data_neg, mz=feature_of_interest, ppm=ppm)
	spectra_of_interest <- ms2_data_neg[ms2_data_neg$peak_id %in% rownames(peaks_of_interest)]
	
	# Continue if feature has MS2 peaks
	if (length(which(ms2_data_neg$peak_id %in% rownames(peaks_of_interest))) > 0) {
		# Filter for retention time (based on parameter rtabs)
		#combined_spectra_of_interest <- filterIntensity(spectra_of_interest, intensity=c(1,Inf), backend=MsBackendDataFrame)
		combined_spectra_of_interest <- Spectra::filterRt(spectra_of_interest, rt=c(ms1_def_neg[i, "rtmin"]-(rtabs/2),ms1_def_neg[i, "rtmax"]+(rtabs/2)))
		if (length(combined_spectra_of_interest) > 0) {
			combined_spectra_of_interest <- setBackend(combined_spectra_of_interest, backend=MsBackendDataFrame())
			
			# Combine spectra
			combined_spectra_of_interest <- Spectra::combineSpectra(combined_spectra_of_interest, FUN=combinePeaks, ppm=ppm, peaks="union", minProp=0.8, intensityFun=median, mzFun=median, backend=MsBackendDataFrame)#f=rownames(peaks_of_interest))
			
			# Remove noise from spectra
			#combined_spectra_of_interest <- pickPeaks(combined_spectra_of_interest, snr=1.0, method="SuperSmoother") #MAD
			#combined_spectra_of_interest <- Spectra::smooth(combined_spectra_of_interest, method="SavitzkyGolay") #(Weighted)MovingAverage
			
			# Plot merged spectrum
			#Spectra::plotSpectra(combined_spectra_of_interest)
			#combined_spectra_peaks <- as.data.frame(Spectra::peaksData(combined_spectra_of_interest)[[1]])
			#plot(x=combined_spectra_peaks[,1], y=combined_spectra_peaks[,2], type="h", xlab="m/z", ylab="intensity", main=paste("Precursor m/z",combined_spectra_of_interest@backend@spectraData$precursorMz[[1]]))
		}
	} else {
		combined_spectra_of_interest <- spectra_of_interest
	}
	
	ms2_spectra_neg[[i]] <- combined_spectra_of_interest
	#return(combined_spectra_of_interest)
}

# Remove empty spectra
names(ms2_spectra_neg) <- rownames(ms1_def_neg)
ms2_spectra_neg <- ms2_spectra_neg[lengths(ms2_spectra_neg) != 0]

# Relate all MS2 spectra to MS1 precursors
ms1_def_neg$has_ms2 <- as.integer(rownames(ms1_def_neg) %in% names(ms2_spectra_neg))
print(paste0("Number of MS2 spectra related to precursor: ", length(which(ms1_def_neg$has_ms2>0))))



# ---------- Annotate MS2 spectra ----------
# Save all MS2 spectra in MSP file
msp_text <- NULL
alignment_id <- 0
for (i in names(ms2_spectra_neg)) {
	alignment_id <- alignment_id + 1
	
	msp_text <- c(msp_text, paste("NAME:", i))
	msp_text <- c(msp_text, paste("AlignmentID:", alignment_id))
	msp_text <- c(msp_text, paste("RETENTIONTIME:", ms1_def_neg[i, "rtmed"]))
	msp_text <- c(msp_text, paste("PRECURSORMZ:", ms1_def_neg[i, "mzmed"]))
	msp_text <- c(msp_text, paste("METABOLITENAME:", i))
	if (polarity == "positive") {
		msp_text <- c(msp_text, paste("ADDUCTIONNAME:", "[M+H]+"))
	} else {
		msp_text <- c(msp_text, paste("ADDUCTIONNAME:", "[M-H]-"))
	}
	msp_text <- c(msp_text, paste("NumPeaks:", nrow(as.data.frame(peaksData(ms2_spectra_neg[[i]])[[1]]))))
	msp_text <- c(msp_text, paste(as.data.frame(peaksData(ms2_spectra_neg[[i]])[[1]])$mz, as.data.frame(peaksData(ms2_spectra_neg[[i]])[[1]])$intensity, sep="\t"))
	msp_text <- c(msp_text, "")
	
}

# Write MSP file
cat(msp_text, file="data/neg_ms2_spectra.msp", sep="\n")


# Save all MS2 spectra in MGF file
mgf_text <- NULL
for (i in names(ms2_spectra_neg)) {
	mgf_text <- c(mgf_text, paste0("COM=", i))
	mgf_text <- c(mgf_text, "BEGIN IONS")
	mgf_text <- c(mgf_text, "MSLEVEL=2")
	mgf_text <- c(mgf_text, paste0("TITLE=", i))
	mgf_text <- c(mgf_text, paste0("RTINSECONDS=", ms1_def_neg[i, "rtmed"]))
	mgf_text <- c(mgf_text, paste0("PEPMASS=", ms1_def_neg[i, "mzmed"]))
	if (polarity == "positive") {
		mgf_text <- c(mgf_text, paste0("CHARGE=", "1+"))
	} else {
		mgf_text <- c(mgf_text, paste0("CHARGE=", "1-"))
	}
	mgf_text <- c(mgf_text, paste(as.data.frame(peaksData(ms2_spectra_neg[[i]])[[1]])$mz, as.data.frame(peaksData(ms2_spectra_neg[[i]])[[1]])$intensity, sep=" "))
	mgf_text <- c(mgf_text, "END IONS")
	mgf_text <- c(mgf_text, "")
}

# Write MGF file
cat(mgf_text, file="data/neg_ms2_spectra.mgf", sep="\n")



# ---------- Annotate MS2 spectra with SIRIUS ----------
# Apply annotated compounds onto feature table
ms1_def_neg$has_id <- 0
ms1_def_neg$smiles <- ""
ms1_def_neg$inchi <- ""
ms1_def_neg$inchikey <- ""
ms1_def_neg$name <- ""
ms1_def_neg$classyfire <- ""

# Read SIRIUS annotation
annotator_sirius_neg <- read.table(file=paste0("data/neg_ms2_compound_identifications.tsv"), header=TRUE, sep="\t", quote="\"", comment.char="", fill=TRUE, dec=".", stringsAsFactors=FALSE)#, encoding="UTF-8", fileEncoding="UTF-8")
annotator_sirius_neg$Metabolite.name <- gsub(x=annotator_sirius_neg$id, pattern=".*_", replacement="")

for (j in unique(annotator_sirius_neg$Metabolite.name)) {
	ms1_def_neg[j, "has_id"] <- 1
	
	obj <- annotator_sirius_neg[annotator_sirius_neg$Metabolite.name %in% j, "smiles"]
	ms1_def_neg[j, "smiles"] <- obj[1]
	
	obj <- annotator_sirius_neg[annotator_sirius_neg$Metabolite.name %in% j, "InChI"]
	ms1_def_neg[j, "inchi"] <- obj[1]
	
	obj <- annotator_sirius_neg[annotator_sirius_neg$Metabolite.name %in% j, "name"]
	if (obj[1] != "null" ) {
		ms1_def_neg[j, "name"] <- obj[1]
	}
}

# Convert InChi to InChiKey
for (i in 1:nrow(ms1_def_neg)) {
	if ((ms1_def_neg[i, "smiles"] != "") & (ms1_def_neg[i, "inchikey"] == "")) {
		ms1_def_neg[i, "inchikey"] <- rinchi::get.inchi.key(ms1_def_neg[i, "smiles"])
	}
}

# ClassyFire of InChiKey
for (i in 1:nrow(ms1_def_neg)) {
	if ((ms1_def_neg[i, "inchikey"] != "") & (ms1_def_neg[i, "classyfire"] =="")) {
		ms1_def_neg[i, "classyfire"] <- paste(as.data.frame(get_classification(ms1_def_neg[i, "inchikey"])@classification)$Classification, collapse="; ")
	}
}

# Annotation
sirius_compounds_neg <- annotator_sirius_neg$smiles
sirius_compound_ids_neg <- gsub(x=annotator_sirius_neg$id, pattern=".*_", replacement="")

# Calculate molecular descriptors with CDK
cdk_descriptors_neg <- NULL
for (i in sirius_compounds_neg) {
	# Get Structure from SMILES
	cdk_mol = parse.smiles(i)[[1]]
	
	# Get simple measures
	cdk_atoms = get.atoms(cdk_mol)
	cdk_bonds = get.bonds(cdk_mol)
	
	# Calculate simple measures
	MolWeight = get.exact.mass(cdk_mol)
	nAtoms = get.atom.count(cdk_mol)
	cdk_num_atoms = as.factor(unlist(lapply(cdk_atoms, get.symbol)))
	cdk_num_atoms = tapply(cdk_num_atoms, cdk_num_atoms, length)
	numC = as.numeric(cdk_num_atoms["C"])
	numN = as.numeric(cdk_num_atoms["N"])
	numP = as.numeric(cdk_num_atoms["P"])
	numO = as.numeric(cdk_num_atoms["O"])
	numHydrogen = get.total.hydrogen.count(cdk_mol)
	CNRatio = as.numeric(numC / numN)
	XLogP = get.xlogp(cdk_mol)
	
	# Calculate descriptors and restrict to only "constitutional" and "topological"
	cdk_mol_des_cats = get.desc.categories()
	cdk_mol_des_names = c(get.desc.names(cdk_mol_des_cats[3]), get.desc.names(cdk_mol_des_cats[4]))
	cdk_mol_des = as.data.frame(eval.desc(cdk_mol, cdk_mol_des_names))
	cdk_descriptors_neg <- plyr::rbind.fill(cdk_descriptors_neg, cbind(data.frame(MolWeight=MolWeight, nAtoms=nAtoms, numC=numC, numN=numN, numP=numP, numO=numO, numHydrogen=numHydrogen, CNRatio=CNRatio, XLogP=XLogP), cdk_mol_des))
}
rownames(cdk_descriptors_neg) <- paste0(sirius_compound_ids_neg,"_",pol)

# Properly format NAs and convert to numeric
for (i in 1:nrow(cdk_descriptors_neg)) {
	cdk_descriptors_neg[i, as.character(cdk_descriptors_neg[i,]) %in% 'list(NULL)'] <- 0
	cdk_descriptors_neg[i, as.character(cdk_descriptors_neg[i,]) %in% 'NA'] <- 0
	cdk_descriptors_neg[i, as.character(cdk_descriptors_neg[i,]) %in% 'NULL'] <- 0
	cdk_descriptors_neg[i, as.character(cdk_descriptors_neg[i,]) %in% 'NaN'] <- 0
	cdk_descriptors_neg[i, ] <- as.numeric(unlist(cdk_descriptors_neg[i,]))
}

# Bug: Remove left-over variables
rm(CNRatio, MolWeight, XLogP, nAtoms, numC, numHydrogen, numN, numO, numP)



# ---------- Classify MS2 spectra with CANOPUS ----------
# Apply classified classes onto feature table
ms1_def_neg$primary_class <- ""
ms1_def_neg$alternative_classes <- ""
ms1_def_neg$NP_class <- ""

# Read SIRIUS/CANOPUS classifier
if (SIRIUS_VERSION == 4) {
	classifier_canopus_neg <- read.table(file=paste0("data/neg_ms2_canopus_compound_summary.tsv"), header=TRUE, sep="\t", quote="\"", fill=FALSE, dec=".", stringsAsFactors=FALSE)#, encoding="UTF-8", fileEncoding="UTF-8")
	classifier_canopus_neg$Metabolite.name <- gsub(x=classifier_canopus_neg$name, pattern=".*_", replacement="")
	classifier_canopus_neg$primary_class <- paste("Organic compounds", classifier_canopus_neg$superclass, classifier_canopus_neg$class, classifier_canopus_neg$subclass, classifier_canopus_neg$level.5, sep="; ")
	classifier_canopus_neg$primary_class <- gsub(x=classifier_canopus_neg$primary_class, pattern="; ; ", replacement="; ", perl=TRUE)
	classifier_canopus_neg$primary_class <- gsub(x=classifier_canopus_neg$primary_class, pattern="; ; ", replacement="; ", perl=TRUE)
	classifier_canopus_neg$primary_class <- gsub(x=classifier_canopus_neg$primary_class, pattern="; $", replacement="", perl=TRUE)
	classifier_canopus_neg$primary_class <- gsub(x=classifier_canopus_neg$primary_class, pattern="(\\'|\\>|\\(|\\))", replacement="", perl=TRUE)
	classifier_canopus_neg$Annotation..putative. <- classifier_canopus_neg$primary_class
	classifier_canopus_neg$alternative_classes <- classifier_canopus_neg$all.classifications
} else {
	classifier_canopus_neg <- read.table(file=paste0("data/neg_ms2_canopus_compound_summary.tsv"), header=TRUE, sep="\t", quote="\"", comment.char="", fill=TRUE, dec=".", stringsAsFactors=FALSE)#, encoding="UTF-8", fileEncoding="UTF-8")
	classifier_canopus_neg$Metabolite.name <- ""
	for (i in 1:length(classifier_canopus_neg$id)) {
		x = which(gsub(x=classifier_canopus_neg$id, pattern=".*?_", replacement="", perl=TRUE) %in% gsub(x=classifier_canopus_neg$id[i], pattern=".*?_", replacement="", perl=TRUE))
		if (length(x) > 0) classifier_canopus_neg$Metabolite.name[i] <- gsub(x=classifier_canopus_neg$id[x], pattern=".*?_", replacement="", perl=TRUE)
	}
	classifier_canopus_neg$Metabolite.name[classifier_canopus_neg$Metabolite.name == "null"] <- ""
	classifier_canopus_neg$primary_class <- paste("Organic compounds", classifier_canopus_neg$ClassyFire.superclass, classifier_canopus_neg$ClassyFire.class, classifier_canopus_neg$ClassyFire.subclass, classifier_canopus_neg$ClassyFire.level.5, sep="; ")
	classifier_canopus_neg$primary_class <- gsub(x=classifier_canopus_neg$primary_class, pattern="; ; ", replacement="; ", perl=TRUE)
	classifier_canopus_neg$primary_class <- gsub(x=classifier_canopus_neg$primary_class, pattern="; ; ", replacement="; ", perl=TRUE)
	classifier_canopus_neg$primary_class <- gsub(x=classifier_canopus_neg$primary_class, pattern="; $", replacement="", perl=TRUE)
	classifier_canopus_neg$primary_class <- gsub(x=classifier_canopus_neg$primary_class, pattern="(\\'|\\>|\\(|\\))", replacement="", perl=TRUE)
	classifier_canopus_neg$Annotation..putative. <- classifier_canopus_neg$primary_class
	classifier_canopus_neg$alternative_classes <- classifier_canopus_neg$all.classifications
}

for (j in unique(classifier_canopus_neg$Metabolite.name)) {
	ms1_def_neg[j, "has_id"] <- 1
	
	# CHEMONT
	obj <- classifier_canopus_neg[classifier_canopus_neg$Metabolite.name %in% j, "Annotation..putative."]
	ms1_def_neg[j, "primary_class"] <- obj[1]
	
	# NPClassifier
	if (SIRIUS_VERSION > 4) {
		obj <- paste(sep="; ", classifier_canopus_neg[classifier_canopus_neg$Metabolite.name %in% j, "NPC.pathway"],
					 classifier_canopus_neg[classifier_canopus_neg$Metabolite.name %in% j, "NPC.superclass"],
					 classifier_canopus_neg[classifier_canopus_neg$Metabolite.name %in% j, "NPC.class"])
		ms1_def_neg[j, "NP_class"] <- obj[1]
	}
}



# ---------- Diversity of MS2 classes ----------
# Create CANOPUS classifier object for each sample
classifiers_neg <- list()
for (i in mzml_names_neg) {
	obj <- names(which(bina_list_neg[rownames(bina_list_neg)==i, colnames(bina_list_neg) %in% classifier_canopus_neg$Metabolite.name] > 0))
	classifiers_neg[[i]] <- classifier_canopus_neg[classifier_canopus_neg$Metabolite.name %in% obj, ]
}

# Diversity of classes per sample
div_classes_samples_neg <- NULL
for (i in mzml_names_neg) {
	obj <- table(classifiers_neg[[i]][,"Annotation..putative."])
	obj <- data.frame(classes=names(obj), frequency=as.numeric(obj))
	if (is.null(div_classes_samples_neg)) {
		div_classes_samples_neg <- obj
	} else {
		div_classes_samples_neg <- merge(div_classes_samples_neg, obj, by="classes", all.x=TRUE, all.y=TRUE)
	}
}
rownames(div_classes_samples_neg) <- div_classes_samples_neg$classes
div_classes_samples_neg <- div_classes_samples_neg[, -which(colnames(div_classes_samples_neg)=="classes")]
colnames(div_classes_samples_neg) <- mzml_names_neg

# Diversity of classes
div_classes_neg <- div_classes_samples_neg
div_classes_neg[is.na(div_classes_neg)] <- 0
div_classes_neg <- apply(X=div_classes_neg, MARGIN=1, FUN=function(x) { sum(x) })
div_classes_neg <- data.frame(row.names=names(div_classes_neg), frequency=as.numeric(div_classes_neg))

# Plot diversity of classes in all samples
pdf(file="plots/neg_ms2_classes_diversity.pdf", encoding="ISOLatin1", pointsize=8, width=6, height=14, family="Helvetica")
par(mfrow=c(1,1), mar=c(4,15,4,1), oma=c(0,0,0,0), cex.axis=0.8, cex=0.9)
barplot(div_classes_neg$frequency, names.arg=gsub('.*; ','',rownames(div_classes_neg)), las=1, horiz=TRUE, xlab="frequency", main="Diversity of compound classes", col=rainbow(n=nrow(div_classes_neg), alpha=0.6))
dev.off()

# Sunburst plot of classes of all samples
pdf(file="plots/neg_ms2_classes_sunburst.pdf", encoding="ISOLatin1", pointsize=8, width=10, height=10, family="Helvetica")
par(mfrow=c(1,1), mar=c(0,0,0,0), oma=c(0,0,0,0), cex.axis=1, cex=1)
sunBurstPlotFromSubstanceClasses(rownames(div_classes_neg), div_classes_neg$frequency, colorStart=0.0, colorAlpha=0.6)
dev.off()
write.csv(data.frame("Compound classes"=rownames(div_classes_neg), "Number of spectra"=div_classes_neg$frequency), file="plots/neg_ms2_classes_sunburst.csv", row.names=FALSE)

# Classes
classes_neg <- rownames(div_classes_neg)
classes_neg <- classes_neg[which(grepl(pattern="^Organic compounds", x=classes_neg))]
classes_neg <- gsub(x=classes_neg, pattern=":", replacement="")
classes_neg <- gsub(x=classes_neg, pattern="/", replacement="; ")



# ---------- Build diversity objects ----------
# Imputation of NA with zeros
div_classes_neg[is.na(div_classes_neg)] <- 0
div_classes_samples_neg[is.na(div_classes_samples_neg)] <- 0

# Classification list for statistics
class_list_neg <- as.data.frame(t(div_classes_samples_neg))
class_list_neg[is.na(class_list_neg)] <- 0

# Log Transform
#class_list_neg <- log2(class_list_neg + 1)

# Only keep class names
colnames(class_list_neg) <- gsub(x=colnames(class_list_neg), pattern='.*; ', replacement='')

# Generate class_int_list_neg with abundances instead of counts
class_int_list_neg <- class_list_neg

for (i in 1:nrow(class_list_neg)) {
	samp <- rownames(class_list_neg)[i]
	for (j in 1:ncol(class_list_neg)) {
		cl <- colnames(class_list_neg)[j]
		ft <- classifier_canopus_neg$Metabolite.name[which(gsub(x=classifier_canopus_neg$primary_class, pattern='.*; ', replacement='') == cl)]
		ints <- as.numeric(feat_list_neg[i, which(colnames(feat_list_neg) %in% ft)])
		class_int_list_neg[i, j] <- sum(ints)
	}
}



# ---------- Classification at CHEMONT level of classes ----------
# Make superclasses at CHEMONT level 2
superclass_level_neg <- 2
div_superclasses_samples_names_neg <- NULL
for (i in c(1:superclass_level_neg)) {
	div_superclasses_samples_names_neg <- c(div_superclasses_samples_names_neg, lapply(X=strsplit(rownames(div_classes_samples_neg), '; '), FUN=function(x) { gsub(x=paste(x[1:i],sep='',collapse='; '),pattern='; NA',replacement='') }))
}
div_superclasses_samples_neg <- data.frame()
for (i in c(1:ncol(div_classes_samples_neg))) div_superclasses_samples_neg <- rbind(div_superclasses_samples_neg, rep(0, length(unique(div_superclasses_samples_names_neg))))
div_superclasses_samples_neg <- t(div_superclasses_samples_neg)
colnames(div_superclasses_samples_neg) <- colnames(div_classes_samples_neg)
rownames(div_superclasses_samples_neg) <- unique(div_superclasses_samples_names_neg)
for (i in rownames(div_superclasses_samples_neg)) {
	for (j in c(1:ncol(div_classes_samples_neg))) {
		div_superclasses_samples_neg[rownames(div_superclasses_samples_neg)==i, j] <- sum(div_classes_samples_neg[grep(x=rownames(div_classes_samples_neg), pattern=i), j])
	}
}

# Diversity of superclasses
div_superclasses_neg <- div_superclasses_samples_neg
div_superclasses_neg[is.na(div_superclasses_neg)] <- 0
div_superclasses_neg <- apply(X=div_superclasses_neg, MARGIN=1, FUN=function(x) { sum(x) })
div_superclasses_neg <- data.frame(row.names=names(div_superclasses_neg), frequency=as.numeric(div_superclasses_neg))

# Sunburst plot
pdf(file="plots/neg_ms2_superclasses_sunburst.pdf", encoding="ISOLatin1", pointsize=8, width=10, height=10, family="Helvetica")
sunBurstPlotFromSubstanceClasses(rownames(div_superclasses_neg), div_superclasses_neg$frequency, colorStart=0.0, colorAlpha=0.6)
dev.off()
write.csv(data.frame("Compound classes"=rownames(div_superclasses_neg), "Number of spectra"=div_superclasses_neg$frequency), file="plots/neg_ms2_superclasses_sunburst.csv", row.names=FALSE)

# Imputation of NA with zeros
div_superclasses_neg[is.na(div_superclasses_neg)] <- 0
div_superclasses_samples_neg[is.na(div_superclasses_samples_neg)] <- 0

# Classification list for statistics
superclass_list_neg <- as.data.frame(t(div_superclasses_samples_neg))
superclass_list_neg[is.na(superclass_list_neg)] <- 0

# Log Transform
#superclass_list_neg <- log2(superclass_list_neg + 1)

# Only keep superclass names
colnames(superclass_list_neg) <- gsub(x=colnames(superclass_list_neg), pattern='.*; ', replacement='')

# Generate superclass_int_list_neg with abundances instead of counts
superclass_int_list_neg <- superclass_list_neg

for (i in 1:nrow(superclass_list_neg)) {
	samp <- rownames(superclass_list_neg)[i]
	for (j in 1:ncol(superclass_list_neg)) {
		cl <- colnames(superclass_list_neg)[j]
		ft <- classifier_canopus_neg$Metabolite.name[which(classifier_canopus_neg$primary_class %in% classifier_canopus_neg$primary_class[grep(x=classifier_canopus_neg$primary_class, pattern=cl)]) ]
		ints <- as.numeric(feat_list_neg[i, which(colnames(feat_list_neg) %in% ft)])
		superclass_int_list_neg[i, j] <- sum(ints)
	}
}



# ---------- Classification at CHEMONT level of subclasses ----------
# Make subclasses at CHEMONT level 3
subclass_level_neg <- 3
div_subclasses_samples_names_neg <- NULL
for (i in c(1:subclass_level_neg)) {
	div_subclasses_samples_names_neg <- c(div_subclasses_samples_names_neg, lapply(X=strsplit(rownames(div_classes_samples_neg), '; '), FUN=function(x) { gsub(x=paste(x[1:i],sep='',collapse='; '),pattern='; NA',replacement='') }))
}
div_subclasses_samples_neg <- data.frame()
for (i in c(1:ncol(div_classes_samples_neg))) div_subclasses_samples_neg <- rbind(div_subclasses_samples_neg, rep(0, length(unique(div_subclasses_samples_names_neg))))
div_subclasses_samples_neg <- t(div_subclasses_samples_neg)
colnames(div_subclasses_samples_neg) <- colnames(div_classes_samples_neg)
rownames(div_subclasses_samples_neg) <- unique(div_subclasses_samples_names_neg)
for (i in rownames(div_subclasses_samples_neg)) {
	for (j in c(1:ncol(div_classes_samples_neg))) {
		div_subclasses_samples_neg[rownames(div_subclasses_samples_neg)==i, j] <- sum(div_classes_samples_neg[grep(x=rownames(div_classes_samples_neg), pattern=i), j])
	}
}

# Diversity of subclasses
div_subclasses_neg <- div_subclasses_samples_neg
div_subclasses_neg[is.na(div_subclasses_neg)] <- 0
div_subclasses_neg <- apply(X=div_subclasses_neg, MARGIN=1, FUN=function(x) { sum(x) })
div_subclasses_neg <- data.frame(row.names=names(div_subclasses_neg), frequency=as.numeric(div_subclasses_neg))

# Sunburst plot
pdf(file="plots/neg_ms2_subclasses_sunburst.pdf", encoding="ISOLatin1", pointsize=8, width=10, height=10, family="Helvetica")
sunBurstPlotFromSubstanceClasses(rownames(div_subclasses_neg), div_subclasses_neg$frequency, colorStart=0.0, colorAlpha=0.6)
dev.off()
write.csv(data.frame("Compound classes"=rownames(div_subclasses_neg), "Number of spectra"=div_subclasses_neg$frequency), file="plots/neg_ms2_subclasses_sunburst.csv", row.names=FALSE)

# Imputation of NA with zeros
div_subclasses_neg[is.na(div_subclasses_neg)] <- 0
div_subclasses_samples_neg[is.na(div_subclasses_samples_neg)] <- 0

# Classification list for statistics
subclass_list_neg <- as.data.frame(t(div_subclasses_samples_neg))
subclass_list_neg[is.na(subclass_list_neg)] <- 0

# Log Transform
#subclass_list_neg <- log2(subclass_list_neg + 1)

# Only keep subclass names
colnames(subclass_list_neg) <- gsub(x=colnames(subclass_list_neg), pattern='.*; ', replacement='')

# Generate subclass_int_list_neg with abundances instead of counts
subclass_int_list_neg <- subclass_list_neg

for (i in 1:nrow(subclass_list_neg)) {
	samp <- rownames(subclass_list_neg)[i]
	for (j in 1:ncol(subclass_list_neg)) {
		cl <- colnames(subclass_list_neg)[j]
		ft <- classifier_canopus_neg$Metabolite.name[which(classifier_canopus_neg$primary_class %in% classifier_canopus_neg$primary_class[grep(x=classifier_canopus_neg$primary_class, pattern=cl)]) ]
		ints <- as.numeric(feat_list_neg[i, which(colnames(feat_list_neg) %in% ft)])
		subclass_int_list_neg[i, j] <- sum(ints)
	}
}



# ---------- Diversity of NP classes ----------
# Create CANOPUS classifier object for each sample
npclassifiers_neg <- list()
for (i in mzml_names_neg) {
	obj <- names(which(bina_list_neg[rownames(bina_list_neg)==i, colnames(bina_list_neg) %in% classifier_canopus_neg$Metabolite.name] > 0))
	npclassifiers_neg[[i]] <- classifier_canopus_neg[classifier_canopus_neg$Metabolite.name %in% obj, ]
}

# Diversity of NPclasses per sample
div_npclasses_samples_neg <- NULL
for (i in mzml_names_neg) {
	obj <- gsub(x=apply(X=npclassifiers_neg[[i]][,c("NPC.pathway","NPC.superclass","NPC.class")], MARGIN=1, FUN=function(x) { paste(x, collapse="; ") }), pattern="(\\(|\\))", replacement="", perl=TRUE)
	obj <- table(obj)
	obj <- data.frame(classes=names(obj), frequency=as.numeric(obj))
	if (is.null(div_npclasses_samples_neg)) {
		div_npclasses_samples_neg <- obj
	} else {
		div_npclasses_samples_neg <- merge(div_npclasses_samples_neg, obj, by="classes", all.x=TRUE, all.y=TRUE)
	}
}
rownames(div_npclasses_samples_neg) <- div_npclasses_samples_neg$classes
div_npclasses_samples_neg <- div_npclasses_samples_neg[, -which(colnames(div_npclasses_samples_neg)=="classes")]
colnames(div_npclasses_samples_neg) <- mzml_names_neg

# Diversity of NPclasses
div_npclasses_neg <- div_npclasses_samples_neg
div_npclasses_neg[is.na(div_npclasses_neg)] <- 0
div_npclasses_neg <- apply(X=div_npclasses_neg, MARGIN=1, FUN=function(x) { sum(x) })
div_npclasses_neg <- data.frame(row.names=names(div_npclasses_neg), frequency=as.numeric(div_npclasses_neg))

# Plot diversity of NPclasses in all samples
pdf(file="plots/neg_ms2_npclasses_diversity.pdf", encoding="ISOLatin1", pointsize=8, width=6, height=14, family="Helvetica")
par(mfrow=c(1,1), mar=c(4,15,4,1), oma=c(0,0,0,0), cex.axis=0.8, cex=0.9)
barplot(div_npclasses_neg$frequency, names.arg=gsub('.*; ','',rownames(div_npclasses_neg)), las=1, horiz=TRUE, xlab="frequency", main="Diversity of compound classes", col=rainbow(n=nrow(div_npclasses_neg), alpha=0.6))
dev.off()

# Sunburst plot of NPclasses of all samples
pdf(file="plots/neg_ms2_npclasses_sunburst.pdf", encoding="ISOLatin1", pointsize=8, width=10, height=10, family="Helvetica")
par(mfrow=c(1,1), mar=c(0,0,0,0), oma=c(0,0,0,0), cex.axis=1, cex=1)
sunBurstPlotFromSubstanceClasses(rownames(div_npclasses_neg), div_npclasses_neg$frequency, colorStart=0.0, colorAlpha=0.6)
dev.off()
write.csv(data.frame("Compound classes"=rownames(div_npclasses_neg), "Number of spectra"=div_npclasses_neg$frequency), file="plots/neg_ms2_npclasses_sunburst.csv", row.names=FALSE)

# NPClasses
npclasses_neg <- rownames(div_npclasses_neg)



# ---------- Build NPclass diversity objects ----------
# Imputation of NA with zeros
div_npclasses_neg[is.na(div_npclasses_neg)] <- 0
div_npclasses_samples_neg[is.na(div_npclasses_samples_neg)] <- 0

# Classification list for statistics
npclass_list_neg <- as.data.frame(t(div_npclasses_samples_neg))
npclass_list_neg[is.na(npclass_list_neg)] <- 0

# Log Transform
#npclass_list_neg <- log2(npclass_list_neg + 1)

# Only keep class names
colnames(npclass_list_neg) <- gsub(x=colnames(npclass_list_neg), pattern='.*; ', replacement='')

# Generate npclass_int_list_neg with abundances instead of counts
npclass_int_list_neg <- npclass_list_neg

for (i in 1:nrow(npclass_list_neg)) {
	samp <- rownames(npclass_list_neg)[i]
	for (j in 1:ncol(npclass_list_neg)) {
		cl <- colnames(npclass_list_neg)[j]
		ft <- classifier_canopus_neg$Metabolite.name[which(gsub(x=classifier_canopus_neg$NPC.class, pattern='.*; ', replacement='') == cl)]
		ints <- as.numeric(feat_list_neg[i, which(colnames(feat_list_neg) %in% ft)])
		npclass_int_list_neg[i, j] <- sum(ints)
	}
}



# ---------- Classification at NP-Pathways ----------
# Extract NP-Pathways
div_nppathway_samples_names_neg <- NULL
div_nppathway_samples_names_neg <- c(div_nppathway_samples_names_neg, lapply(X=strsplit(rownames(div_npclasses_samples_neg), '; '), FUN=function(x) { gsub(x=paste(x[1],sep='',collapse='; '),pattern='; NA',replacement='') }))

div_nppathway_samples_neg <- data.frame()
for (i in c(1:ncol(div_npclasses_samples_neg))) div_nppathway_samples_neg <- rbind(div_nppathway_samples_neg, rep(0, length(unique(div_nppathway_samples_names_neg))))
div_nppathway_samples_neg <- t(div_nppathway_samples_neg)
colnames(div_nppathway_samples_neg) <- colnames(div_npclasses_samples_neg)
rownames(div_nppathway_samples_neg) <- unique(div_nppathway_samples_names_neg)
for (i in rownames(div_nppathway_samples_neg)) {
	for (j in c(1:ncol(div_npclasses_samples_neg))) {
		div_nppathway_samples_neg[rownames(div_nppathway_samples_neg)==i, j] <- sum(div_npclasses_samples_neg[grep(x=rownames(div_npclasses_samples_neg), pattern=i), j])
	}
}

# Diversity of nppathway
div_nppathway_neg <- div_nppathway_samples_neg
div_nppathway_neg[is.na(div_nppathway_neg)] <- 0
div_nppathway_neg <- apply(X=div_nppathway_neg, MARGIN=1, FUN=function(x) { sum(x) })
div_nppathway_neg <- data.frame(row.names=names(div_nppathway_neg), frequency=as.numeric(div_nppathway_neg))

# Sunburst plot
pdf(file="plots/neg_ms2_nppathway_sunburst.pdf", encoding="ISOLatin1", pointsize=8, width=10, height=10, family="Helvetica")
sunBurstPlotFromSubstanceClasses(rownames(div_nppathway_neg), div_nppathway_neg$frequency, colorStart=0.0, colorAlpha=0.6)
dev.off()
write.csv(data.frame("Compound classes"=rownames(div_nppathway_neg), "Number of spectra"=div_nppathway_neg$frequency), file="plots/neg_ms2_nppathway_sunburst.csv", row.names=FALSE)

# Imputation of NA with zeros
div_nppathway_neg[is.na(div_nppathway_neg)] <- 0
div_nppathway_samples_neg[is.na(div_nppathway_samples_neg)] <- 0

# Classification list for statistics
nppathway_list_neg <- as.data.frame(t(div_nppathway_samples_neg))
nppathway_list_neg[is.na(nppathway_list_neg)] <- 0

# Log Transform
#nppathway_list_neg <- log2(nppathway_list_neg + 1)

# Only keep nppathway names
colnames(nppathway_list_neg) <- gsub(x=colnames(nppathway_list_neg), pattern='.*; ', replacement='')

# Generate nppathway_int_list_neg with abundances instead of counts
nppathway_int_list_neg <- nppathway_list_neg

for (i in 1:nrow(nppathway_list_neg)) {
	samp <- rownames(nppathway_list_neg)[i]
	for (j in 1:ncol(nppathway_list_neg)) {
		cl <- colnames(nppathway_list_neg)[j]
		ft <- classifier_canopus_neg$Metabolite.name[which(classifier_canopus_neg$NPC.pathway %in% classifier_canopus_neg$NPC.pathway[grep(x=classifier_canopus_neg$NPC.pathway, pattern=cl)]) ]
		ints <- as.numeric(feat_list_neg[i, which(colnames(feat_list_neg) %in% ft)])
		nppathway_int_list_neg[i, j] <- sum(ints)
	}
}


