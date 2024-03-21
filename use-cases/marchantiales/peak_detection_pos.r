


# ############################## MS1 ##############################



# ---------- MS1 Preparations ----------
# MS1 variables
polarity <- "positive"
pol <- substr(x=polarity, start=1, stop=3)

# General variables
mzml_files_pos <- NULL
mzml_names_pos <- NULL
mzml_times_pos <- NULL

# Load files
mzml_files_pos <- list.files(mzml_dir, pattern="*.mzML", recursive=T, full.names=T)
mzml_files_pos <- mzml_files_pos[grep(pol, mzml_files_pos, invert=FALSE)]

# Basenames of files without path and without extension
mzml_names_pos <- gsub('(.*)\\..*', '\\1', gsub('( |-|,)', '.', basename(mzml_files_pos)))

mzml_names_pos <- gsub(x=mzml_names_pos, pattern="^1\\.", replacement="R.cavernosa.SWE.", perl=TRUE)
mzml_names_pos <- gsub(x=mzml_names_pos, pattern="^11\\.", replacement="R.canaliculata.GOT.", perl=TRUE)
mzml_names_pos <- gsub(x=mzml_names_pos, pattern="^12\\.", replacement="R.bifurca.GOT.", perl=TRUE)
mzml_names_pos <- gsub(x=mzml_names_pos, pattern="^13\\.", replacement="R.ciliifera.GOT.", perl=TRUE)
mzml_names_pos <- gsub(x=mzml_names_pos, pattern="^17\\.", replacement="R.subbifurca.OE1.", perl=TRUE)
mzml_names_pos <- gsub(x=mzml_names_pos, pattern="^2\\.", replacement="R.sorocarpa.SWE.", perl=TRUE)
mzml_names_pos <- gsub(x=mzml_names_pos, pattern="^20\\.", replacement="R.beyrichiana.OEL.", perl=TRUE)
mzml_names_pos <- gsub(x=mzml_names_pos, pattern="^21\\.", replacement="R.subbifurca.OE2.", perl=TRUE)
mzml_names_pos <- gsub(x=mzml_names_pos, pattern="^25\\.", replacement="R.ciliifera.HAL.", perl=TRUE)
mzml_names_pos <- gsub(x=mzml_names_pos, pattern="^26\\.", replacement="R.gougetiana.HAL.", perl=TRUE)
mzml_names_pos <- gsub(x=mzml_names_pos, pattern="^45\\.", replacement="R.huebeneriana.GOT.", perl=TRUE)
mzml_names_pos <- gsub(x=mzml_names_pos, pattern="^5\\.", replacement="A.gracilis.SWE.", perl=TRUE)
mzml_names_pos <- gsub(x=mzml_names_pos, pattern="^7\\.", replacement="R.gothica.GOT.", perl=TRUE)
mzml_names_pos <- gsub(x=mzml_names_pos, pattern="^8\\.", replacement="M.fragrans.GOT.", perl=TRUE)
mzml_names_pos <- gsub(x=mzml_names_pos, pattern="^8c\\.", replacement="R.hemisphaerica.GOT.", perl=TRUE)
mzml_names_pos <- gsub(x=mzml_names_pos, pattern="^9\\.", replacement="A.hyalina.GOT.", perl=TRUE)

# Create phenodata based on species
mzml_pheno_pos <- data.frame(sample_name=mzml_names_pos, sample_group=as.factor(gsub(pattern="\\.\\d.*", replacement="", x=mzml_names_pos, perl=TRUE)))
mzml_pheno_pos$location <- as.factor(gsub(x=gsub(x=mzml_pheno_pos$sample_group, pattern=".*\\.", replacement=""), pattern="\\d", replacement="L"))
mzml_pheno_pos$species <- as.factor(gsub(x=mzml_pheno_pos$sample_group, pattern="\\.([^\\.]*)$", replacement=""))
mzml_pheno_samples_pos <- as.factor(mzml_pheno_pos$sample_group)
#mzml_pheno_colors_pos <- c("lemonchiffon4","orange4","mediumaquamarine","palegreen4","limegreen","royalblue","turquoise3","limegreen","palegreen4","seagreen4","tomato4","orchid3","green3","magenta3","violetred","mediumpurple3")
mzml_pheno_colors_pos <- c("orchid3","violetred","magenta3","turquoise3","mediumaquamarine","orange4","lemonchiffon4","palegreen4","steelblue3","seagreen4","mediumpurple3","tomato4","royalblue","limegreen")
#mzml_pheno_colors_pos <- RColorBrewer::brewer.pal(n=nlevels(mzml_pheno_samples_pos), name="Set1")
mzml_pheno_colors_samples_pos <- sapply(mzml_pheno_samples_pos, function(x) { x <- mzml_pheno_colors_pos[which(x==levels(mzml_pheno_samples_pos))] } )

# Save timestamps of samples
for (i in 1:length(mzml_files_pos)) {
	fl <- mzR::openMSfile(mzml_files_pos[i])
	run_info <- mzR::runInfo(fl)
	mzR::close(fl)
	mzml_times_pos <- c(mzml_times_pos, run_info$startTimeStamp)
}

# Display MSn levels
mzml_msn_pos <- NULL
for (i in 1:length(mzml_files_pos)) {
	mzml_data_pos <- readMSData(mzml_files_pos[i], mode="onDisk")
	mzml_msn_pos <- rbind(mzml_msn_pos, t(as.matrix(table(msLevel(mzml_data_pos)))))
}
colnames(mzml_msn_pos) <- c("MS1", "MS2")
rownames(mzml_msn_pos) <- mzml_names_pos

# Plot MSn levels
pdf(file="plots/pos_msn_levels.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=16, family="Helvetica")
par(mfrow=c(2,1), mar=c(16,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.8)
boxplot(mzml_msn_pos, main="Number of spectra")

model_boxplot <- boxplot(t(mzml_msn_pos[,2]), main="Number of MS2 spectra per sample", xaxt="n")
tick <- seq_along(model_boxplot$names)
axis(1, at=tick, labels=F)
text(tick, par("usr")[3]-par("usr")[3]/10, model_boxplot$names, adj=0, srt=270, xpd=T)
dev.off()



# ---------- Peak detection ----------
# Import raw data as MSnbase object
raw_data_pos <- readMSData(files=mzml_files_pos, pdata=new("NAnnotatedDataFrame", mzml_pheno_pos), mode="onDisk", centroided=TRUE)
table(msLevel(raw_data_pos))
head(fData(raw_data_pos)[, c("isolationWindowTargetMZ", "isolationWindowLowerOffset",
						     "isolationWindowUpperOffset", "msLevel", "retentionTime")])
write.csv(fData(raw_data_pos), file="data/pos_raw_data.csv", row.names=FALSE)

# Restrict data to 1020 seconds (17 minutes)
raw_data_pos <- filterRt(raw_data_pos, c(min.rt, max.rt))

# Inspect mz values per file
raw_mz_pos <- mz(raw_data_pos)
raw_mz_pos <- split(raw_mz_pos, f = fromFile(raw_data_pos))
print(length(raw_mz_pos))

# Plot base peak chromatograms
chromas_bpc_pos <- chromatogram(raw_data_pos, aggregationFun="max")

pdf(file="plots/pos_chromas_bpc.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
plot(chromas_bpc_pos, main="Raw base peak chromatograms", xlab="retention time [s]", ylab="intensity", col=mzml_pheno_colors_samples)
legend("topleft", bty="n", pt.cex=0.5, cex=0.7, y.intersp=0.7, text.width=0.5, pch=20, col=mzml_pheno_colors_pos, legend=levels(mzml_pheno_samples))
dev.off()

# Plot TICs
chromas_tic_pos <- chromatogram(raw_data_pos, aggregationFun="sum")

pdf(file="plots/pos_chromas_tic.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
plot(chromas_tic_pos, main="Raw TIC chromatograms", xlab="retention time [s]", ylab="intensity", col=mzml_pheno_colors_samples_pos)
legend("topleft", bty="n", pt.cex=0.5, cex=0.7, y.intersp=0.7, text.width=0.5, pch=20, col=mzml_pheno_colors_pos, legend=levels(mzml_pheno_samples_pos))
dev.off()

# Get TICs
pdf(file="plots/pos_tics.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
par(mfrow=c(1,1), mar=c(19,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
tics_pos <- split(tic(raw_data_pos), f=fromFile(raw_data_pos))
boxplot(tics_pos, names=mzml_names_pos, las=2, col=mzml_pheno_colors_samples_pos, ylab="intensity", main="Total ion current")
dev.off()

# Grouping/binning the samples based on similarity of their base peak chromatogram to spot potentially problematic samples
chromas_bin_bpc_pos <- MSnbase::bin(chromas_bpc_pos, binSize=2)
chromas_bin_bpc_cor_pos <- cor(do.call(cbind, lapply(chromas_bin_bpc_pos, intensity)))
colnames(chromas_bin_bpc_cor_pos) <- rownames(chromas_bin_bpc_cor_pos) <- raw_data_pos$sample_name
pdf(file="plots/pos_chromas_bpc_bin_cor.pdf", encoding="ISOLatin1", pointsize=10, width=8, height=8, family="Helvetica")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
heatmap(chromas_bin_bpc_cor_pos, margins=c(15,15))
dev.off()

# Grouping/binning the samples based on similarity of their TIC to spot potentially problematic samples
chromas_bin_tic_pos <- MSnbase::bin(chromas_tic_pos, binSize=2)
chromas_bin_tic_cor_pos <- cor(do.call(cbind, lapply(chromas_bin_tic_pos, intensity)))
colnames(chromas_bin_tic_cor_pos) <- rownames(chromas_bin_tic_cor_pos) <- raw_data_pos$sample_name
pdf(file="plots/pos_chromas_tic_bin_cor.pdf", encoding="ISOLatin1", pointsize=10, width=8, height=8, family="Helvetica")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
heatmap(chromas_bin_tic_cor_pos, margins=c(15,15))
dev.off()

# Assess retention times and intensities of first file
head(rtime(chromas_pos[1, 1]))
head(intensity(chromas_pos[1, 1]))

# Peak detection in MS1 data
if (polarity=="negative") {
	ms1_params_neg <- CentWaveParam(ppm=25, mzCenterFun="mean", peakwidth=c(7.2, 36), prefilter=c(2, 50), mzdiff=0.0012, snthresh=7, noise=0, integrate=1,
									firstBaselineCheck=TRUE, verboseColumns=TRUE, fitgauss=FALSE, roiList=list(), roiScales=numeric())
} else {
	ms1_params_pos <- CentWaveParam(ppm=25, mzCenterFun="mean", peakwidth=c(9.4, 32), prefilter=c(6, 51), mzdiff=-0.0043, snthresh=2, noise=0, integrate=1,
									firstBaselineCheck=TRUE, verboseColumns=FALSE, fitgauss=FALSE, roiList=list(), roiScales=numeric())
}
ms1_data_pos <- findChromPeaks(raw_data_pos, param=ms1_params_pos)

# Per file summary
ms1_summary_pos <- lapply(split.data.frame(chromPeaks(ms1_data_pos), f=chromPeaks(ms1_data_pos)[, "sample"]), FUN=function(z) { c(peak_count=nrow(z), rt=quantile(z[, "rtmax"] - z[, "rtmin"])) } )
ms1_summary_pos <- do.call(rbind, ms1_summary_pos)
rownames(ms1_summary_pos) <- basename(fileNames(ms1_data_pos))
print(ms1_summary_pos)
table(msLevel(ms1_data_pos))
write.csv(as.data.frame(table(msLevel(ms1_data_pos))), file="data/pos_ms1_data.csv", row.names=FALSE)

# To get a global overview of the peak detection we can plot the frequency of identified peaks per file along the retention time axis. This allows to identify time periods along the MS run in which a higher number of peaks was identified and evaluate whether this is consistent across files.
pdf(file="plots/pos_ms1_data.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
par(mfrow=c(1,1), mar=c(4,18,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
plotChromPeakImage(ms1_data_pos, main="Frequency of identified peaks per RT")
dev.off()

# Group peaks
if (polarity=="negative") {
	ms1_data_neg <- groupChromPeaks(ms1_data_neg, param=PeakDensityParam(sampleGroups=ms1_data_neg$sample_group, minFraction=0.7, bw=0.25, minSamples=1, binSize=0.5))
} else {
	ms1_data_pos <- groupChromPeaks(ms1_data_pos, param=PeakDensityParam(sampleGroups=ms1_data_pos$sample_group, minFraction=0.7, bw=0.25, minSamples=1, binSize=0.5))
}

# RT correction
if (polarity=="negative") {
	ms1_data_neg <- adjustRtime(ms1_data_neg, param=PeakGroupsParam(minFraction=0.7, smooth="loess", span=0.5, extra=1, family="gaussian"))
} else {
	ms1_data_pos <- adjustRtime(ms1_data_pos, param=PeakGroupsParam(minFraction=0.7, smooth="loess", span=0.5, extra=1, family="gaussian"))
}

# Plot the difference of raw and adjusted retention times
pdf(file="plots/pos_ms1_raw_adjusted.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=8, family="Helvetica")
par(mfrow=c(2,1), mar=c(4.5,4.2,4,1), cex=0.8)
plot(chromas_pos, main="Raw chromatograms", col=mzml_pheno_colors_samples_pos)
legend("topleft", bty="n", pt.cex=0.5, cex=0.7, y.intersp=0.7, text.width=0.5, pch=20, col=mzml_pheno_colors_pos, legend=levels(mzml_pheno_samples_pos))
plotAdjustedRtime(ms1_data_pos, lwd=2, ylim=c(-2.2,2.2), main="Retention Time correction", col=mzml_pheno_colors_samples_pos)
legend("topleft", bty="n", pt.cex=0.5, cex=0.7, y.intersp=0.7, text.width=0.5, pch=20, col=mzml_pheno_colors_pos, legend=levels(mzml_pheno_samples_pos))
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.8)
dev.off()

# Group peaks
if (polarity=="negative") {
	ms1_data_neg <- groupChromPeaks(ms1_data_neg, param=PeakDensityParam(sampleGroups=ms1_data_neg$sample_group, minFraction=0.7, bw=0.25, minSamples=1, binSize=0.5))
} else {
	ms1_data_pos <- groupChromPeaks(ms1_data_pos, param=PeakDensityParam(sampleGroups=ms1_data_pos$sample_group, minFraction=0.7, bw=0.25, minSamples=1, binSize=0.5))
}

# Get integrated peak intensity per feature/sample
print(head(featureValues(ms1_data_pos, value="into")))

# Fill peaks
#ms1_data_pos <- fillChromPeaks(ms1_data_pos, param=FillChromPeaksParam(ppm=ppm, fixedRt=0, expandRt=5))
#head(featureValues(ms1_data_pos))
#head(featureSummary(ms1_data_pos, group=ms1_data_pos$sample_group))

# Evaluate grouping
pdf(file="plots/pos_ms1_grouping.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
ms1_pca_pos <- prcomp(t(na.omit(log2(featureValues(ms1_data_pos, value="into")))), center=TRUE)
plot(ms1_pca_pos$x[, 1], ms1_pca_pos$x[, 2], pch=19, main="PCA: Grouping of samples",
	 xlab=paste0("PC1: ", format(summary(ms1_pca_pos)$importance[2, 1] * 100, digits=3), " % variance"),
	 ylab=paste0("PC2: ", format(summary(ms1_pca_pos)$importance[2, 2] * 100, digits=3), " % variance"),
	 col=mzml_pheno_colors_samples_pos, cex=2)
grid()
text(ms1_pca_pos$x[, 1], ms1_pca_pos$x[, 2], labels=ms1_data_pos$sample_name, col=mzml_pheno_colors_samples_pos, pos=3, cex=0.5)
dev.off()

# Evaluate QC Deviations
pdf(file="plots/pos_ms1_mzrt_deviations.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=8, family="Helvetica")
par(mfrow=c(2,1), mar=c(14.5,4,4,1), cex=0.8, cex.axis=0.6)
plotQC(ms1_data_pos, what="mzdevsample", sampColors=mzml_pheno_colors_samples_pos, sampNames=mzml_pheno_pos$sample_name)
title(main="Median mz deviation for each sample")
plotQC(ms1_data_pos, what="rtdevsample", sampColors=mzml_pheno_colors_samples_pos, sampNames=mzml_pheno_pos$sample_name)
title(main="Median RT deviation for each sample")
dev.off()

# Show peaks
tail(chromPeaks(ms1_data_pos))
tail(chromPeakData(ms1_data_pos))

# Show process history
processHistory(ms1_data_pos)



# ---------- Build MS1 feature tables ----------
# Build feature matrix
ms1_matrix_pos <- featureValues(ms1_data_pos, method="medret", value="into")
colnames(ms1_matrix_pos) <- mzml_names_pos
dim(ms1_matrix_pos)
feat_list_pos <- t(ms1_matrix_pos)

# Build feature summary
ms1_summary_pos <- featureSummary(ms1_data_pos)
ms1_def_pos <- featureDefinitions(ms1_data_pos)

# Transform data
feat_list_pos <- log2(feat_list_pos)

# Missing value imputation
#feat_list_pos[which(is.na(feat_list_pos))] <- median(na.omit(as.numeric(unlist(feat_list_pos))))
feat_list_pos <- scale(feat_list_pos, scale=FALSE, center=FALSE)
feat_list_pos[is.na(feat_list_pos)] <- 0
feat_list_pos[which(feat_list_pos < 0)] <- 0
feat_list_pos[is.infinite(feat_list_pos)] <- 0
feat_list_pos <- feat_list_pos[!apply(feat_list_pos, MARGIN=1, function(x) max(x,na.rm=TRUE) == min(x,na.rm=TRUE)),]

# Plot histogram
pdf(file="plots/pos_feat_list_hist.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
hist(as.numeric(feat_list_pos), main="Histogram of feature table")
dev.off()

# PCA
pdf(file="plots/pos_ms1_feature_table_pca.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
ms1_pca_pos <- prcomp(feat_list_pos, center=TRUE)
plot(ms1_pca_pos$x[, 1], ms1_pca_pos$x[, 2], pch=19, main="PCA of feature table",
	 xlab=paste0("PC1: ", format(summary(ms1_pca_pos)$importance[2, 1] * 100, digits=3), " % variance"),
	 ylab=paste0("PC2: ", format(summary(ms1_pca_pos)$importance[2, 2] * 100, digits=3), " % variance"),
	 col=mzml_pheno_colors_samples_pos, cex=2)
grid()
text(ms1_pca_pos$x[, 1], ms1_pca_pos$x[, 2], labels=ms1_data_pos$sample_name, col=mzml_pheno_colors_samples_pos, pos=3, cex=0.5)

plot(ms1_pca_pos$x[, 3], ms1_pca_pos$x[, 4], pch=19, main="PCA of feature table",
	 xlab=paste0("PC3: ", format(summary(ms1_pca_pos)$importance[2, 3] * 100, digits=3), " % variance"),
	 ylab=paste0("PC4: ", format(summary(ms1_pca_pos)$importance[2, 4] * 100, digits=3), " % variance"),
	 col=mzml_pheno_colors_samples_pos, cex=2)
grid()
text(ms1_pca_pos$x[, 3], ms1_pca_pos$x[, 4], labels=ms1_data_pos$sample_name, col=mzml_pheno_colors_samples_pos, pos=3, cex=0.5)

plot(ms1_pca_pos$x[, 5], ms1_pca_pos$x[, 6], pch=19, main="PCA of feature table",
	 xlab=paste0("PC5: ", format(summary(ms1_pca_pos)$importance[2, 5] * 100, digits=3), " % variance"),
	 ylab=paste0("PC6: ", format(summary(ms1_pca_pos)$importance[2, 6] * 100, digits=3), " % variance"),
	 col=mzml_pheno_colors_samples_pos, cex=2)
grid()
text(ms1_pca_pos$x[, 5], ms1_pca_pos$x[, 6], labels=ms1_data_pos$sample_name, col=mzml_pheno_colors_samples_pos, pos=3, cex=0.5)

plot(ms1_pca_pos$x[, 7], ms1_pca_pos$x[, 8], pch=19, main="PCA of feature table",
	 xlab=paste0("PC7: ", format(summary(ms1_pca_pos)$importance[2, 7] * 100, digits=3), " % variance"),
	 ylab=paste0("PC8: ", format(summary(ms1_pca_pos)$importance[2, 8] * 100, digits=3), " % variance"),
	 col=mzml_pheno_colors_samples_pos, cex=2)
grid()
text(ms1_pca_pos$x[, 7], ms1_pca_pos$x[, 8], labels=ms1_data_pos$sample_name, col=mzml_pheno_colors_samples_pos, pos=3, cex=0.5)

plot(ms1_pca_pos$x[, 9], ms1_pca_pos$x[, 10], pch=19, main="PCA of feature table",
	 xlab=paste0("PC9: ", format(summary(ms1_pca_pos)$importance[2, 9] * 100, digits=3), " % variance"),
	 ylab=paste0("PC10: ", format(summary(ms1_pca_pos)$importance[2, 10] * 100, digits=3), " % variance"),
	 col=mzml_pheno_colors_samples_pos, cex=2)
grid()
text(ms1_pca_pos$x[, 9], ms1_pca_pos$x[, 10], labels=ms1_data_pos$sample_name, col=mzml_pheno_colors_samples_pos, pos=3, cex=0.5)

screeplot(ms1_pca_pos, bstick=TRUE, type=c("barplot", "lines"), npcs=min(20, length(ms1_pca_pos$sdev)), ptype="o", bst.col="red", bst.lty="solid", xlab="Component", ylab="Inertia", main="Broken stick test of PCA", legend=FALSE)
dev.off()

# Create single 0/1 matrix
bina_list_pos <- t(ms1_matrix_pos)
bina_list_pos[is.na(bina_list_pos)] <- 1
bina_list_pos <- log2(bina_list_pos)
bina_list_pos[bina_list_pos < log2(ms1_intensity_cutoff)] <- 0
bina_list_pos[bina_list_pos != 0] <- 1

# Only unique compounds in group mzml_pheno$ and not the others
uniq_list_pos <- apply(X=bina_list_pos, MARGIN=2, FUN=function(x) { if (length(unique(mzml_pheno_pos$sample_group[grepl("1", x)])) == 1) x else rep(0, length(x)) } )
colnames(uniq_list_pos) <- colnames(bina_list_pos)
rownames(uniq_list_pos) <- rownames(bina_list_pos)

# Create data frame
model_div_pos             <- data.frame(features=apply(X=bina_list_pos, MARGIN=1, FUN=function(x) { sum(x) } ))
model_div_pos$richness    <- apply(X=bina_list_pos, MARGIN=1, FUN=function(x) { sum(x) } )
#model_div_pos$menhinick   <- apply(X=bina_list_pos, MARGIN=1, FUN=function(x) { menhinick.diversity(x) } )
model_div_pos$shannon     <- apply(X=feat_list_pos, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") })
model_div_pos$pielou      <- apply(X=scale(feat_list_pos, center=F), MARGIN=1, FUN=function(x) { vegan::diversity(x, index="shannon") / log(vegan::specnumber(x)) })
#model_div_pos$chao        <- vegan::specpool2vect(X=vegan::specpool(feat_list_pos, species), index="chao")
model_div_pos$simpson     <- apply(X=feat_list_pos, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="simpson") })
model_div_pos$inverse     <- apply(X=feat_list_pos, MARGIN=1, FUN=function(x) { vegan::diversity(x, index="inv") })
model_div_pos$fisher      <- apply(X=feat_list_pos, MARGIN=1, FUN=function(x) { fisher.alpha(round(x,0)) })
model_div_pos$unique      <- apply(X=uniq_list_pos, MARGIN=1, FUN=function(x) { sum(x) })
model_div_pos$hillfunc    <- as.numeric(unlist(calcDiv(feat_list_pos, compDisMat=scales::rescale(as.matrix(dist(t(feat_list_pos)), diag=TRUE, upper=TRUE)), q=1, type="FuncHillDiv")))

# Remove NAs if present
model_div_pos[is.na(model_div_pos)] <- 0



# ---------- QC compounds ----------
# Check whether QC compounds have been picked by XCMS
QC_pos <- NULL
QC_pos <- rbind(QC_pos, data.frame(compound="N-(3-Indolylacetyl)-L-Valine", rt=462, mz=275.1343))
QC_pos <- rbind(QC_pos, data.frame(compound="Biochanin A", rt=585, mz=285.0711))
QC_pos <- rbind(QC_pos, data.frame(compound="Kinetin", rt=281, mz=216.086))
QC_pos$FT <- ""
QC_pos$median_into <- 0
QC_pos$sd_into <- 0

for (i in 1:nrow(QC_pos)) {
	QC_pos$FT[i] <- rownames(ms1_def_pos[ms1_def_pos$rtmed>=QC_pos$rt[i]-rtabs & ms1_def_pos$rtmed<=QC_pos$rt[i]+rtabs & ms1_def_pos$mzmed>=QC_pos$mz[i]-mzabs & ms1_def_pos$mzmed<=QC_pos$mz[i]+mzabs, ])[1]
	QC_pos$median_into[i] <- median(ms1_matrix_pos[which(rownames(ms1_matrix_pos)==QC_pos$FT[i]), ], na.rm=TRUE)
	QC_pos$sd_into[i] <- sd(ms1_matrix_pos[which(rownames(ms1_matrix_pos)==QC_pos$FT[i]), ], na.rm=TRUE)
}

# Plot variation of MM8 compounds in the samples
pdf(file="plots/pos_ms1_qc_compounds_into.pdf", encoding="ISOLatin1", pointsize=12, width=5, height=5, family="Helvetica")
par(mar=c(7,4,1,1))
boxplot(t(ms1_matrix_pos[rownames(ms1_matrix_pos) %in% QC_pos$FT,]), main="", xlab="", ylab="log(into)", names=NA)
text(1:nrow(QC_pos), par("usr")[3]-(par("usr")[4]-par("usr")[3])/24, srt=-40, adj=c(0,1), labels=QC_pos$compound, xpd=TRUE, cex=0.8)
dev.off()

# Plot extracted ion chromatograms
chromas_eic_pos <- chromatogram(ms1_data_pos, aggregationFun="sum", adjustedRtime=TRUE, include="any", msLevel=1L)

# Plot chromatograms of QC compounds
QC_pos[QC_pos$compound=="Biochanin A", "rt"] <- 590
pdf(file="plots/pos_ms1_qc_compounds_chromatogram.pdf", encoding="ISOLatin1", pointsize=10, width=3, height=9, family="Helvetica")
par(mfrow=c(3,1), mar=c(4,4,3,1), cex=0.8)
for (j in 1:nrow(QC_pos)) {
	ymax = NULL
	for (i in 1:ncol(ms1_matrix_pos)) ymax = c(ymax, chromas_eic_pos@.Data[[i]]@intensity[chromas_eic_pos@.Data[[i]]@rtime>=QC_pos$rt[j]-rtabs & chromas_eic_pos@.Data[[i]]@rtime<=QC_pos$rt[j]+rtabs])
	#for (i in 1:ncol(ms1_matrix_pos)) ymax = c(ymax, chromas_tic_pos@.Data[[i]]@intensity[chromas_tic_pos@.Data[[i]]@rtime>=QC_pos$rt[j]-rtabs & chromas_tic_pos@.Data[[i]]@rtime<=QC_pos$rt[j]+rtabs])
	#for (i in 1:ncol(ms1_matrix_pos)) ymax = c(ymax, chromas_bpc_pos@.Data[[i]]@intensity[chromas_bpc_pos@.Data[[i]]@rtime>=QC_pos$rt[j]-rtabs & chromas_bpc_pos@.Data[[i]]@rtime<=QC_pos$rt[j]+rtabs])
	ymin = floor(min(ymax))
	ymax = ceiling(max(ymax))
	plot(0, 0, xlim=c(QC_pos[j,"rt"]-rtabs,QC_pos[j,"rt"]+rtabs), ylim=c(ymin,ymax), main=QC_pos[j,"compound"], xlab="RT [s]", ylab="log(into)")
	for (i in 1:ncol(ms1_matrix_pos)) {
		lines(chromas_eic_pos@.Data[[i]]@rtime, chromas_eic_pos@.Data[[i]]@intensity, lwd=2, col=mzml_pheno_colors_samples_pos[i])
		#lines(chromas_tic_pos@.Data[[i]]@rtime, chromas_tic_pos@.Data[[i]]@intensity, lwd=2, col=mzml_pheno_colors_samples_pos[i])
		#lines(chromas_bpc_pos@.Data[[i]]@rtime, chromas_bpc_pos@.Data[[i]]@intensity, lwd=2, col=mzml_pheno_colors_samples_pos[i])
	}
}
dev.off()



# ############################## MS2 ##############################



# ---------- MS2 spectra detection ----------
# Estimate precursor intensity
precursor_intensity_pos <- xcms::estimatePrecursorIntensity(ms1_data_pos)
print(head(na.omit(precursor_intensity_pos)))

# Reconstruct MS2 spectra from MS1 data
ms2_data_pos <- chromPeakSpectra(ms1_data_pos, msLevel=2L, return.type="Spectra")
print(ms2_data_pos)
print(length(ms2_data_pos$peak_id))

# Extract all usable MS2 spectra
ms2_spectra_pos <- list()
for (i in 1:nrow(ms1_def_pos)) {
	#ms2_spectra_pos <- foreach(i=1:nrow(ms1_def_pos)) %dopar% {
	# Extract existing MS2 spectra for feature
	feature_of_interest <- ms1_def_pos[i, "mzmed"]
	peaks_of_interest <- chromPeaks(ms1_data_pos, mz=feature_of_interest, ppm=ppm)
	spectra_of_interest <- ms2_data_pos[ms2_data_pos$peak_id %in% rownames(peaks_of_interest)]
	
	# Continue if feature has MS2 peaks
	if (length(which(ms2_data_pos$peak_id %in% rownames(peaks_of_interest))) > 0) {
		# Filter for retention time (based on parameter rtabs)
		#combined_spectra_of_interest <- filterIntensity(spectra_of_interest, intensity=c(1,Inf), backend=MsBackendDataFrame)
		combined_spectra_of_interest <- Spectra::filterRt(spectra_of_interest, rt=c(ms1_def_pos[i, "rtmin"]-(rtabs/2),ms1_def_pos[i, "rtmax"]+(rtabs/2)))
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
	
	ms2_spectra_pos[[i]] <- combined_spectra_of_interest
	#return(combined_spectra_of_interest)
}

# Remove empty spectra
names(ms2_spectra_pos) <- rownames(ms1_def_pos)
ms2_spectra_pos <- ms2_spectra_pos[lengths(ms2_spectra_pos) != 0]

# Relate all MS2 spectra to MS1 precursors
ms1_def_pos$has_ms2 <- as.integer(rownames(ms1_def_pos) %in% names(ms2_spectra_pos))
print(paste0("Number of MS2 spectra related to precursor: ", length(which(ms1_def_pos$has_ms2>0))))



# ---------- Annotate MS2 spectra ----------
# Save all MS2 spectra in MSP file
msp_text <- NULL
alignment_id <- 0
for (i in names(ms2_spectra_pos)) {
	alignment_id <- alignment_id + 1
	
	msp_text <- c(msp_text, paste("NAME:", i))
	msp_text <- c(msp_text, paste("AlignmentID:", alignment_id))
	msp_text <- c(msp_text, paste("RETENTIONTIME:", ms1_def_pos[i, "rtmed"]))
	msp_text <- c(msp_text, paste("PRECURSORMZ:", ms1_def_pos[i, "mzmed"]))
	msp_text <- c(msp_text, paste("METABOLITENAME:", i))
	if (polarity == "positive") {
		msp_text <- c(msp_text, paste("ADDUCTIONNAME:", "[M+H]+"))
	} else {
		msp_text <- c(msp_text, paste("ADDUCTIONNAME:", "[M-H]-"))
	}
	msp_text <- c(msp_text, paste("NumPeaks:", nrow(as.data.frame(peaksData(ms2_spectra_pos[[i]])[[1]]))))
	msp_text <- c(msp_text, paste(as.data.frame(peaksData(ms2_spectra_pos[[i]])[[1]])$mz, as.data.frame(peaksData(ms2_spectra_pos[[i]])[[1]])$intensity, sep="\t"))
	msp_text <- c(msp_text, "")
	
}

# Write MSP file
cat(msp_text, file="data/pos_ms2_spectra.msp", sep="\n")


# Save all MS2 spectra in MGF file
mgf_text <- NULL
for (i in names(ms2_spectra_pos)) {
	mgf_text <- c(mgf_text, paste0("COM=", i))
	mgf_text <- c(mgf_text, "BEGIN IONS")
	mgf_text <- c(mgf_text, "MSLEVEL=2")
	mgf_text <- c(mgf_text, paste0("TITLE=", i))
	mgf_text <- c(mgf_text, paste0("RTINSECONDS=", ms1_def_pos[i, "rtmed"]))
	mgf_text <- c(mgf_text, paste0("PEPMASS=", ms1_def_pos[i, "mzmed"]))
	if (polarity == "positive") {
		mgf_text <- c(mgf_text, paste0("CHARGE=", "1+"))
	} else {
		mgf_text <- c(mgf_text, paste0("CHARGE=", "1-"))
	}
	mgf_text <- c(mgf_text, paste(as.data.frame(peaksData(ms2_spectra_pos[[i]])[[1]])$mz, as.data.frame(peaksData(ms2_spectra_pos[[i]])[[1]])$intensity, sep=" "))
	mgf_text <- c(mgf_text, "END IONS")
	mgf_text <- c(mgf_text, "")
}

# Write MGF file
cat(mgf_text, file="data/pos_ms2_spectra.mgf", sep="\n")



# ---------- Annotate MS2 spectra with SIRIUS ----------
# Apply annotated compounds onto feature table
ms1_def_pos$has_id <- 0
ms1_def_pos$smiles <- ""
ms1_def_pos$inchi <- ""
ms1_def_pos$inchikey <- ""
ms1_def_pos$name <- ""
ms1_def_pos$classyfire <- ""

# Read SIRIUS annotation
annotator_sirius_pos <- read.table(file=paste0("data/pos_ms2_compound_identifications.tsv"), header=TRUE, sep="\t", quote="\"", comment.char="", fill=TRUE, dec=".", stringsAsFactors=FALSE)#, encoding="UTF-8", fileEncoding="UTF-8")
annotator_sirius_pos$Metabolite.name <- gsub(x=annotator_sirius_pos$id, pattern=".*_", replacement="")

for (j in unique(annotator_sirius_pos$Metabolite.name)) {
	ms1_def_pos[j, "has_id"] <- 1
	
	obj <- annotator_sirius_pos[annotator_sirius_pos$Metabolite.name %in% j, "smiles"]
	ms1_def_pos[j, "smiles"] <- obj[1]
	
	obj <- annotator_sirius_pos[annotator_sirius_pos$Metabolite.name %in% j, "InChI"]
	ms1_def_pos[j, "inchi"] <- obj[1]
	
	obj <- annotator_sirius_pos[annotator_sirius_pos$Metabolite.name %in% j, "name"]
	if (obj[1] != "null" ) {
		ms1_def_pos[j, "name"] <- obj[1]
	}
}

# Convert InChi to InChiKey
for (i in 1:nrow(ms1_def_pos)) {
	if ((ms1_def_pos[i, "smiles"] != "") & (ms1_def_pos[i, "inchikey"] == "")) {
		ms1_def_pos[i, "inchikey"] <- rinchi::get.inchi.key(ms1_def_pos[i, "smiles"])
	}
}

# ClassyFire of InChiKey
for (i in 1:nrow(ms1_def_pos)) {
	if ((ms1_def_pos[i, "inchikey"] != "") & (ms1_def_pos[i, "classyfire"] =="")) {
		ms1_def_pos[i, "classyfire"] <- paste(as.data.frame(get_classification(ms1_def_pos[i, "inchikey"])@classification)$Classification, collapse="; ")
	}
}

# Annotation
sirius_compounds_pos <- annotator_sirius_pos$smiles
sirius_compound_ids_pos <- gsub(x=annotator_sirius_pos$id, pattern=".*_", replacement="")

# Calculate molecular descriptors with CDK
cdk_descriptors_pos <- NULL
for (i in sirius_compounds_pos) {
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
	cdk_descriptors_pos <- plyr::rbind.fill(cdk_descriptors_pos, cbind(data.frame(MolWeight=MolWeight, nAtoms=nAtoms, numC=numC, numN=numN, numP=numP, numO=numO, numHydrogen=numHydrogen, CNRatio=CNRatio, XLogP=XLogP), cdk_mol_des))
}
rownames(cdk_descriptors_pos) <- paste0(sirius_compound_ids_pos,"_",pol)

# Properly format NAs and convert to numeric
for (i in 1:nrow(cdk_descriptors_pos)) {
	cdk_descriptors_pos[i, as.character(cdk_descriptors_pos[i,]) %in% 'list(NULL)'] <- 0
	cdk_descriptors_pos[i, as.character(cdk_descriptors_pos[i,]) %in% 'NA'] <- 0
	cdk_descriptors_pos[i, as.character(cdk_descriptors_pos[i,]) %in% 'NULL'] <- 0
	cdk_descriptors_pos[i, as.character(cdk_descriptors_pos[i,]) %in% 'NaN'] <- 0
	cdk_descriptors_pos[i, ] <- as.numeric(unlist(cdk_descriptors_pos[i,]))
}

# Bug: Remove left-over variables
rm(CNRatio, MolWeight, XLogP, nAtoms, numC, numHydrogen, numN, numO, numP)



# ---------- Classify MS2 spectra with CANOPUS ----------
# Apply classified classes onto feature table
ms1_def_pos$primary_class <- ""
ms1_def_pos$alternative_classes <- ""
ms1_def_pos$NP_class <- ""

# Read SIRIUS/CANOPUS classifier
if (SIRIUS_VERSION == 4) {
	classifier_canopus_pos <- read.table(file=paste0("data/pos_ms2_canopus_compound_summary.tsv"), header=TRUE, sep="\t", quote="\"", fill=FALSE, dec=".", stringsAsFactors=FALSE)#, encoding="UTF-8", fileEncoding="UTF-8")
	classifier_canopus_pos$Metabolite.name <- gsub(x=classifier_canopus_pos$name, pattern=".*_", replacement="")
	classifier_canopus_pos$primary_class <- paste("Organic compounds", classifier_canopus_pos$superclass, classifier_canopus_pos$class, classifier_canopus_pos$subclass, classifier_canopus_pos$level.5, sep="; ")
	classifier_canopus_pos$primary_class <- gsub(x=classifier_canopus_pos$primary_class, pattern="; ; ", replacement="; ", perl=TRUE)
	classifier_canopus_pos$primary_class <- gsub(x=classifier_canopus_pos$primary_class, pattern="; ; ", replacement="; ", perl=TRUE)
	classifier_canopus_pos$primary_class <- gsub(x=classifier_canopus_pos$primary_class, pattern="; $", replacement="", perl=TRUE)
	classifier_canopus_pos$primary_class <- gsub(x=classifier_canopus_pos$primary_class, pattern="(\\'|\\>|\\(|\\))", replacement="", perl=TRUE)
	classifier_canopus_pos$Annotation..putative. <- classifier_canopus_pos$primary_class
	classifier_canopus_pos$alternative_classes <- classifier_canopus_pos$all.classifications
} else {
	classifier_canopus_pos <- read.table(file=paste0("data/pos_ms2_canopus_compound_summary.tsv"), header=TRUE, sep="\t", quote="\"", comment.char="", fill=TRUE, dec=".", stringsAsFactors=FALSE)#, encoding="UTF-8", fileEncoding="UTF-8")
	classifier_canopus_pos$Metabolite.name <- ""
	for (i in 1:length(classifier_canopus_pos$id)) {
		x = which(gsub(x=classifier_canopus_pos$id, pattern=".*?_", replacement="", perl=TRUE) %in% gsub(x=classifier_canopus_pos$id[i], pattern=".*?_", replacement="", perl=TRUE))
		if (length(x) > 0) classifier_canopus_pos$Metabolite.name[i] <- gsub(x=classifier_canopus_pos$id[x], pattern=".*?_", replacement="", perl=TRUE)
	}
	classifier_canopus_pos$Metabolite.name[classifier_canopus_pos$Metabolite.name == "null"] <- ""
	classifier_canopus_pos$primary_class <- paste("Organic compounds", classifier_canopus_pos$ClassyFire.superclass, classifier_canopus_pos$ClassyFire.class, classifier_canopus_pos$ClassyFire.subclass, classifier_canopus_pos$ClassyFire.level.5, sep="; ")
	classifier_canopus_pos$primary_class <- gsub(x=classifier_canopus_pos$primary_class, pattern="; ; ", replacement="; ", perl=TRUE)
	classifier_canopus_pos$primary_class <- gsub(x=classifier_canopus_pos$primary_class, pattern="; ; ", replacement="; ", perl=TRUE)
	classifier_canopus_pos$primary_class <- gsub(x=classifier_canopus_pos$primary_class, pattern="; $", replacement="", perl=TRUE)
	classifier_canopus_pos$primary_class <- gsub(x=classifier_canopus_pos$primary_class, pattern="(\\'|\\>|\\(|\\))", replacement="", perl=TRUE)
	classifier_canopus_pos$Annotation..putative. <- classifier_canopus_pos$primary_class
	classifier_canopus_pos$alternative_classes <- classifier_canopus_pos$all.classifications
}

for (j in unique(classifier_canopus_pos$Metabolite.name)) {
	ms1_def_pos[j, "has_id"] <- 1
	
	# CHEMONT
	obj <- classifier_canopus_pos[classifier_canopus_pos$Metabolite.name %in% j, "Annotation..putative."]
	ms1_def_pos[j, "primary_class"] <- obj[1]
	
	# NPClassifier
	if (SIRIUS_VERSION > 4) {
		obj <- paste(sep="; ", classifier_canopus_pos[classifier_canopus_pos$Metabolite.name %in% j, "NPC.pathway"],
					           classifier_canopus_pos[classifier_canopus_pos$Metabolite.name %in% j, "NPC.superclass"],
					           classifier_canopus_pos[classifier_canopus_pos$Metabolite.name %in% j, "NPC.class"])
		ms1_def_pos[j, "NP_class"] <- obj[1]
	}
}



# ---------- Diversity of MS2 classes ----------
# Create CANOPUS classifier object for each sample
classifiers_pos <- list()
for (i in mzml_names_pos) {
	obj <- names(which(bina_list_pos[rownames(bina_list_pos)==i, colnames(bina_list_pos) %in% classifier_canopus_pos$Metabolite.name] > 0))
	classifiers_pos[[i]] <- classifier_canopus_pos[classifier_canopus_pos$Metabolite.name %in% obj, ]
}

# Diversity of classes per sample
div_classes_samples_pos <- NULL
for (i in mzml_names_pos) {
	obj <- table(classifiers_pos[[i]][,"Annotation..putative."])
	obj <- data.frame(classes=names(obj), frequency=as.numeric(obj))
	if (is.null(div_classes_samples_pos)) {
		div_classes_samples_pos <- obj
	} else {
		div_classes_samples_pos <- merge(div_classes_samples_pos, obj, by="classes", all.x=TRUE, all.y=TRUE)
	}
}
rownames(div_classes_samples_pos) <- div_classes_samples_pos$classes
div_classes_samples_pos <- div_classes_samples_pos[, -which(colnames(div_classes_samples_pos)=="classes")]
colnames(div_classes_samples_pos) <- mzml_names_pos

# Diversity of classes
div_classes_pos <- div_classes_samples_pos
div_classes_pos[is.na(div_classes_pos)] <- 0
div_classes_pos <- apply(X=div_classes_pos, MARGIN=1, FUN=function(x) { sum(x) })
div_classes_pos <- data.frame(row.names=names(div_classes_pos), frequency=as.numeric(div_classes_pos))

# Plot diversity of classes in all samples
pdf(file="plots/pos_ms2_classes_diversity.pdf", encoding="ISOLatin1", pointsize=8, width=6, height=14, family="Helvetica")
par(mfrow=c(1,1), mar=c(4,15,4,1), oma=c(0,0,0,0), cex.axis=0.8, cex=0.9)
barplot(div_classes_pos$frequency, names.arg=gsub('.*; ','',rownames(div_classes_pos)), las=1, horiz=TRUE, xlab="frequency", main="Diversity of compound classes", col=rainbow(n=nrow(div_classes_pos), alpha=0.6))
dev.off()

# Sunburst plot of classes of all samples
pdf(file="plots/pos_ms2_classes_sunburst.pdf", encoding="ISOLatin1", pointsize=8, width=10, height=10, family="Helvetica")
par(mfrow=c(1,1), mar=c(0,0,0,0), oma=c(0,0,0,0), cex.axis=1, cex=1)
sunBurstPlotFromSubstanceClasses(rownames(div_classes_pos), div_classes_pos$frequency, colorStart=0.0, colorAlpha=0.6)
dev.off()
write.csv(data.frame("Compound classes"=rownames(div_classes_pos), "Number of spectra"=div_classes_pos$frequency), file="plots/pos_ms2_classes_sunburst.csv", row.names=FALSE)

# Classes
classes_pos <- rownames(div_classes_pos)
classes_pos <- classes_pos[which(grepl(pattern="^Organic compounds", x=classes_pos))]
classes_pos <- gsub(x=classes_pos, pattern=":", replacement="")
classes_pos <- gsub(x=classes_pos, pattern="/", replacement="; ")



# ---------- Build diversity objects ----------
# Imputation of NA with zeros
div_classes_pos[is.na(div_classes_pos)] <- 0
div_classes_samples_pos[is.na(div_classes_samples_pos)] <- 0

# Classification list for statistics
class_list_pos <- as.data.frame(t(div_classes_samples_pos))
class_list_pos[is.na(class_list_pos)] <- 0

# Log Transform
#class_list_pos <- log2(class_list_pos + 1)

# Only keep class names
colnames(class_list_pos) <- gsub(x=colnames(class_list_pos), pattern='.*; ', replacement='')

# Generate class_int_list_pos with abundances instead of counts
class_int_list_pos <- class_list_pos

for (i in 1:nrow(class_list_pos)) {
	samp <- rownames(class_list_pos)[i]
	for (j in 1:ncol(class_list_pos)) {
		cl <- colnames(class_list_pos)[j]
		ft <- classifier_canopus_pos$Metabolite.name[which(gsub(x=classifier_canopus_pos$primary_class, pattern='.*; ', replacement='') == cl)]
		ints <- as.numeric(feat_list_pos[i, which(colnames(feat_list_pos) %in% ft)])
		class_int_list_pos[i, j] <- sum(ints)
	}
}



# ---------- Classification at CHEMONT level of classes ----------
# Make superclasses at CHEMONT level 2
superclass_level_pos <- 2
div_superclasses_samples_names_pos <- NULL
for (i in c(1:superclass_level_pos)) {
	div_superclasses_samples_names_pos <- c(div_superclasses_samples_names_pos, lapply(X=strsplit(rownames(div_classes_samples_pos), '; '), FUN=function(x) { gsub(x=paste(x[1:i],sep='',collapse='; '),pattern='; NA',replacement='') }))
}
div_superclasses_samples_pos <- data.frame()
for (i in c(1:ncol(div_classes_samples_pos))) div_superclasses_samples_pos <- rbind(div_superclasses_samples_pos, rep(0, length(unique(div_superclasses_samples_names_pos))))
div_superclasses_samples_pos <- t(div_superclasses_samples_pos)
colnames(div_superclasses_samples_pos) <- colnames(div_classes_samples_pos)
rownames(div_superclasses_samples_pos) <- unique(div_superclasses_samples_names_pos)
for (i in rownames(div_superclasses_samples_pos)) {
	for (j in c(1:ncol(div_classes_samples_pos))) {
		div_superclasses_samples_pos[rownames(div_superclasses_samples_pos)==i, j] <- sum(div_classes_samples_pos[grep(x=rownames(div_classes_samples_pos), pattern=i), j])
	}
}

# Diversity of superclasses
div_superclasses_pos <- div_superclasses_samples_pos
div_superclasses_pos[is.na(div_superclasses_pos)] <- 0
div_superclasses_pos <- apply(X=div_superclasses_pos, MARGIN=1, FUN=function(x) { sum(x) })
div_superclasses_pos <- data.frame(row.names=names(div_superclasses_pos), frequency=as.numeric(div_superclasses_pos))

# Sunburst plot
pdf(file="plots/pos_ms2_superclasses_sunburst.pdf", encoding="ISOLatin1", pointsize=8, width=10, height=10, family="Helvetica")
sunBurstPlotFromSubstanceClasses(rownames(div_superclasses_pos), div_superclasses_pos$frequency, colorStart=0.0, colorAlpha=0.6)
dev.off()
write.csv(data.frame("Compound classes"=rownames(div_superclasses_pos), "Number of spectra"=div_superclasses_pos$frequency), file="plots/pos_ms2_superclasses_sunburst.csv", row.names=FALSE)

# Imputation of NA with zeros
div_superclasses_pos[is.na(div_superclasses_pos)] <- 0
div_superclasses_samples_pos[is.na(div_superclasses_samples_pos)] <- 0

# Classification list for statistics
superclass_list_pos <- as.data.frame(t(div_superclasses_samples_pos))
superclass_list_pos[is.na(superclass_list_pos)] <- 0

# Log Transform
#superclass_list_pos <- log2(superclass_list_pos + 1)

# Only keep superclass names
colnames(superclass_list_pos) <- gsub(x=colnames(superclass_list_pos), pattern='.*; ', replacement='')

# Generate superclass_int_list_pos with abundances instead of counts
superclass_int_list_pos <- superclass_list_pos

for (i in 1:nrow(superclass_list_pos)) {
	samp <- rownames(superclass_list_pos)[i]
	for (j in 1:ncol(superclass_list_pos)) {
		cl <- colnames(superclass_list_pos)[j]
		ft <- classifier_canopus_pos$Metabolite.name[which(classifier_canopus_pos$primary_class %in% classifier_canopus_pos$primary_class[grep(x=classifier_canopus_pos$primary_class, pattern=cl)]) ]
		ints <- as.numeric(feat_list_pos[i, which(colnames(feat_list_pos) %in% ft)])
		superclass_int_list_pos[i, j] <- sum(ints)
	}
}



# ---------- Classification at CHEMONT level of subclasses ----------
# Make subclasses at CHEMONT level 3
subclass_level_pos <- 3
div_subclasses_samples_names_pos <- NULL
for (i in c(1:subclass_level_pos)) {
	div_subclasses_samples_names_pos <- c(div_subclasses_samples_names_pos, lapply(X=strsplit(rownames(div_classes_samples_pos), '; '), FUN=function(x) { gsub(x=paste(x[1:i],sep='',collapse='; '),pattern='; NA',replacement='') }))
}
div_subclasses_samples_pos <- data.frame()
for (i in c(1:ncol(div_classes_samples_pos))) div_subclasses_samples_pos <- rbind(div_subclasses_samples_pos, rep(0, length(unique(div_subclasses_samples_names_pos))))
div_subclasses_samples_pos <- t(div_subclasses_samples_pos)
colnames(div_subclasses_samples_pos) <- colnames(div_classes_samples_pos)
rownames(div_subclasses_samples_pos) <- unique(div_subclasses_samples_names_pos)
for (i in rownames(div_subclasses_samples_pos)) {
	for (j in c(1:ncol(div_classes_samples_pos))) {
		div_subclasses_samples_pos[rownames(div_subclasses_samples_pos)==i, j] <- sum(div_classes_samples_pos[grep(x=rownames(div_classes_samples_pos), pattern=i), j])
	}
}

# Diversity of subclasses
div_subclasses_pos <- div_subclasses_samples_pos
div_subclasses_pos[is.na(div_subclasses_pos)] <- 0
div_subclasses_pos <- apply(X=div_subclasses_pos, MARGIN=1, FUN=function(x) { sum(x) })
div_subclasses_pos <- data.frame(row.names=names(div_subclasses_pos), frequency=as.numeric(div_subclasses_pos))

# Sunburst plot
pdf(file="plots/pos_ms2_subclasses_sunburst.pdf", encoding="ISOLatin1", pointsize=8, width=10, height=10, family="Helvetica")
sunBurstPlotFromSubstanceClasses(rownames(div_subclasses_pos), div_subclasses_pos$frequency, colorStart=0.0, colorAlpha=0.6)
dev.off()
write.csv(data.frame("Compound classes"=rownames(div_subclasses_pos), "Number of spectra"=div_subclasses_pos$frequency), file="plots/pos_ms2_subclasses_sunburst.csv", row.names=FALSE)

# Imputation of NA with zeros
div_subclasses_pos[is.na(div_subclasses_pos)] <- 0
div_subclasses_samples_pos[is.na(div_subclasses_samples_pos)] <- 0

# Classification list for statistics
subclass_list_pos <- as.data.frame(t(div_subclasses_samples_pos))
subclass_list_pos[is.na(subclass_list_pos)] <- 0

# Log Transform
#subclass_list_pos <- log2(subclass_list_pos + 1)

# Only keep subclass names
colnames(subclass_list_pos) <- gsub(x=colnames(subclass_list_pos), pattern='.*; ', replacement='')

# Generate subclass_int_list_pos with abundances instead of counts
subclass_int_list_pos <- subclass_list_pos

for (i in 1:nrow(subclass_list_pos)) {
	samp <- rownames(subclass_list_pos)[i]
	for (j in 1:ncol(subclass_list_pos)) {
		cl <- colnames(subclass_list_pos)[j]
		ft <- classifier_canopus_pos$Metabolite.name[which(classifier_canopus_pos$primary_class %in% classifier_canopus_pos$primary_class[grep(x=classifier_canopus_pos$primary_class, pattern=cl)]) ]
		ints <- as.numeric(feat_list_pos[i, which(colnames(feat_list_pos) %in% ft)])
		subclass_int_list_pos[i, j] <- sum(ints)
	}
}



# ---------- Diversity of NP classes ----------
# Create CANOPUS classifier object for each sample
npclassifiers_pos <- list()
for (i in mzml_names_pos) {
	obj <- names(which(bina_list_pos[rownames(bina_list_pos)==i, colnames(bina_list_pos) %in% classifier_canopus_pos$Metabolite.name] > 0))
	npclassifiers_pos[[i]] <- classifier_canopus_pos[classifier_canopus_pos$Metabolite.name %in% obj, ]
}

# Diversity of NPclasses per sample
div_npclasses_samples_pos <- NULL
for (i in mzml_names_pos) {
	obj <- gsub(x=apply(X=npclassifiers_pos[[i]][,c("NPC.pathway","NPC.superclass","NPC.class")], MARGIN=1, FUN=function(x) { paste(x, collapse="; ") }), pattern="(\\(|\\))", replacement="", perl=TRUE)
	obj <- table(obj)
	obj <- data.frame(classes=names(obj), frequency=as.numeric(obj))
	if (is.null(div_npclasses_samples_pos)) {
		div_npclasses_samples_pos <- obj
	} else {
		div_npclasses_samples_pos <- merge(div_npclasses_samples_pos, obj, by="classes", all.x=TRUE, all.y=TRUE)
	}
}
rownames(div_npclasses_samples_pos) <- div_npclasses_samples_pos$classes
div_npclasses_samples_pos <- div_npclasses_samples_pos[, -which(colnames(div_npclasses_samples_pos)=="classes")]
colnames(div_npclasses_samples_pos) <- mzml_names_pos

# Diversity of NPclasses
div_npclasses_pos <- div_npclasses_samples_pos
div_npclasses_pos[is.na(div_npclasses_pos)] <- 0
div_npclasses_pos <- apply(X=div_npclasses_pos, MARGIN=1, FUN=function(x) { sum(x) })
div_npclasses_pos <- data.frame(row.names=names(div_npclasses_pos), frequency=as.numeric(div_npclasses_pos))

# Plot diversity of NPclasses in all samples
pdf(file="plots/pos_ms2_npclasses_diversity.pdf", encoding="ISOLatin1", pointsize=8, width=6, height=14, family="Helvetica")
par(mfrow=c(1,1), mar=c(4,15,4,1), oma=c(0,0,0,0), cex.axis=0.8, cex=0.9)
barplot(div_npclasses_pos$frequency, names.arg=gsub('.*; ','',rownames(div_npclasses_pos)), las=1, horiz=TRUE, xlab="frequency", main="Diversity of compound classes", col=rainbow(n=nrow(div_npclasses_pos), alpha=0.6))
dev.off()

# Sunburst plot of NPclasses of all samples
pdf(file="plots/pos_ms2_npclasses_sunburst.pdf", encoding="ISOLatin1", pointsize=8, width=10, height=10, family="Helvetica")
par(mfrow=c(1,1), mar=c(0,0,0,0), oma=c(0,0,0,0), cex.axis=1, cex=1)
sunBurstPlotFromSubstanceClasses(rownames(div_npclasses_pos), div_npclasses_pos$frequency, colorStart=0.0, colorAlpha=0.6)
dev.off()
write.csv(data.frame("Compound classes"=rownames(div_npclasses_pos), "Number of spectra"=div_npclasses_pos$frequency), file="plots/pos_ms2_npclasses_sunburst.csv", row.names=FALSE)

# NPClasses
npclasses_pos <- rownames(div_npclasses_pos)



# ---------- Build NPclass diversity objects ----------
# Imputation of NA with zeros
div_npclasses_pos[is.na(div_npclasses_pos)] <- 0
div_npclasses_samples_pos[is.na(div_npclasses_samples_pos)] <- 0

# Classification list for statistics
npclass_list_pos <- as.data.frame(t(div_npclasses_samples_pos))
npclass_list_pos[is.na(npclass_list_pos)] <- 0

# Log Transform
#npclass_list_pos <- log2(npclass_list_pos + 1)

# Only keep class names
colnames(npclass_list_pos) <- gsub(x=colnames(npclass_list_pos), pattern='.*; ', replacement='')

# Generate npclass_int_list_pos with abundances instead of counts
npclass_int_list_pos <- npclass_list_pos

for (i in 1:nrow(npclass_list_pos)) {
	samp <- rownames(npclass_list_pos)[i]
	for (j in 1:ncol(npclass_list_pos)) {
		cl <- colnames(npclass_list_pos)[j]
		ft <- classifier_canopus_pos$Metabolite.name[which(gsub(x=classifier_canopus_pos$NPC.class, pattern='.*; ', replacement='') == cl)]
		ints <- as.numeric(feat_list_pos[i, which(colnames(feat_list_pos) %in% ft)])
		npclass_int_list_pos[i, j] <- sum(ints)
	}
}



# ---------- Classification at NP-Pathways ----------
# Extract NP-Pathways
div_nppathway_samples_names_pos <- NULL
div_nppathway_samples_names_pos <- c(div_nppathway_samples_names_pos, lapply(X=strsplit(rownames(div_npclasses_samples_pos), '; '), FUN=function(x) { gsub(x=paste(x[1],sep='',collapse='; '),pattern='; NA',replacement='') }))

div_nppathway_samples_pos <- data.frame()
for (i in c(1:ncol(div_npclasses_samples_pos))) div_nppathway_samples_pos <- rbind(div_nppathway_samples_pos, rep(0, length(unique(div_nppathway_samples_names_pos))))
div_nppathway_samples_pos <- t(div_nppathway_samples_pos)
colnames(div_nppathway_samples_pos) <- colnames(div_npclasses_samples_pos)
rownames(div_nppathway_samples_pos) <- unique(div_nppathway_samples_names_pos)
for (i in rownames(div_nppathway_samples_pos)) {
	for (j in c(1:ncol(div_npclasses_samples_pos))) {
		div_nppathway_samples_pos[rownames(div_nppathway_samples_pos)==i, j] <- sum(div_npclasses_samples_pos[grep(x=rownames(div_npclasses_samples_pos), pattern=i), j])
	}
}

# Diversity of nppathway
div_nppathway_pos <- div_nppathway_samples_pos
div_nppathway_pos[is.na(div_nppathway_pos)] <- 0
div_nppathway_pos <- apply(X=div_nppathway_pos, MARGIN=1, FUN=function(x) { sum(x) })
div_nppathway_pos <- data.frame(row.names=names(div_nppathway_pos), frequency=as.numeric(div_nppathway_pos))

# Sunburst plot
pdf(file="plots/pos_ms2_nppathway_sunburst.pdf", encoding="ISOLatin1", pointsize=8, width=10, height=10, family="Helvetica")
sunBurstPlotFromSubstanceClasses(rownames(div_nppathway_pos), div_nppathway_pos$frequency, colorStart=0.0, colorAlpha=0.6)
dev.off()
write.csv(data.frame("Compound classes"=rownames(div_nppathway_pos), "Number of spectra"=div_nppathway_pos$frequency), file="plots/pos_ms2_nppathway_sunburst.csv", row.names=FALSE)

# Imputation of NA with zeros
div_nppathway_pos[is.na(div_nppathway_pos)] <- 0
div_nppathway_samples_pos[is.na(div_nppathway_samples_pos)] <- 0

# Classification list for statistics
nppathway_list_pos <- as.data.frame(t(div_nppathway_samples_pos))
nppathway_list_pos[is.na(nppathway_list_pos)] <- 0

# Log Transform
#nppathway_list_pos <- log2(nppathway_list_pos + 1)

# Only keep nppathway names
colnames(nppathway_list_pos) <- gsub(x=colnames(nppathway_list_pos), pattern='.*; ', replacement='')

# Generate nppathway_int_list_pos with abundances instead of counts
nppathway_int_list_pos <- nppathway_list_pos

for (i in 1:nrow(nppathway_list_pos)) {
	samp <- rownames(nppathway_list_pos)[i]
	for (j in 1:ncol(nppathway_list_pos)) {
		cl <- colnames(nppathway_list_pos)[j]
		ft <- classifier_canopus_pos$Metabolite.name[which(classifier_canopus_pos$NPC.pathway %in% classifier_canopus_pos$NPC.pathway[grep(x=classifier_canopus_pos$NPC.pathway, pattern=cl)]) ]
		ints <- as.numeric(feat_list_pos[i, which(colnames(feat_list_pos) %in% ft)])
		nppathway_int_list_pos[i, j] <- sum(ints)
	}
}


