#!/usr/bin/env Rscript



# ---------- Preparations ----------
# Load libraries
library(parallel)               # Detect number of cpu cores
library(foreach)                # For multicore parallel
library(doMC)                   # For multicore parallel
library(RColorBrewer)           # For colors
library(multtest)               # For diffreport
library(MSnbase)                # MS features
library(xcms)                   # Swiss army knife for metabolomics
library(CAMERA)                 # Metabolite Profile Annotation
library(Spectra)                # Spectra package needed for XCMS3
library(mixOmics)               # Statistics for various Omics
# ERROR: vegan conflicts with mda and klaR, unload packages before using any of the analyses !!!
if ("package:mda" %in% search()) detach(package:mda, unload=TRUE)
if ("package:klaR" %in% search()) detach(package:klaR, unload=TRUE)
library(vegan)
library(multcomp)               # For Tukey test
library(Hmisc)                  # For correlation test
library(gplots)                 # For fancy heatmaps
library(lme4)                   # Linear Models with Random Effects
library(lmerTest)               # Create p-values for LMER
library(ape)                    # Phylogeny
library(pvclust)                # Phylogeny
library(dendextend)             # Phylogeny
library(cba)                    # Phylogeny
library(phangorn)               # Phylogeny
library(ontologyIndex)          # Reading obo ontology files
library(webchem)                # Converting InChI to InChIKey
library(circlize)               # For sunburst plot
library(plotrix)                # For sunburst plot
library(Boruta)                 # Random-Forest BORUTA
library(rpart)                  # Regression trees
library(caret)                  # Swiss-army knife for statistics
library(pROC)                   # Evaluation metrics
library(PRROC)                  # Evaluation metrics
library(multiROC)               # Evaluation metrics

# Setup R error handling to go to stderr
#options(show.error.messages=F, error=function() { cat(geterrmessage(), file=stderr()); q("no",1,F) } )

# Set locales and encoding
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
loc <- Sys.setlocale(category="LC_ALL", locale="C")
options(encoding="UTF-8")

# Set options
options(stringAsfactors=FALSE, useFancyQuotes=FALSE)

# Multicore parallel
#nSlaves <- detectCores(all.tests=FALSE, logical=FALSE)
#nSlaves <- 16
#registerDoMC(nSlaves)



# ---------- User variables ----------
# Take in trailing command line arguments
#args <- commandArgs(trailingOnly=TRUE)
#if (length(args) < 4) {
#    print("Error! No or not enough arguments given.")
#    print("Usage: $0 metfrag_results.csv metfrag_hits_limit metfrag_structure_method output.pdf")
#    quit(save="no", status=1, runLast=FALSE)
#}



# Working directory
setwd("~/Desktop/Projekte/Habilitation/Mosses/iestimate_mtbls709/")

# Data directory
mzml_dir <- "~/Desktop/Projekte/Habilitation/Mosses/ms2-dda/data/raw/"

# MS1 variables
polarity <- "positive"
pol <- substr(x=polarity, start=1, stop=3)
ppm <- 35
peakwidth <- c(4,21)
prefilter <- c(5,50)
#fitgauss <- FALSE
#verbose.columns <- FALSE
snthresh <- 11
mzdiff <- 0.01
noise <- 0
integrate <- 1
minfrac <- 0.7
bwindow <- 0.25
minfracspan <- 0.2
ms1_intensity_cutoff <- 500	          # approx. 0.01% (= log 100)

# Peak detection variables

# MS2 variables
mzabs <- 0.01                         # Absolute mass error (in seconds) used for merging MS/MS spectra
mzppm <- 5 #5                         # ppm error used for merging MS/MS spectra
rtabs <- 5 #20                        # Retention time error (in seconds) used for merging MS/MS spectra
max.rt.range <- 20                    # Permitted retention time window (in seconds) of grouped MS1 precursors
max.mz.range <- 0.01                  # Permitted m/z window of grouped MS1 precursors
min.rt <- 10                          # Minimum retention time for selected precursors
max.rt <- 1020                        # Maximum retention time for selected precursors
min.mz <- 50                          # Minimum m/z value for selected precursors
max.mz <- 1500                        # Maximum m/z value for selected precursors
msms.intensity.threshold <- 100       # Minimum intensity value for MS/MS peaks

# Classifier options
mzabs <- 0.01                        # Absolute mass error (in seconds) used for merging MS/MS spectra
mzppm <- 5 #35                       # ppm error used for merging MS/MS spectra
msms.intensity.threshold <- 100 #10  # Minimum intensity value for MS/MS peaks

CANOPUS <- TRUE                      # Use CANOPUS classification instead of MetFamily
SIRIUS_VERSION <- 5

# Preparations for plotting
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.8)

# Save and load RData
#save.image()
load(".RData")

# Preparations
source("_functions.r")

setMSnbaseFastLoad(TRUE)



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
mzml_files_pos <- mzml_files_pos[grep("MM8", mzml_files_pos, invert=TRUE)]

# Let's focus on one species
mzml_files_pos <- mzml_files_pos[grep("Marpol", mzml_files_pos, invert=FALSE, perl=TRUE)]

# Basenames of files without path and without extension
mzml_names_pos <- gsub('(.*)\\..*', '\\1', gsub('( |-|,)', '.', basename(mzml_files_pos)))

# Create phenodata based on species
mzml_pheno_pos <- as.factor(sapply(strsplit(as.character(mzml_names_pos), "_"), function(x) {
	nam <- x[3];
	nam;
}))
mzml_pheno_pos <- data.frame(sample_name=mzml_pheno_pos, sample_group=mzml_pheno_pos)
mzml_pheno_samples_pos <- as.factor(mzml_pheno_pos$sample_group)
mzml_pheno_colors_pos <- c("firebrick3", "chartreuse3", "darkgoldenrod3", "deepskyblue3")
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

# Restrict data to max.rt seconds
raw_data_pos <- filterRt(raw_data_pos, c(0, max.rt))

# Inspect mz values per file
raw_mz_pos <- mz(raw_data_pos)
raw_mz_pos <- split(raw_mz_pos, f = fromFile(raw_data_pos))
print(length(raw_mz_pos))

# Get base peak chromatograms
chromas_pos <- chromatogram(raw_data_pos, aggregationFun="max")

# Plot chromatograms based on phenodata groups
pdf(file="plots/pos_chromas.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
plot(chromas_pos, main="Raw chromatograms", xlab="retention time [s]", ylab="intensity", col=mzml_pheno_colors_samples_pos)
legend("topleft", bty="n", pt.cex=0.5, cex=0.7, y.intersp=0.7, text.width=0.5, pch=20, col=mzml_pheno_colors_pos, legend=mzml_pheno_samples_pos)
dev.off()

# Get TICs
pdf(file="plots/pos_tics.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
par(mfrow=c(1,1), mar=c(5,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
tics_pos <- split(tic(raw_data_pos), f=fromFile(raw_data_pos))
boxplot(tics_pos, names=mzml_names_pos, las=2, col=mzml_pheno_colors_pos, ylab="intensity", main="Total ion current")
dev.off()

# Grouping/binning the samples based on similarity of their base peak chromatogram to spot potentially problematic samples
chromas_bin_pos <- MSnbase::bin(chromas_pos, binSize=2)
chromas_bin_cor_pos <- cor(log2(do.call(cbind, lapply(chromas_bin_pos, intensity))))
colnames(chromas_bin_cor_pos) <- rownames(chromas_bin_cor_pos) <- raw_data_pos$sample_name
pdf(file="plots/pos_chromas_bin_cor.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
heatmap(chromas_bin_cor_pos)
dev.off()

# Assess retention times and intensities of first file
head(rtime(chromas_pos[1, 1]))
head(intensity(chromas_pos[1, 1]))

# Inspect peak width of standard compound for defining base peakwidth parameter below
pdf(file="plots/pos_chromas_standard_peakwidth.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
par(mfrow=c(2,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.6)
plot(chromatogram(raw_data_pos, mz=c(332.9, 335.1), rt=c(244, 251)), col=mzml_pheno_colors_samples_pos, main="EIC of perrottenine")
dev.off()

# Peak detection in MS1 data
if (polarity=="negative") {
	ms1_params_neg <- CentWaveParam(ppm=ppm, mzCenterFun="mean", peakwidth=peakwidth, prefilter=prefilter, mzdiff=mzdiff, snthresh=snthresh, noise=noise, integrate=integrate,
									firstBaselineCheck=TRUE, verboseColumns=TRUE, fitgauss=FALSE, roiList=list(), roiScales=numeric())
} else {
	ms1_params_pos <- CentWaveParam(ppm=ppm, mzCenterFun="mean", peakwidth=peakwidth, prefilter=prefilter, mzdiff=mzdiff, snthresh=snthresh, noise=noise, integrate=integrate,
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
	ms1_data_neg <- groupChromPeaks(ms1_data_neg, param=PeakDensityParam(sampleGroups=ms1_data_neg$sample_group, minFraction=minfrac, bw=bwindow))
} else {
	ms1_data_pos <- groupChromPeaks(ms1_data_pos, param=PeakDensityParam(sampleGroups=ms1_data_pos$sample_group, minFraction=minfrac, bw=bwindow))
}

# RT correction
if (polarity=="negative") {
	ms1_data_neg <- adjustRtime(ms1_data_neg, param=PeakGroupsParam(minFraction=minfrac,smooth="loess",span=minfracspan,family="gaussian"))
} else {
	ms1_data_pos <- adjustRtime(ms1_data_pos, param=PeakGroupsParam(minFraction=minfrac,smooth="loess",span=minfracspan,family="gaussian"))
}

# Plot the difference of raw and adjusted retention times
pdf(file="plots/pos_ms1_raw_adjusted.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=8, family="Helvetica")
par(mfrow=c(2,1), mar=c(4.5,4.2,4,1), cex=0.8)
plot(chromas_pos, peakType="none", main="Raw chromatograms", col=mzml_pheno_colors_samples_pos)
plotAdjustedRtime(ms1_data_pos, lwd=2, main="Retention Time correction", col=mzml_pheno_colors_samples_pos)
par(mfrow=c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0), cex.axis=0.9, cex=0.8)
dev.off()

# Group peaks
if (polarity=="negative") {
	ms1_data_neg <- groupChromPeaks(ms1_data_neg, param=PeakDensityParam(sampleGroups=ms1_data_neg$sample_group, minFraction=minfrac, bw=bwindow))
} else {
	ms1_data_pos <- groupChromPeaks(ms1_data_pos, param=PeakDensityParam(sampleGroups=ms1_data_pos$sample_group, minFraction=minfrac, bw=bwindow))
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
plot(ms1_pca_pos$x[, 1], ms1_pca_pos$x[,2], pch=19, main="PCA: Grouping of samples",
	 xlab=paste0("PC1: ", format(summary(ms1_pca_pos)$importance[2, 1] * 100, digits=3), " % variance"),
	 ylab=paste0("PC2: ", format(summary(ms1_pca_pos)$importance[2, 2] * 100, digits=3), " % variance"),
	 col=mzml_pheno_colors_samples_pos, cex=2)
grid()
text(ms1_pca_pos$x[, 1], ms1_pca_pos$x[,2], labels=ms1_data_pos$sample_name, col=mzml_pheno_colors_samples_pos, pos=3, cex=0.5)
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

# Missing value imputation
#feat_list_pos[is.na(feat_list_pos)] <- median(na.omit(as.numeric(unlist(feat_list_pos))))
#feat_list_pos[is.na(feat_list_pos)] <- 0

# Transform data
feat_list_pos <- log(feat_list_pos)

# Missing value imputation
#feat_list_pos[which(is.na(feat_list_pos))] <- median(na.omit(as.numeric(unlist(feat_list_pos))))

#feat_list_pos <- scale(feat_list_pos, scale=FALSE, center=FALSE)
feat_list_pos[is.na(feat_list_pos)] <- 0
feat_list_pos[which(feat_list_pos < 0)] <- 0
feat_list_pos[is.infinite(feat_list_pos)] <- 0
feat_list_pos <- feat_list_pos[!apply(feat_list_pos, MARGIN=1, function(x) max(x,na.rm=TRUE) == min(x,na.rm=TRUE)),]

# Plot histogram
pdf(file="plots/pos_feat_list_hist.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
hist(as.numeric(feat_list_pos), main="Histogram of feature table")
dev.off()

# PCA
pdf(file="plots/pos_ms1_feature_table_pca12.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
ms1_pca_pos <- prcomp(feat_list_pos, center=TRUE)
plot(ms1_pca_pos$x[, 1], ms1_pca_pos$x[,2], pch=19, main="PCA of feature table",
	 xlab=paste0("PC1: ", format(summary(ms1_pca_pos)$importance[2, 1] * 100, digits=3), " % variance"),
	 ylab=paste0("PC2: ", format(summary(ms1_pca_pos)$importance[2, 2] * 100, digits=3), " % variance"),
	 col=mzml_pheno_colors_samples_pos, cex=2)
grid()
text(ms1_pca_pos$x[, 1], ms1_pca_pos$x[, 2], labels=ms1_data_pos$sample_name, col=mzml_pheno_colors_samples_pos, pos=3, cex=0.5)
dev.off()

pdf(file="plots/pos_ms1_feature_table_pca23.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
ms1_pca_pos <- prcomp(feat_list_pos, center=TRUE)
plot(ms1_pca_pos$x[, 2], ms1_pca_pos$x[,3], pch=19, main="PCA of feature table",
	 xlab=paste0("PC3: ", format(summary(ms1_pca_pos)$importance[2, 2] * 100, digits=3), " % variance"),
	 ylab=paste0("PC4: ", format(summary(ms1_pca_pos)$importance[2, 3] * 100, digits=3), " % variance"),
	 col=mzml_pheno_colors_samples_pos, cex=2)
grid()
text(ms1_pca_pos$x[, 2], ms1_pca_pos$x[, 3], labels=ms1_data_pos$sample_name, col=mzml_pheno_colors_samples_pos, pos=3, cex=0.5)
dev.off()

# Create single 0/1 matrix
#bina_list_pos <- t(ms1_matrix_pos)
#bina_list_pos[is.na(bina_list_pos)] <- 1
#bina_list_pos <- log(bina_list_pos)
bina_list_pos <- feat_list_pos
bina_list_pos[bina_list_pos < log(ms1_intensity_cutoff)] <- 0
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

# Remove NAs if present
model_div_pos[is.na(model_div_pos)] <- 0



# ############################## MS2 ##############################



# ---------- MS2 spectra detection ----------
# Estimate precursor intensity
precursor_intensity_pos <- estimatePrecursorIntensity(ms1_data_pos)
print(head(na.omit(precursor_intensity_pos)))

# Reconstruct MS2 spectra from MS1 data
ms2_data_pos <- chromPeakSpectra(ms1_data_pos, msLevel=2L, return.type="Spectra")
print(ms2_data_pos)
print(length(ms2_data_pos$peak_id))

# Extract all usable MS2 spectra
ms2_spectra_pos <- list()
#for (i in 1:nrow(ms1_def_pos)) {
ms2_spectra_pos <- foreach(i=1:nrow(ms1_def_pos)) %dopar% {
	# Extract existing MS2 spectra for feature
	feature_of_interest <- ms1_def_pos[i, "mzmed"]
	peaks_of_interest <- chromPeaks(ms1_data_pos, mz=feature_of_interest, ppm=ppm)
	
	# Continue if feature has MS2 peaks
	if (length(which(ms2_data_pos$peak_id %in% rownames(peaks_of_interest))) > 0) {
		# Extract spectra
		spectra_of_interest <- ms2_data_pos[ms2_data_pos$peak_id %in% rownames(peaks_of_interest)]
		combined_spectra_of_interest <- filterIntensity(spectra_of_interest, intensity=c(1,Inf), backend=MsBackendDataFrame)
		combined_spectra_of_interest <- setBackend(combined_spectra_of_interest, backend=MsBackendDataFrame())
		
		# Combine spectra
		combined_spectra_of_interest <- Spectra::combineSpectra(combined_spectra_of_interest, FUN=combinePeaks, ppm=ppm, peaks="union", minProp=0.8, intensityFun=median, mzFun=median, backend=MsBackendDataFrame)#f=rownames(peaks_of_interest))
		
		# Remove noise from spectra
		#combined_spectra_of_interest <- pickPeaks(combined_spectra_of_interest, snr=1.0, method="SuperSmoother") #MAD
		#combined_spectra_of_interest <- Spectra::smooth(combined_spectra_of_interest, method="SavitzkyGolay") #(Weighted)MovingAverage
		
		# Only keep spectral data
		combined_spectra_peaks <- as.data.frame(peaksData(combined_spectra_of_interest)[[1]])
		
		# Plot merged spectrum
		#Spectra::plotSpectra(combined_spectra_of_interest)
		#plot(x=combined_spectra_peaks[,1], y=combined_spectra_peaks[,2], type="h", xlab="m/z", ylab="intensity", main=paste("Precursor m/z",combined_spectra_of_interest@backend@spectraData$precursorMz[[1]]))
		#length(spectra_of_interest$peak_id)
		
		#ms2_spectra_pos[[i]] <- combined_spectra_of_interest
		return(combined_spectra_of_interest)
	}
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



# ---------- Classify MS2 spectra with CANOPUS ----------
# Apply classified classes onto feature table
ms1_def_pos$primary_class <- ""
ms1_def_pos$alternative_classes <- ""

# Read SIRIUS/CANOPUS classifier
classifier_canopus_pos <- read.table(file=paste0("data/pos_ms2_canopus_summary.tsv"), header=TRUE, sep="\t", quote="\"", fill=FALSE, dec=".", stringsAsFactors=FALSE)#, encoding="UTF-8", fileEncoding="UTF-8")
if (SIRIUS_VERSION == 4) {
	classifier_canopus_pos$Metabolite.name <- gsub(x=classifier_canopus_pos$name, pattern=".*_", replacement="")
	classifier_canopus_pos$primary_class <- paste("Organic compounds", classifier_canopus_pos$superclass, classifier_canopus_pos$class, classifier_canopus_pos$subclass, classifier_canopus_pos$level.5, sep="; ")
	classifier_canopus_pos$primary_class <- gsub(x=classifier_canopus_pos$primary_class, pattern="; ; ", replacement="; ", perl=TRUE)
	classifier_canopus_pos$primary_class <- gsub(x=classifier_canopus_pos$primary_class, pattern="; ; ", replacement="; ", perl=TRUE)
	classifier_canopus_pos$primary_class <- gsub(x=classifier_canopus_pos$primary_class, pattern="; $", replacement="", perl=TRUE)
	classifier_canopus_pos$primary_class <- gsub(x=classifier_canopus_pos$primary_class, pattern="(\\'|\\>|\\(|\\))", replacement="", perl=TRUE)
	classifier_canopus_pos$Annotation..putative. <- classifier_canopus_pos$primary_class
	classifier_canopus_pos$alternative_classes <- classifier_canopus_pos$all.classifications
} else {
	sirius_identification_pos <- read.table(file=paste0("data/pos_ms2_compound_identifications.tsv"), header=TRUE, sep="\t", quote="\"", fill=FALSE, dec=".", stringsAsFactors=FALSE)#, encoding="UTF-8", fileEncoding="UTF-8")
	classifier_canopus_pos$Metabolite.name <- ""
	for (i in 1:length(classifier_canopus_pos$id)) {
		x = which(gsub(x=sirius_identification_pos$id, pattern=".*?_", replacement="", perl=TRUE) %in% gsub(x=classifier_canopus_pos$id[i], pattern=".*?_", replacement="", perl=TRUE))
		if (length(x) > 0) classifier_canopus_pos$Metabolite.name[i] <- gsub(x=sirius_identification_pos$id[x], pattern=".*?_", replacement="", perl=TRUE)
	}
	classifier_canopus_pos$Metabolite.name[classifier_canopus_pos$Metabolite.name == "null"] <- ""
	classifier_canopus_pos$primary_class <- paste("Organic compounds", classifier_canopus_pos$ClassyFire_superclass, classifier_canopus_pos$ClassyFire_class, classifier_canopus_pos$ClassyFire_subclass, classifier_canopus_pos$ClassyFire_level.5, sep="; ")
	classifier_canopus_pos$primary_class <- gsub(x=classifier_canopus_pos$primary_class, pattern="; ; ", replacement="; ", perl=TRUE)
	classifier_canopus_pos$primary_class <- gsub(x=classifier_canopus_pos$primary_class, pattern="; ; ", replacement="; ", perl=TRUE)
	classifier_canopus_pos$primary_class <- gsub(x=classifier_canopus_pos$primary_class, pattern="; $", replacement="", perl=TRUE)
	classifier_canopus_pos$primary_class <- gsub(x=classifier_canopus_pos$primary_class, pattern="(\\'|\\>|\\(|\\))", replacement="", perl=TRUE)
	classifier_canopus_pos$Annotation..putative. <- classifier_canopus_pos$primary_class
	classifier_canopus_pos$alternative_classes <- classifier_canopus_pos$all.classifications
}

# Store primary class in ms1_def_pos
for (j in unique(classifier_canopus_pos$Metabolite.name)) {
	obj <- classifier_canopus_pos[classifier_canopus_pos$Metabolite.name %in% j, "Annotation..putative."]
	ms1_def_pos[j, "primary_class"] <- obj[1]
}

# Build comp_list
comp_list_pos <- feat_list_pos[, c(na.omit(rownames(ms1_def_pos)[ms1_def_pos$has_ms2==1]))]

# Store identification in ms1_def_pos
ms1_def_pos[, "identification"] <- ""
for (i in 1:nrow(sirius_identification_pos)) {
	if ((sirius_identification_pos$name[i] != "") & (sirius_identification_pos$name[i] != "null")) {
		ms1_def_pos[gsub(x=sirius_identification_pos$id, pattern=".*?_", replacement="", perl=TRUE)[i], "identification"] <- sirius_identification_pos$name[i]
	}
}



# ---------- Annotate MS2 spectra for samples ----------
# Write MSP for all samples
foreach (j=mzml_names_pos) %dopar% {
	#for (j in mzml_names_pos) {
	# Iterate through all MS2 spectra found in the sample
	msp_text <- NULL
	alignment_id <- 0
	for (i in names(ms2_spectra_pos)) {
		if (bina_list_pos[j,i] != 0) {
			alignment_id <- alignment_id + 1
			
			msp_text <- c(msp_text, paste("NAME:", i))
			msp_text <- c(msp_text, paste("AlignmentID:", alignment_id))
			msp_text <- c(msp_text, paste("RETENTIONTIME:", ms1_def_pos[i, "rtmed"]))
			msp_text <- c(msp_text, paste("PRECURSORMZ:", ms1_def_pos[i, "mzmed"]))
			msp_text <- c(msp_text, paste("METABOLITENAME:", i))
			msp_text <- c(msp_text, paste("COMPOUND_CLASS:", ms1_def_pos[i, "primary_class"]))
			msp_text <- c(msp_text, paste("ALTERNATIVE_CLASSES:", ms1_def_pos[i, "alternative_classes"]))
			if (polarity == "positive") {
				msp_text <- c(msp_text, paste("ADDUCTIONNAME:", "[M+H]+"))
			} else {
				msp_text <- c(msp_text, paste("ADDUCTIONNAME:", "[M-H]-"))
			}
			msp_text <- c(msp_text, paste("NumPeaks:", nrow(as.data.frame(peaksData(ms2_spectra_pos[[i]])[[1]]))))
			msp_text <- c(msp_text, paste(as.data.frame(peaksData(ms2_spectra_pos[[i]])[[1]])$mz, as.data.frame(peaksData(ms2_spectra_pos[[i]])[[1]])$intensity, sep="\t"))
			msp_text <- c(msp_text, "")
		}
	}
	
	# Write MSP file
	cat(msp_text, file=paste0("msp/",j,".msp"), sep="\n")
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
		class_int_list_pos[i, j] <- sum(ints) * as.numeric(class_list_pos[i, j])
	}
}



# ---------- Classification at level of superclasses ----------
if (FALSE) {
	# Use CANOPUS Metabolite.name classification
	classifier_canopus_pos$primary_class <- paste("Organic compounds", classifier_canopus_pos$superclass, classifier_canopus_pos$most.specific.class, sep="; ")
	classifier_canopus_pos$primary_class <- gsub(x=classifier_canopus_pos$primary_class, pattern="; ; ", replacement="; ", perl=TRUE)
	classifier_canopus_pos$primary_class <- gsub(x=classifier_canopus_pos$primary_class, pattern="; $", replacement="", perl=TRUE)
	classifier_canopus_pos$primary_class <- gsub(x=classifier_canopus_pos$primary_class, pattern="(\\'|\\>|\\(|\\))", replacement="", perl=TRUE)
	classifier_canopus_pos$Annotation..putative. <- classifier_canopus_pos$primary_class
	
	# Create CANOPUS classifier object for each sample
	superclassifiers_pos <- list()
	for (i in mzml_names_pos) {
		obj <- names(which(bina_list_pos[rownames(bina_list_pos)==i, colnames(bina_list_pos) %in% classifier_canopus_pos$Metabolite.name] > 0))
		superclassifiers_pos[[i]] <- classifier_canopus_pos[classifier_canopus_pos$Metabolite.name %in% obj, ]
	}
	
	# Diversity of classes per sample
	div_superclasses_samples_pos <- NULL
	for (i in mzml_names_pos) {
		obj <- table(superclassifiers_pos[[i]][,"Annotation..putative."])
		obj <- data.frame(classes=names(obj), frequency=as.numeric(obj))
		if (is.null(div_superclasses_samples_pos)) {
			div_superclasses_samples_pos <- obj
		} else {
			div_superclasses_samples_pos <- merge(div_superclasses_samples_pos, obj, by="classes", all.x=TRUE, all.y=TRUE)
		}
	}
	rownames(div_superclasses_samples_pos) <- div_superclasses_samples_pos$classes
	div_superclasses_samples_pos <- div_superclasses_samples_pos[, -which(colnames(div_superclasses_samples_pos)=="classes")]
	colnames(div_superclasses_samples_pos) <- mzml_names_pos
} else {
	# Make superclasses at level 3
	superclass_level_pos <- 3
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
		superclass_int_list_pos[i, j] <- sum(ints) * as.numeric(superclass_list_pos[i, j])
	}
}



# ############################## MS1 statistics ##############################



# ---------- Histogram ----------
pdf(file="plots/pos_ms1_hist.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
hist(as.numeric(feat_list_pos))
dev.off()



# ---------- Variation partitioning ----------
model_varpart_pos <- varpart(scale(feat_list_pos), ~ mzml_pheno_samples_pos, ~ mzml_pheno_samples_pos)

# Plot results
pdf(file="plots/pos_ms1_varpart.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
plot(model_varpart_pos, Xnames=c("samples","species"), cutoff=0, cex=1.2, id.size=1.2, digits=1, bg=c("blue","green"))
dev.off()



# ---------- PCA ----------
model_pca_pos <- prcomp(feat_list_pos)

pdf(paste("plots/pos_ms1_pca.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
plot(model_pca_pos$x[, 1], model_pca_pos$x[,2], pch=19, main="PCA",
	 xlab=paste0("PC1: ", format(summary(model_pca_pos)$importance[2, 2] * 100, digits=3), " % variance"),
	 ylab=paste0("PC2: ", format(summary(model_pca_pos)$importance[2, 3] * 100, digits=3), " % variance"),
	 col=mzml_pheno_colors_samples_pos, bg=mzml_pheno_colors_samples_pos, cex=2)
grid()
text(model_pca_pos$x[, 1], model_pca_pos$x[,2], labels=mzml_names_pos, col=mzml_pheno_colors_samples_pos, pos=3, cex=0.5)
dev.off()



# ---------- Diversity ----------
# Plot unique features
pdf(paste("plots/pos_ms1_richness_unique.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(model_div_pos$unique ~ mzml_pheno_samples_pos, col=mzml_pheno_colors_pos, names=NA, main="Number of unique features", xlab="treatment", ylab="number of unique features")
text(1:length(unique(mzml_pheno_samples_pos)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=levels(mzml_pheno_samples_pos), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=model_div_pos$unique, term=as.factor(mzml_pheno_samples_pos))
text(1:length(unique(mzml_pheno_samples_pos)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# Plot total features
pdf(paste("plots/pos_ms1_richness_total.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(model_div_pos$features ~ mzml_pheno_samples_pos, col=mzml_pheno_colors_pos, names=NA, main="Number of features", xlab="treatment", ylab="number of features")
text(1:length(unique(mzml_pheno_samples_pos)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=levels(mzml_pheno_samples_pos), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=model_div_pos$features, term=as.factor(mzml_pheno_samples_pos))
text(1:length(unique(mzml_pheno_samples_pos)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# Plot Shannon index
pdf(paste("plots/pos_ms1_diversity_shannon.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(model_div_pos$shannon ~ mzml_pheno_samples_pos, col=mzml_pheno_colors_pos, names=NA, main="Shannon diversity (H\')", xlab="treatment", ylab="Shannon diversity index (H\')")
text(1:length(unique(mzml_pheno_samples_pos)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=levels(mzml_pheno_samples_pos), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=model_div_pos$shannon, term=as.factor(mzml_pheno_samples_pos))
text(1:length(unique(mzml_pheno_samples_pos)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()

# Plot Pielou evenness for species
pdf(paste("plots/pos_ms1_diversity_pielou.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
boxplot(model_div_pos$pielou ~ mzml_pheno_samples_pos, col=mzml_pheno_colors_pos, names=NA, main="Pielou\'s evenness", xlab="treatment", ylab="Pielou diversity index (J)")
text(1:length(unique(mzml_pheno_samples_pos)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/14, srt=-22.5, adj=0.5, labels=levels(mzml_pheno_samples_pos), xpd=TRUE, cex=0.9)
div_tukey <- tukey.test(response=model_div_pos$pielou, term=as.factor(mzml_pheno_samples_pos))
text(1:length(unique(mzml_pheno_samples_pos)), par("usr")[4]+(par("usr")[4]-par("usr")[3])/40, adj=0.5, labels=div_tukey[,1], xpd=TRUE, cex=0.8)
dev.off()



# ---------- Variable Selection ----------
# Random Forest
sel_rf_pos <- f.select_features_random_forest(feat_matrix=feat_list_pos, sel_factor=as.factor(mzml_pheno_samples_pos), sel_colors=mzml_pheno_colors_pos, tune_length=10, quantile_threshold=0.90, plot_roc_filename="plots/pos_ms1_select_rf_roc.pdf")
print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_rf_pos$`_selected_variables_`)))
f.heatmap.selected_features(feat_list=feat_list_pos, sel_feat=sel_rf_pos$`_selected_variables_`, sel_names=paste0("     ",sel_rf_pos$`_selected_variables_`), filename="plots/pos_ms1_select_rf.pdf", main="Random Forest", plot_width=6, plot_height=5, cex_col=0.5, cex_row=0.4)

# PLS
sel_pls_pos <- f.select_features_pls(feat_matrix=feat_list_pos, sel_factor=as.factor(mzml_pheno_samples_pos), sel_colors=mzml_pheno_colors_pos, components=length(unique(mzml_pheno_samples_pos))-1, tune_length=10, quantile_threshold=0.982, plot_roc_filename="plots/pos_ms1_select_pls_roc.pdf")
print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_pls_pos$`_selected_variables_`)))
f.heatmap.selected_features(feat_list=feat_list_pos, sel_feat=sel_pls_pos$`_selected_variables_`, sel_names=paste0("     ",sel_pls_pos$`_selected_variables_`), filename="plots/pos_ms1_select_pls.pdf", main="PLS", plot_width=6, plot_height=5, cex_col=0.5, cex_row=0.4)

# PLS with comp_list
sel_pls_pos <- f.select_features_pls(feat_matrix=comp_list_pos, sel_factor=as.factor(mzml_pheno_samples_pos), sel_colors=mzml_pheno_colors_pos, components=length(unique(mzml_pheno_samples_pos))-1, tune_length=10, quantile_threshold=0.982, plot_roc_filename="plots/pos_ms1_comp_select_pls_roc.pdf")
print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_pls_pos$`_selected_variables_`)))
sel_var_pos <- sel_pls_pos$`_selected_variables_`
for (i in 1:length(sel_var_pos)) {
	if (ms1_def_pos[sel_var_pos[i], "identification"] != "") sel_var_pos[i] <- ms1_def_pos[sel_var_pos[i], "identification"]
	if (nchar(sel_var_pos[i]) > 20) sel_var_pos[i] <- substr(sel_var_pos[i], 1, 20)
}
f.heatmap.selected_features(feat_list=comp_list_pos, sel_feat=sel_pls_pos$`_selected_variables_`, sel_names=sel_var_pos, filename="plots/pos_ms1_comp_select_pls.pdf", main="PLS", plot_width=6, plot_height=5, cex_col=0.3, cex_row=0.4)



# ############################## MS2 Classes statistics ##############################



# ---------- Histogram ----------
pdf(file="plots/pos_ms2_classes_hist.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
hist(as.numeric(unlist(class_list_pos)), main="Histogram", xlab="class_list")
dev.off()



# ---------- Variation partitioning ----------
model_varpart_pos <- varpart(class_list_pos, ~ mzml_pheno_samples_pos, ~ mzml_pheno_samples_pos)

# Plot results
pdf(file="plots/pos_ms2_classes_varpart.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
plot(model_varpart_pos, Xnames=c("sample group","samples"), cutoff=0, cex=1.2, id.size=1.2, digits=1, bg=c("blue","green"))
dev.off()



# ---------- PCA ----------
model_pca_pos <- prcomp(class_list_pos)

pdf(paste("plots/pos_ms2_classes_pca.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
plot(model_pca_pos$x[, 1], model_pca_pos$x[, 2], pch=19, main="PCA",
	 xlab=paste0("PC2: ", format(summary(model_pca_pos)$importance[2, 1] * 100, digits=3), " % variance"),
	 ylab=paste0("PC3: ", format(summary(model_pca_pos)$importance[2, 2] * 100, digits=3), " % variance"),
	 col=mzml_pheno_colors_samples_pos, bg=mzml_pheno_colors_samples_pos, cex=2)
grid()
text(model_pca_pos$x[, 1], model_pca_pos$x[, 2], labels=mzml_names_pos, col=mzml_pheno_colors_samples_pos, pos=3, cex=0.5)
dev.off()


# ---------- Feature selection of classes ----------
# Tukey on ANOVA
#sel_aov_pos <- f.select_features_aov(feat_class=mzml_pheno_samples_pos, feat_list=class_list_pos, conf_level=0.95, pvalue=0.99)
#print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_aov_pos)))
#f.heatmap.selected_features(feat_list=class_list_pos, sel_feat=sel_aov_pos, filename="plots/pos_ms2_classes_select_aov.pdf", main="Tukey on ANOVA")

# Linear Model
#sel_lm_pos <- f.select_features_linear_model(feat_matrix=class_list_pos, sel_factor=mzml_pheno_samples_pos, sel_formula=feature ~ mzml_pheno_samples_pos, pvalue=0.05)
#print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_lm_pos)))
#f.heatmap.selected_features(feat_list=class_list_pos, sel_feat=sel_lm_pos, filename="plots/pos_ms2_classes_select_lm.pdf", main="Linear Model")

# Random Forest
sel_rf_pos <- f.select_features_random_forest(feat_matrix=class_list_pos, sel_factor=mzml_pheno_samples_pos, sel_colors=mzml_pheno_colors_pos, tune_length=10, quantile_threshold=0.95, plot_roc_filename="plots/pos_ms2_classes_select_rf_roc.pdf")
print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_rf_pos$`_selected_variables_`)))
f.heatmap.selected_features(feat_list=class_list_pos, sel_feat=sel_rf_pos$`_selected_variables_`, sel_names=paste0("     ",sel_rf_pos$`_selected_variables_`), filename="plots/pos_ms2_classes_select_rf.pdf", main="Random Forest", plot_width=6, plot_height=5, cex_col=0.5, cex_row=0.4)

# PLS
sel_pls_pos <- f.select_features_pls(feat_matrix=class_list_pos, sel_factor=mzml_pheno_samples_pos, sel_colors=mzml_pheno_colors_pos, components=length(unique(mzml_pheno_samples_pos))-1, tune_length=10, quantile_threshold=0.9, plot_roc_filename="plots/pos_ms2_classes_select_pls_roc.pdf")
print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_pls_pos$`_selected_variables_`)))
f.heatmap.selected_features(feat_list=class_list_pos, sel_feat=sel_pls_pos$`_selected_variables_`, sel_names=paste0("     ",sel_pls_pos$`_selected_variables_`), filename="plots/pos_ms2_classes_select_pls.pdf", main="PLS", plot_width=6, plot_height=5, cex_col=0.5, cex_row=0.4)



# ############################## MS2 Superclasses statistics ##############################



# ---------- Histogram ----------
pdf(file="plots/pos_ms2_superclasses_hist.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
hist(as.numeric(unlist(superclass_list_pos)), main="Histogram", xlab="superclass_list")
dev.off()



# ---------- Variation partitioning ----------
model_varpart_pos <- varpart(superclass_list_pos, ~ mzml_pheno_samples_pos, ~ mzml_pheno_samples_pos)

# Plot results
pdf(file="plots/pos_ms2_superclasses_varpart.pdf", encoding="ISOLatin1", pointsize=10, width=6, height=4, family="Helvetica")
plot(model_varpart_pos, Xnames=c("sample group","samples"), cutoff=0, cex=1.2, id.size=1.2, digits=1, bg=c("blue","green"))
dev.off()



# ---------- PCA ----------
model_pca_pos <- prcomp(superclass_list_pos)

pdf(paste("plots/pos_ms2_superclasses_pca.pdf",sep=""), encoding="ISOLatin1", pointsize=10, width=5, height=5, family="Helvetica")
plot(model_pca_pos$x[, 1], model_pca_pos$x[, 2], pch=19, main="PCA",
	 xlab=paste0("PC2: ", format(summary(model_pca_pos)$importance[2, 1] * 100, digits=3), " % variance"),
	 ylab=paste0("PC3: ", format(summary(model_pca_pos)$importance[2, 2] * 100, digits=3), " % variance"),
	 col=mzml_pheno_colors_samples_pos, bg=mzml_pheno_colors_samples_pos, cex=2)
grid()
text(model_pca_pos$x[, 1], model_pca_pos$x[, 2], labels=mzml_names_pos, col=mzml_pheno_colors_samples_pos, pos=3, cex=0.5)
dev.off()



# ---------- Feature selection of superclasses ----------
# Tukey on ANOVA
#sel_aov_pos <- f.select_features_aov(feat_class=mzml_pheno_samples_pos, feat_list=superclass_list_pos, conf_level=0.95, pvalue=0.99)
#print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_aov_pos)))
#f.heatmap.selected_features(feat_list=superclass_list_pos, sel_feat=sel_aov_pos, filename="plots/pos_ms2_superclasses_select_aov.pdf", main="Tukey on ANOVA")

# Linear Model
#sel_lm_pos <- f.select_features_linear_model(feat_matrix=superclass_list_pos, sel_factor=mzml_pheno_samples_pos, sel_formula=feature ~ mzml_pheno_samples_pos, pvalue=0.05)
#print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_lm_pos)))
#f.heatmap.selected_features(feat_list=superclass_list_pos, sel_feat=sel_lm_pos, filename="plots/pos_ms2_superclasses_select_lm.pdf", main="Linear Model")

# Random Forest
sel_rf_pos <- f.select_features_random_forest(feat_matrix=superclass_list_pos, sel_factor=mzml_pheno_samples_pos, sel_colors=mzml_pheno_colors_pos, tune_length=10, quantile_threshold=0.95, plot_roc_filename="plots/pos_ms2_superclasses_select_rf_roc.pdf")
print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_rf_pos$`_selected_variables_`)))
f.heatmap.selected_features(feat_list=superclass_list_pos, sel_feat=sel_rf_pos$`_selected_variables_`, filename="plots/pos_ms2_superclasses_select_rf.pdf", main="Random Forest", plot_width=6, plot_height=5, cex_col=0.5, cex_row=0.4)

# PLS
sel_pls_pos <- f.select_features_pls(feat_matrix=superclass_list_pos, sel_factor=mzml_pheno_samples_pos, sel_colors=mzml_pheno_colors_pos, components=length(unique(mzml_pheno_samples_pos))-1, tune_length=10, quantile_threshold=0.95, plot_roc_filename="plots/pos_ms2_superclasses_select_pls_roc.pdf")
print(paste("Number of selected features:", f.count.selected_features(sel_feat=sel_pls_pos$`_selected_variables_`)))
f.heatmap.selected_features(feat_list=superclass_list_pos, sel_feat=sel_pls_pos$`_selected_variables_`, filename="plots/pos_ms2_superclasses_select_pls.pdf", main="PLS", plot_width=6, plot_height=5, cex_col=0.5, cex_row=0.4)



