


# ---------- Shannon Diversity ----------
#' Function to calculate the Shannon diversity measure
#'
#' @param p vector with response variables of one sample
#' @export
#' @import vegan MASS
#' @examples
#' shannon.diversity(p=c(4,8))
shannon.diversity <- function(p) {
	# Based on Li et al. (2016)
	# Function is obsolete, as it returns same values than vegan::diversity(x, index="shannon")
	# Since we do not expect features with negative intensity,
	# we exclude negative values and take the natural logarithm
	if (min(p) < 0 || sum(p) <= 0)
		return(NA)
	pij <- p[p>0] / sum(p)
	-sum(pij * log(pij))
}



# ---------- Menhinick Diversity ----------
#' Function to calculate the Menhinick diversity measure
#'
#' @param p vector with response variables of one sample
#' @export
#' @import vegan MASS
#' @examples
#' menhinick.diversity(p=c(4,8))
menhinick.diversity <- function(p) {
	# Based on: http://www.coastalwiki.org/wiki/Measurements_of_biodiversity#Species_richness_indices
	D_Mn <- length(p) / sqrt(vegan::specnumber(p))
}



# ---------- Tukey-Test ----------
#' Function to perform the Tukey post-hoc HSD test on response and model terms
#'
#' @param response vector with response variables
#' @param term vector with factorized terms
#' @export
#' @import stats multcomp
###_ @examples
###_ tukey.test(response=model_div$unique, term=as.factor(mzml_pheno_samples))
tukey.test <- function(response, term) {
	model_anova <- stats::aov(formula(response ~ term))
	model_mc <- multcomp::glht(model_anova, multcomp::mcp(term="Tukey"))
	model_cld <- multcomp::cld(summary(model_mc), decreasing=TRUE)
	model_tukey <- data.frame("tukey_groups"=model_cld$mcletters$Letters)
	return(model_tukey)
}



# ---------- P-Values ----------
print_p.values <- function(p.values) {
	p.values[1] <- 1
	p.values[p.values < 0.001] <- '***'
	p.values[(p.values >= 0.001 & p.values < 0.01)] <- '**'
	p.values[(p.values >= 0.01 & p.values < 0.05)] <- '*'
	p.values[(p.values >= 0.05 & p.values < 0.1)] <- '.'
	p.values[p.values >= 0.1] <- ' '
	return(p.values)
}



# ---------- Calculate R-squared ----------
f.r2 <- function(actual, predicted) {
	actual <- as.numeric(actual)
	predicted <- as.numeric(predicted)

	R2 <- caret::postResample(actual, predicted)[2]
	#R2 <- 1 - (sum((actual-predicted)^2)/sum((actual-mean(actual))^2))

	return(R2)
}



# ---------- Calculate Classification Rate ----------
f.classification_rate <- function(sel_factor, predicted) {
	accuracy <- table(predicted, sel_factor)
	cr <- sum(diag(accuracy))/sum(accuracy)

	#print(paste0("Classification rate: ", round(cr ,3)))

	return(cr)
}



# ---------- Calculate Accuracy ----------
f.accuracy <- function(sel_factor, predicted) {
	return(length(which(predicted == sel_factor)) / length(sel_factor))
}



# ---------- Calculate Weighted Accuracy ----------
f.weighted_accuracy <- function(sel_factor, predicted) {
	# Calculate weights
	weights <- rep(1 / length(levels(sel_factor)), length(levels(sel_factor)))

	# Check
	sel_levels <- levels(sel_factor)
	if (length(weights) != length(sel_levels)) {
		stop("Error! Number of weights should have some length as the number of classes.")
	}
	if (sum(weights) != 1) {
		stop("Error! Weights do not sum to 1.")
	}

	# Calculate
	accuracies <- lapply(sel_levels, function(x) {
		idx <- which(sel_factor==x)
		return(f.accuracy(predicted[idx], sel_factor[idx]))
	})
	accuracies <- mean(unlist(accuracies))

	return(accuracies)
}



# ---------- Performance Measures ----------
f.performance_measures_caret <- function(model, sel_factor, sel_colors) {
	# List of measures
	sel <- list()

	# Select index with highest accuracy from prediction model, only when not using savePredictions="final"
	#sel_ind_acc <- which(model$results$Accuracy==max(model$results$Accuracy))
	#if (! is.null(names(model$pred$mtry))) {
	#	sel_ind <- model$pred$mtry == unique(model$pred$mtry)[sel_ind_acc]
	#} else {
	#	sel_ind <- c(1:nrow(model$pred))
	#}

	# Performance measures
	sel_levels <- levels(as.factor(sel_factor))

	# Build prediction matrix
	model_pred <- as.data.frame(model$pred[, which(colnames(model$pred) %in% sel_levels)])
	colnames(model_pred) <- paste0(colnames(model_pred), "_pred_")
	#model_pred <- scales::rescale(x=as.matrix(model_pred), to=c(0,1))
	pred_true <- data.frame(dummies::dummy(model$pred$obs))
	colnames(pred_true) <- paste0(sel_levels, "_true")
	model_pred <- cbind(pred_true, model_pred)

	# Create classification objects
	sel_obs <- model$pred$obs
	sel_pred <- model$pred$pred
	sel_prob <- as.numeric(as.factor(sel_obs)); for (i in sel_levels) sel_prob[which(sel_obs==i)] <- as.numeric(model$pred[which(model$pred$obs==i), which(colnames(model$pred)==i)])

	# Multi-class AUC as defined by Hand & Till (2001)
	if (nlevels(sel_factor) > 2) {
		multiclass_roc <- pROC::multiclass.roc(response=sel_obs, predictor=sel_prob, levels=sel_levels, percent=FALSE, print.auc=TRUE)
		multiclass_auc <- round(as.numeric(multiclass_roc$auc),3)
	}

	# Multi-class classification metrics
	sel_pred <- as.factor(sel_pred)
	#levels(sel_pred) <- levels(sel_obs)
	model_metrics <- mltest::ml_test(predicted=sel_pred, true=sel_obs)

	# Multi-class Classification rate (= 1 - error rate)
	mcr <- f.classification_rate(sel_factor=sel_obs, predicted=sel_pred)

	# Balanced Error Rate (BER) (= 1 - multi-class classification rate)
	ber <- mlr::measureBER(truth=sel_obs, response=sel_pred)

	# Accuracy
	accuracy <- f.accuracy(sel_factor=sel_obs, predicted=sel_pred)
	weighted_accuracy <- f.weighted_accuracy(sel_factor=as.factor(sel_obs), predicted=as.factor(sel_pred))

	# R2
	R2 <- f.r2(actual=as.numeric(as.factor(sel_obs)), predicted=as.numeric(as.factor(sel_pred)))

	# ROC curve (one vs. rest)
	model_roc <- plot.roc(as.numeric(as.factor(sel_obs)), as.numeric(as.factor(sel_pred)), main="Confidence intervals", percent=TRUE, ci=TRUE, print.auc=TRUE)
	model_roc_confidence <- ci.se(model_roc, specificities=seq(0, 100, 1), of="thresholds", thresholds="local maximas") # Confidence
	model_roc_sensitivity <- ci.se(model_roc, specificities=seq(0, 100, 10), of="thresholds", thresholds="local maximas") # Sensitivity
	model_roc_specificity <- ci.sp(model_roc, sensitivities=seq(0, 100, 10), of="thresholds", thresholds="local maximas") # Specificity

	# PR curve
	#model_pr <- pr.curve(scores.class0=sel_prob[sel_obs == sel_pred], scores.class1=sel_prob[sel_obs != sel_pred], curve=TRUE)

	# Save metrics
	sel[["_model_roc_"]] <- model_roc
	sel[["_model_roc_confidence_"]] <- model_roc_confidence
	sel[["_model_roc_sensitivity_"]] <- model_roc_sensitivity
	sel[["_model_roc_specificity_"]] <- model_roc_specificity
	sel[["_multiclass_metrics_"]] <- model_metrics
	if (nlevels(sel_factor) > 2) { sel[["_model_multiclass_roc_"]] <- multiclass_roc }
	if (nlevels(sel_factor) > 2) { sel[["_multiclass_auc_"]] <- multiclass_auc }
	sel[["_multi_class_rate_"]] <- mcr
	sel[["_balanced_error_rate_"]] <- ber
	sel[["_accuracy_"]] <- accuracy
	sel[["_weighted_accuracy_"]] <- weighted_accuracy
	sel[["_model_r2_"]] <- R2
	sel[["_AUC_"]] <- round(as.numeric(model_roc$auc)/100, 3)
	#sel[["_AUCPR_"]] <- round(model_pr$auc.integral, 3)

	# Calculate ROC and PR with multiROC
	model_roc <- multi_roc(model_pred, force_diag=TRUE)
	model_pr <- multi_pr(model_pred, force_diag=TRUE)

	# Plot ROC using multiROC
	plot_roc <- plot_roc_data(model_roc)
	plot(x=c(0,1), y=c(0,1), type="l", xlim=c(0,1), ylim=c(0,1), xlab="1 - Specificity", ylab="Sensitivity", main="ROC")
	for (i in sel_levels) {
		lines(x=1-plot_roc[which(plot_roc$Group==i),"Specificity"], y=plot_roc[which(plot_roc$Group==i),"Sensitivity"], col=sel_colors[which(sel_levels==i)], lwd=2)
	}
	legend("bottomright", legend=paste0(model_roc$Groups, " (AUC=", round(unlist(model_roc$AUC)[model_roc$Groups]*100, 2), "%)"), col=sel_colors, lwd=2, cex=0.75)

	# Plot ROC using pROC
	model_proc <- list()
	for (i in c(1:length(sel_levels))) {
		if (i==1) add=FALSE else add=TRUE
		model_proc[[as.character(sel_levels[i])]] <- plot.roc(model_pred[,i], model_pred[,i+length(sel_levels)], main="ROC curves of levels", col=sel_colors[i], percent=TRUE, ci=TRUE, print.auc=FALSE, add=add)
		plot(ci.se(model_proc[[as.character(sel_levels[i])]], specificities=seq(0, 100, 1), of="thresholds", thresholds="local maximas"), type="shape", col=rgb(col2rgb(sel_colors[i])[1]/255, col2rgb(sel_colors[i])[2]/255, col2rgb(sel_colors[i])[3]/255, alpha=0.3), lty=0, no.roc=TRUE)
	}
	legend("bottomright", legend=paste0(sel_levels, ", AUC: ", round(as.numeric(unlist(lapply(model_proc, function(x) x$auc))), 1), "%"), col=sel_colors, lwd=2, cex=0.75)

	# Plot Precision and Recall Curve using multiROC
	plot_pr <- plot_pr_data(model_pr)
	plot(x=c(0,1), y=c(1,0), type="n", xlim=c(0,1), ylim=c(0,1), xlab="Recall", ylab="Precision", main="PR")
	for (i in sel_levels) {
		lines(x=plot_pr[which(plot_pr$Group==i),"Recall"], y=plot_pr[which(plot_pr$Group==i),"Precision"], col=sel_colors[which(sel_levels==i)], lwd=2)
	}
	legend("bottomleft", legend=paste0(model_pr$Groups, " (AUC-PR=", round(unlist(model_pr$AUC)[model_pr$Groups]*100, 2), "%)"), col=sel_colors, lwd=2, cex=0.75)

	# Plot Precision and Recall Curve using pROC
	model_ppr <- list()
	plot(xlim=c(0,1), ylim=c(0,1), x=NULL, y=NULL, xlab="Recall", ylab="Precision", main="Precision Recall Curves of levels")
	for (i in 1:length(sel_levels)) {
		model_ppr[[as.character(sel_levels[i])]] <- pr.curve(scores.class0=model_pred[as.logical(model_pred[,i]),i+length(sel_levels)], scores.class1=model_pred[!as.logical(model_pred[,i]),i+length(sel_levels)], curve=TRUE)
		lines(x=model_ppr[[as.character(sel_levels[i])]]$curve[,1], y=model_ppr[[as.character(sel_levels[i])]]$curve[,2], xlab="Recall",ylab="Precision", t="l", col=sel_colors[i], lwd=2)
	}
	legend("bottomleft", legend=paste0(sel_levels, ", AUC-PR: ", round(as.numeric(unlist(lapply(model_ppr, function(x) x$auc.integral*100))), 1), "%"), col=sel_colors, lwd=2, cex=0.75)

	# Save metrics
	sel[["_AUCPR_"]] <- model_ppr[[1]]$auc.integral
	sel[["_model_pr_"]] <- model_pr
	sel[["_model_factor_roc_"]] <- model_proc
	sel[["_model_factor_pr_"]] <- model_ppr

	return(sel)
}



# ---------- Select features from model ----------
f.select_features_from_model <- function(feat_list, feat_class, model_varimp, keepx_min, confidence=0.95) {
	if (ncol(as.data.frame(model_varimp)) < 2 ) {
		model_varimp <- as.data.frame(model_varimp)
		model_varimp <- cbind(model_varimp, model_varimp[,1])
		colnames(model_varimp) <- as.character(unique(feat_class))
	}

	# Create list
	sel_list <- list()

	# Limit selected features between keepx_max and keepx_min
	for (i in unique(feat_class)) {
		elements <- colnames(feat_list)[which(model_varimp[,i] >= confidence * max(model_varimp[,i]))]

		#if (length(elements) > keepx_max) {
		#  sel_list[[i]] <- colnames(feat_list)[base::order(model_varimp[, i], decreasing=TRUE)[1:keepx_max]]
		#} else
		if (length(elements) < keepx_min) {
			sel_list[[i]] <- colnames(feat_list)[base::order(model_varimp[, i], decreasing=TRUE)[1:keepx_min]]
		} else {
			sel_list[[i]] <- elements
		}

		sel_list[[i]] <- sort(sel_list[[i]])
	}

	# Return selected features
	return(sel_list)
}



# ---------- Draw heatmap of selected features ----------
f.count.selected_features <- function(sel_feat) {
	# Return number of selected features
	return(length(unique(unlist(sel_feat))))
}



# ---------- PLS ----------
#' Function for variables selection using Partial Least Squares Discriminant Analysis
#'
#' @param feat_matrix the matrix containing the molecular information, with samples in rows and molecular features in columns
#' @param sel_factor the factor which is used to select the variables
#' @param sel_colors color vector for the factor
#' @param components number of components, defaults to 2
#' @param tune_length how often repeated cross validation is performed, defaults to 10
#' @param quantile_threshold confidence level, defaults to 0.95
#' @param plot_roc_filename plot metrics to external pdf, defaults to NULL
#' @export
#' @import caret multiROC PRROC pROC
#_' @examples
#_' select_features_pls(feat_matrix=tenerife$comp_list, sel_factor=tenerife$species, sel_colors=tenerife$colors, components=5)
select_features_pls <- function(feat_matrix, sel_factor, sel_colors, components=2, tune_length=10, quantile_threshold=0.95, plot_roc_filename=NULL) {
	# Make factors readible by R
	sel_factor <- as.factor(make.names(sel_factor))

	# Handle plotting
	if (is.null(plot_roc_filename)) {
		pdf(file=NULL)
	} else {
		pdf(file=as.character(plot_roc_filename), encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
	}

	# Detach conflicting mixOmics package
	#if ("package:mixOmics" %in% search()) detach(package:mixOmics, unload=TRUE)

	# Train PLS model
	#tuneGrid=data.frame(ncomp=2)
	model_pls <- caret::train(x=as.matrix(feat_matrix), y=sel_factor, method="pls",
							  preProcess=c("center", "scale"),
							  tuneGrid=data.frame(ncomp=components),
							  tuneLength=tune_length, trControl=caret::trainControl(method="repeatedcv", number=10, repeats=5, classProbs=TRUE, savePredictions="final"))
	print(paste("Number of chosen components:",as.numeric(model_pls$bestTune)))

	# Get variable importances
	imp_pls <- varImp(object=model_pls)
	rownames(imp_pls$importance) <- as.character(rownames(imp_pls$importance))

	# Names of selected features
	sel_pls <- f.select_features_from_model(feat_list=feat_matrix, feat_class=sel_factor, model_varimp=imp_pls$importance, confidence=quantile_threshold, keepx_min=10)

	# Save selected variables
	sel_pls[["_selected_variables_"]] <- unique(unlist(f.select_features_from_model(feat_list=feat_matrix, feat_class=sel_factor, model_varimp=imp_pls$importance, confidence=quantile_threshold, keepx_min=10)))

	# Performance measures
	sel_pls <- do.call(c, list(sel_pls, f.performance_measures_caret(model=model_pls, sel_factor=sel_factor, sel_colors=sel_colors)))

	dev.off()

	return(sel_pls)
}



# ---------- Random Forest ----------
#' Function for variables selection using Random Forest
#'
#' @param feat_matrix the matrix containing the molecular information, with samples in rows and molecular features in columns
#' @param sel_factor the factor which is used to select the variables
#' @param sel_colors color vector for the factor
#' @param tune_length how often repeated cross validation is performed, defaults to 10
#' @param quantile_threshold confidence level, defaults to 0.95
#' @param plot_roc_filename plot metrics to external pdf, defaults to NULL
#' @export
#' @import caret randomForest multiROC PRROC pROC
#_' @examples
#_' select_features_random_forest(feat_matrix=tenerife$comp_list, sel_factor=tenerife$species, sel_colors=tenerife$colors, components=5)
select_features_random_forest <- function(feat_matrix, sel_factor, sel_colors, tune_length=10, quantile_threshold=0.95, plot_roc_filename=NULL) {
	# Make factors readible by R
	sel_factor <- as.factor(make.names(sel_factor))

	# Handle plotting
	if (is.null(plot_roc_filename)) {
		pdf(file=NULL)
	} else {
		pdf(file=as.character(plot_roc_filename), encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
	}

	# Train RF model
	model_rf <- caret::train(x=as.matrix(feat_matrix), y=sel_factor, method="rf", importance=TRUE, proximity=TRUE,
							 tuneLength=tune_length, trControl=caret::trainControl(method="repeatedcv", number=10, repeats=5, classProbs=TRUE, savePredictions="final"))

	# Get variable importances
	imp_rf <- varImp(object=model_rf)
	rownames(imp_rf$importance) <- as.character(rownames(imp_rf$importance))

	# Save names of selected features
	sel_rf <- f.select_features_from_model(feat_list=feat_matrix, feat_class=sel_factor, model_varimp=imp_rf$importance, confidence=quantile_threshold, keepx_min=10)

	# Save selected variables
	sel_rf[["_selected_variables_"]] <- unique(unlist(f.select_features_from_model(feat_list=feat_matrix, feat_class=sel_factor, model_varimp=imp_rf$importance, confidence=quantile_threshold, keepx_min=10)))

	# Performance measures
	sel_rf <- do.call(c, list(sel_rf, f.performance_measures_caret(model=model_rf, sel_factor=sel_factor, sel_colors=sel_colors)))

	dev.off()

	return(sel_rf)
}



# ---------- Draw heatmap of selected variables ----------
#' Function to draw a heatmap of selected variables
#'
#' @param feat_list the matrix containing the molecular information, with samples in rows and molecular features in columns
#' @param sel_feat vector of selected variables
#' @param sel_names vector of names of selected variables
#' @param sample_colors color vector for the factor
#' @param filename plot to external pdf, use NULL for internal plotting
#' @param main title of the plot
#' @param scale scaling of columns or rows, defaults to "col"
#' @param plot_width width of plot, defaults to 6
#' @param plot_height height of plot, defaults to 5
#' @param cex_col size for column, defaults to 0.5
#' @param cex_row size for row, defaults to 0.7
#' @export
#' @import gplots
#_' @examples
#_' heatmap.selected_features(feat_list=tenerife$comp_list, sel_feat=sel_pls_subclass_int_list$`_selected_variables_`, sample_colors=mzml_pheno_colors_samples, plot_width=7, plot_height=7, cex_col=0.5, cex_row=0.4, filename="plots_species/ms2_subclass_int_list_select_pls.pdf", main="PLS")
heatmap.selected_features <- function(feat_list, sel_feat, sel_names=NULL, sample_colors=NULL, filename, main, scale="col", plot_width=6, plot_height=5, cex_col=0.5, cex_row=0.7) {
	# Use existing dendrogram for rows
	if (any(names(sel_feat) %in% '_dendrogram_row_')) {
		print("Using existing dendrogram for clustering of rows.")
		rowv = sel_feat$'_dendrogram_row_'
		sel_feat$'_dendrogram_row_' <- NULL
	} else {
		rowv = NULL
	}

	# Use existing dendrogram for columns
	if (any(names(sel_feat) %in% '_dendrogram_col_')) {
		print("Using existing dendrogram for clustering of columns.")
		colv = sel_feat$'_dendrogram_col_'
		sel_feat$'_dendrogram_col_' <- NULL
	} else {
		colv = NULL
	}

	# Plot '_selected_variables_' if they exist
	if (any(names(sel_feat) %in% '_selected_variables_')) {
		print("Using matrix with selected variables.")
		sel_list <- scale((feat_list[, which(colnames(feat_list) %in% sel_feat[["_selected_variables_"]])]),scale=T,center=T)
	} else {
		# Data frame with only selected features
		sel_list <- as.data.frame(feat_list[,c(sort(as.character(unique(unlist(sel_feat))))), drop=FALSE])
	}

	# Use sel_names
	if (! is.null(sel_names)) {
		colnames(sel_list) <- sel_names
	}

	# Use colors for samples (row)
	if (is.null(sample_colors)) {
		sample_colors <- "black"
	}

	# Clustering of rows and columns
	if (is.null(rowv)) {
		rowv = as.dendrogram(hclust(dist(scale(sel_list),method="euclidean"),method="complete"), center=T)
		#rowv = as.dendrogram(hclust(dist(as.numeric(as.factor(mymetadata$SpecCode))),method="ward.D"), center=T), offsetRow=0, #<< samples need to be sorted by SpecCode (not clustered)
	}
	if (is.null(colv)) {
		colv = as.dendrogram(hclust(dist(t(scale(sel_list)),method="euclidean"),method="ward.D"), center=T)
	}

	# Draw heatmap
	if (! is.null(filename)) pdf(file=as.character(filename), encoding="ISOLatin1", pointsize=10, width=plot_width, height=plot_height, family="Helvetica")
	heatmap.2(x=as.matrix(sel_list), scale=scale, cexRow=cex_row, cexCol=cex_col, main=main,
			  Rowv=rowv, offsetRow=0, colRow=sample_colors,
			  Colv=colv, offsetCol=0,
			  col=colorRampPalette(c('darkblue','white','darkred'), alpha=0.1, bias=1)(256),
			  #trace="none", margins=c(2,4),
			  trace="none", margins=c(max(nchar(colnames(sel_list)))/5, max(nchar(rownames(sel_list)))/3),
			  key=TRUE, key.title="Color key", density.info='density', denscol="black")
	if (! is.null(filename)) dev.off()
}



# ---------- Sunburst plot ----------
#' Function to draw a sunburst plot
#'
#' @param classifierClasses compound class names
#' @param numberOfSpectra frequency of classes
#' @param colorStart color start, defaults to 0.0
#' @param colorAlpha color alpha, defaults to 0.6
#' @export
#' @import circlize plotrix
#_' @examples
#_' sunBurstPlotFromSubstanceClasses(rownames(tenerife$div_classes), tenerife$div_classes$frequency, colorStart=0.0, colorAlpha=0.6)
sunBurstPlotFromSubstanceClasses <- function(classifierClasses, numberOfSpectra, colorStart=0.0, colorAlpha=0.6){
	level <- unlist(lapply(X = strsplit(x = classifierClasses, split = "; "), FUN = length))

	## all class levels
	classesAndSubClasses <- lapply(X = strsplit(x = classifierClasses, split = "; "), FUN = function(x){
		sapply(X = seq_along(x), FUN = function(y){paste(x[1:y], collapse = "; ")})
	})
	classesByLevel <- list()
	labelsByLevel <- list()
	for(levelHere in seq_len(max(level))){
		classesByLevel[[levelHere]] <- sort(unique(unlist(lapply(X = classesAndSubClasses, FUN = function(y){
			if(length(y) < levelHere) return(NULL)
			else return(y[[levelHere]])
		}))))
		labelsByLevel[[levelHere]] <- unlist(lapply(X = strsplit(x = classesByLevel[[levelHere]], split = "; "), FUN = tail, n=1))
	}

	## class counts
	countsByLevel <- list()
	for(levelHere in rev(seq_len(max(level)))){
		countsByLevel[[levelHere]] <- unlist(lapply(X = classesByLevel[[levelHere]], FUN = function(class){
			newSpectra <- ifelse(test = class %in% classifierClasses, yes = numberOfSpectra[[which(class == classifierClasses)]], no = 0)
			oldSpectra <- ifelse(test = levelHere < max(level), yes = sum(countsByLevel[[levelHere+1]][grepl(x = classesByLevel[[levelHere+1]], pattern = paste("^", class, sep = ""))]), no = 0)
			return(newSpectra + oldSpectra)
		}))
	}
	rootCount <- sum(countsByLevel[[1]])

	## coordinates
	colors <- rainbow(start = colorStart, alpha = colorAlpha, n = 1000)
	startDegreeByLevel <- list()
	spanDegreeByLevel <- list()
	colorByLevel <- list()

	for(levelHere in seq_len(max(level))){
		startDegreeByLevel[[levelHere]] <- list()
		spanDegreeByLevel[[levelHere]] <- list()
		colorByLevel[[levelHere]] <- list()

		classesToProcess <- classesByLevel[[levelHere]]
		precursorClasses <- NULL
		if(levelHere == 1)  precursorClasses <- ""
		else                precursorClasses <- classesByLevel[[levelHere-1]]

		for(precursorClassIdx in seq_along(precursorClasses)){
			precursorClass <- precursorClasses[[precursorClassIdx]]
			classesToProcessHereSelection <- grepl(x = classesToProcess, pattern = precursorClass)
			classesToProcessHere <- classesToProcess[classesToProcessHereSelection]
			startDegree <- ifelse(test = levelHere == 1, yes = 0, no = startDegreeByLevel[[levelHere-1]][[precursorClassIdx]])
			scalingFactor <- ifelse(test = levelHere == 1, yes = 1, no = countsByLevel[[levelHere-1]][[precursorClassIdx]] / sum(countsByLevel[[levelHere]][classesToProcessHereSelection])) ## ambiguous classes
			#startColor  <- ifelse(test = levelHere == 1, yes = 0, no = startDegreeByLevel[[levelHere-1]][[precursorClassIdx]])
			for(classToProcessHere in classesToProcessHere){
				classIdx <- which(classesByLevel[[levelHere]] == classToProcessHere)
				degreeSpan <- 360 * countsByLevel[[levelHere]][[classIdx]] / rootCount * scalingFactor
				startDegreeByLevel[[levelHere]][[classIdx]] <- startDegree
				spanDegreeByLevel [[levelHere]][[classIdx]] <- degreeSpan
				colorByLevel      [[levelHere]][[classIdx]] <- colors[[(floor(startDegree + degreeSpan / 2) / 360 * length(colors)) + ifelse(test = (floor(startDegree + degreeSpan / 2) / 360 * length(colors))==length(colors), yes = 0, no = 1) ]]
				startDegree <- startDegree + degreeSpan
			}
		}
	}

	thereIsNextLevelByLevel <- list()
	for(levelHere in seq_len(max(level))){
		thereIsNextLevelByLevel[[levelHere]] <- list()
		if(levelHere == max(level)){
			thereIsNextLevelByLevel[[levelHere]] <- rep(x = FALSE, times = length(classesByLevel[[levelHere]]))
		} else {
			for(classIdx in seq_along(classesByLevel[[levelHere]]))
				thereIsNextLevelByLevel[[levelHere]][[classIdx]] <- any(grepl(x = classesByLevel[[levelHere+1]], pattern = classesByLevel[[levelHere]][[classIdx]]))
		}
	}

	degreeThresholdForDrawing <- 0.5
	#maxLevel <- max(which(unlist(lapply(X = spanDegreeByLevel, FUN = function(x){any(unlist(x) >= degreeThresholdForDrawing)}))))

	plotRadius <- 8
	plotCex1 <- 1
	plotCex2 <- 1

	plot(1, type="n", xlab="", ylab="", xlim=c(-plotRadius, plotRadius), ylim=c(-plotRadius, plotRadius), axes = FALSE)

	## circle segments
	tmp <- sapply(X = seq_along(classesByLevel), FUN = function(levelHere){
		sapply(X = seq_along(classesByLevel[[levelHere]]), FUN = function(classIdx){
			if(spanDegreeByLevel[[levelHere]][[classIdx]] < degreeThresholdForDrawing)  return()
			draw.sector(
				start.degree = startDegreeByLevel[[levelHere]][[classIdx]],
				end.degree = startDegreeByLevel[[levelHere]][[classIdx]] + spanDegreeByLevel[[levelHere]][[classIdx]],
				rou1 = levelHere - 1, rou2 = levelHere, center = c(0,0), clock.wise = FALSE, col = colorByLevel [[levelHere]][[classIdx]], border = "white"
			)
		})
	})
	## segment text
	minimumAngleToShowSegmentText <- 15
	tmp <- sapply(X = seq_along(classesByLevel), FUN = function(levelHere){
		sapply(X = seq_along(classesByLevel[[levelHere]]), FUN = function(classIdx){
			if(spanDegreeByLevel[[levelHere]][[classIdx]] < minimumAngleToShowSegmentText)  return()
			textTokens <- strwrap(x = labelsByLevel[[levelHere]][[classIdx]], width = max(nchar(strsplit(x = "Some text", split = " ")[[1]])))
			#firstOffset <- 1 / (length(textTokens) * 2)

			for(idx in seq_along(textTokens)){
				#offset <- firstOffset * (2 * idx - 1)
				middle <- (startDegreeByLevel[[levelHere]][[classIdx]] + spanDegreeByLevel[[levelHere]][[classIdx]]/2)/360 * 2*pi
				isSwitched <- middle > pi/2 & middle < 3 * pi/2
				isSwitched <- middle > pi
				offset <-  ifelse(test = middle > pi, yes = (length(textTokens) - idx + 1) / (length(textTokens) + 1), no = idx / (length(textTokens) + 1))
				#if(isSwitched) textTokens[[idx]] <- rev(textTokens[[idx]])
				#label <- ifelse(test = isSwitched, yes = paste(strsplit(x = labelsByLevel[[levelHere]][[classIdx]], split = " ")[[1]], collapse = " "), no = labelsByLevel[[levelHere]][[classIdx]])
				#middle <- ifelse(test = middle < pi, yes = middle + pi, no = middle)
				arctext(
					x = textTokens[[idx]],
					center = c(0, 0), radius = levelHere - offset - 0.04,
					middle = middle,
					cex = plotCex1, stretch = 1,
					clockwise = !isSwitched
				)
			}
		})
	})

	## outer text
	levelMaxHere <- max(level)
	tmp <- sapply(X = seq_along(classesByLevel), FUN = function(levelHere){
		if(levelHere > levelMaxHere)  return()
		sapply(X = seq_along(classesByLevel[[levelHere]]), FUN = function(classIdx){
			if(thereIsNextLevelByLevel[[levelHere]][[classIdx]] | spanDegreeByLevel[[levelHere]][[classIdx]] >= minimumAngleToShowSegmentText)  return()
			#radius <- maxLevel + 0.2
			radius <- levelHere + 0.2
			angle <- (startDegreeByLevel[[levelHere]][[classIdx]] + spanDegreeByLevel[[levelHere]][[classIdx]]/2)/360 * 2*pi
			x <- radius * sin(angle+pi/2)
			y <- radius * cos(angle+pi/2)
			srt <- startDegreeByLevel[[levelHere]][[classIdx]] + spanDegreeByLevel[[levelHere]][[classIdx]]/2
			isSwitched <- srt > 90 & srt < 270
			if(isSwitched) adj <- c(1,0) else adj <- c(0,0.5)
			srt <- ifelse(test = isSwitched, yes = srt + 180, no = srt)
			text(
				x = x, y = -y, labels = labelsByLevel[[levelHere]][[classIdx]], adj = adj,
				srt=srt,
				cex = plotCex2
			)
		})
	})
}


