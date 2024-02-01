


# ---------- Shannon Diversity ----------
#' Function to calculate the Shannon diversity measure
#'
#' This function is now obsolete, as it returns the same values than vegan::diversity(x, index="shannon"). Since we do not expect variables with negative intensity, we exclude negative values and take the natural logarithm.
#'
#' @param p vector with response variables of one sample
#' @export
#' @examples
#' shannon.diversity(p=c(4,8))
shannon.diversity <- function(p) {
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
#' @importFrom vegan specnumber
#' @examples
#' menhinick.diversity(p=c(4,8))
menhinick.diversity <- function(p) {
	D_Mn <- length(p) / sqrt(specnumber(p))
}



# ---------- Tukey-Test ----------
#' Function to perform the Tukey post-hoc HSD test on response and model terms
#'
#' @param response vector with response variables
#' @param term vector with factorized terms
#' @export
#' @importFrom stats aov
#' @importFrom multcomp glht mcp cld
#' @examples
#' tukey.test(response=marchantiales$model_div$unique, term=as.factor(marchantiales$metadata$species))
tukey.test <- function(response, term) {
	model_anova <- aov(formula(response ~ term))
	model_mc <- glht(model_anova, mcp(term="Tukey"))
	model_cld <- cld(summary(model_mc), decreasing=TRUE)
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
#' @importFrom caret postResample
f.r2 <- function(actual, predicted) {
	actual <- as.numeric(actual)
	predicted <- as.numeric(predicted)

	R2 <- postResample(actual, predicted)[2]
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
	if (round(sum(weights),0) != 1) {
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
#' @importFrom scales rescale
#' @importFrom dummies dummy
#' @importFrom pROC multiclass.roc
#' @importFrom pROC plot.roc
#' @importFrom pROC ci.se
#' @importFrom pROC ci.sp
#' @importFrom mltest ml_test
#' @importFrom mlr measureBER
#' @importFrom PRROC pr.curve
#' @importFrom multiROC multi_roc
#' @importFrom multiROC multi_pr
#' @importFrom multiROC plot_roc_data
#' @importFrom multiROC plot_pr_data
#' @importFrom grDevices col2rgb pdf dev.off rgb
#' @importFrom graphics legend lines text
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
	#model_pred <- rescale(x=as.matrix(model_pred), to=c(0,1))
	pred_true <- data.frame(dummy(model$pred$obs))
	colnames(pred_true) <- paste0(sel_levels, "_true")
	model_pred <- cbind(pred_true, model_pred)

	# Create classification objects
	sel_obs <- model$pred$obs
	sel_pred <- model$pred$pred
	sel_prob <- as.numeric(as.factor(sel_obs)); for (i in sel_levels) sel_prob[which(sel_obs==i)] <- as.numeric(model$pred[which(model$pred$obs==i), which(colnames(model$pred)==i)])

	# Multi-class AUC as defined by Hand & Till (2001)
	if (nlevels(sel_factor) > 2) {
		multiclass_roc <- multiclass.roc(response=sel_obs, predictor=sel_prob, levels=sel_levels, percent=FALSE, print.auc=TRUE)
		multiclass_auc <- round(as.numeric(multiclass_roc$auc),3)
	}

	# Multi-class classification metrics
	sel_pred <- as.factor(sel_pred)
	#levels(sel_pred) <- levels(sel_obs)
	model_metrics <- ml_test(predicted=sel_pred, true=sel_obs)

	# Multi-class Classification rate (= 1 - error rate)
	mcr <- f.classification_rate(sel_factor=sel_obs, predicted=sel_pred)

	# Balanced Error Rate (BER) (= 1 - multi-class classification rate)
	ber <- measureBER(truth=sel_obs, response=sel_pred)

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
		class0 <- model_pred[as.logical(model_pred[,i]),i+length(sel_levels)]
		class0[is.na(class0)] <- 0
		class1 <- model_pred[!(as.logical(model_pred[,i])),i+length(sel_levels)]
		class1[is.na(class1)] <- 0
		model_ppr[[as.character(sel_levels[i])]] <- pr.curve(scores.class0=class0, scores.class1=class1, curve=TRUE)
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
#_' @importFrom base order
f.select_features_from_model <- function(feat_list, feat_class, model_varimp, keepx_min=10, keepx_max=0, confidence=0.95) {
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

		if (length(elements) < keepx_min) {
			select_elements <- keepx_min
		} else if ((length(elements) > keepx_max) & (keepx_max > 0)) {
			select_elements <- keepx_max
		} else {
			select_elements <- length(elements)
		}

		sel_list[[i]] <- colnames(feat_list)[order(model_varimp[, i], decreasing=TRUE)[1:select_elements]]

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
#' @param keepx_min Minimum number of variables to keep for level (defaults to 10)
#' @param keepx_max Maximum number of variables to keep for level (defaults to 0, which means no limit)
#' @param plot_roc_filename plot metrics to external pdf, defaults to NULL
#' @export
#' @importFrom pls mvr plsr cppls
#' @importFrom caret train
#' @importFrom caret trainControl
#' @importFrom caret varImp
#' @importFrom grDevices pdf dev.off
#' @examples
#' select_features_pls(feat_matrix=marchantiales$comp_list,
#' sel_factor=as.factor(marchantiales$metadata$species),
#' sel_colors=marchantiales$metadata$color,
#' components=(nlevels(as.factor(marchantiales$metadata$species))-1))
select_features_pls <- function(feat_matrix, sel_factor, sel_colors, components=2, tune_length=10, quantile_threshold=0.95, keepx_min=10, keepx_max=0, plot_roc_filename=NULL) {
	# Make factors readible by R
	sel_factor <- as.factor(make.names(sel_factor))

	# Handle plotting
	if (is.null(plot_roc_filename)) {
		pdf(file=NULL)
	} else {
		pdf(file=as.character(plot_roc_filename), encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
	}

	# Detach conflicting mixOmics package
	if ("package:mixOmics" %in% search()) detach(package:mixOmics, unload=TRUE)

	# Train PLS model
	#tuneGrid=data.frame(ncomp=2)
	model_pls <- train(x=as.matrix(feat_matrix), y=sel_factor, method="pls",
							  preProcess=c("center", "scale"),
							  tuneGrid=data.frame(ncomp=components),
							  tuneLength=tune_length, trControl=trainControl(method="repeatedcv", number=10, repeats=5, classProbs=TRUE, savePredictions="final"))
	print(paste("Number of chosen components:",as.numeric(model_pls$bestTune)))

	# Get variable importances
	imp_pls <- varImp(object=model_pls)
	rownames(imp_pls$importance) <- as.character(rownames(imp_pls$importance))

	# Names of selected features
	sel_pls <- f.select_features_from_model(feat_list=feat_matrix, feat_class=sel_factor, model_varimp=imp_pls$importance, confidence=quantile_threshold, keepx_min=keepx_min, keepx_max=keepx_max)

	# Save selected variables
	sel_pls[["_selected_variables_"]] <- unique(unlist(f.select_features_from_model(feat_list=feat_matrix, feat_class=sel_factor, model_varimp=imp_pls$importance, confidence=quantile_threshold, keepx_min=keepx_min, keepx_max=keepx_max)))
	sel_pls[["_varImp_"]] <- as.data.frame(imp_pls$importance)

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
#' @param keepx_min Minimum number of variables to keep for level (defaults to 10)
#' @param keepx_max Maximum number of variables to keep for level (defaults to 0, which means no limit)
#' @param plot_roc_filename plot metrics to external pdf, defaults to NULL
#' @export
#' @importFrom randomForest randomForest
#' @importFrom caret train
#' @importFrom caret trainControl
#' @importFrom caret varImp
#_' @examples
#_' select_features_random_forest(feat_matrix=marchantiales$comp_list,
#_' sel_factor=as.factor(marchantiales$metadata$species),
#_' sel_colors=marchantiales$metadata$color)
select_features_random_forest <- function(feat_matrix, sel_factor, sel_colors, tune_length=10, quantile_threshold=0.95, keepx_min=10, keepx_max=0, plot_roc_filename=NULL) {
	# Make factors readible by R
	sel_factor <- as.factor(make.names(sel_factor))

	# Handle plotting
	if (is.null(plot_roc_filename)) {
		pdf(file=NULL)
	} else {
		pdf(file=as.character(plot_roc_filename), encoding="ISOLatin1", pointsize=10, width=6, height=6, family="Helvetica")
	}

	# Train RF model (method="rf", method="ranger" is parallel)
	model_rf <- train(x=as.matrix(feat_matrix), y=sel_factor, method="rf", importance=TRUE, proximity=TRUE,
							 tuneLength=tune_length, trControl=trainControl(method="repeatedcv", number=10, repeats=5, classProbs=TRUE, savePredictions="final"))

	# Get variable importances
	imp_rf <- varImp(object=model_rf)
	rownames(imp_rf$importance) <- as.character(rownames(imp_rf$importance))

	# Save names of selected features
	sel_rf <- f.select_features_from_model(feat_list=feat_matrix, feat_class=sel_factor, model_varimp=imp_rf$importance, confidence=quantile_threshold, keepx_min=keepx_min, keepx_max=keepx_max)

	# Save selected variables
	sel_rf[["_selected_variables_"]] <- unique(unlist(f.select_features_from_model(feat_list=feat_matrix, feat_class=sel_factor, model_varimp=imp_rf$importance, confidence=quantile_threshold, keepx_min=keepx_min, keepx_max=keepx_max)))
	sel_rf[["_varImp_"]] <- as.data.frame(imp_rf$importance)

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
#_' @importFrom base scale
#_' @importFrom stats hclust
#_' @importFrom stats dist
#' @importFrom gplots heatmap.2
#' @importFrom grDevices colorRampPalette pdf dev.off
#' @importFrom stats as.dendrogram dist formula hclust
#' @examples
#' heatmap.selected_features(feat_list=marchantiales$comp_list,
#' sel_feat=select_features_pls(feat_matrix=marchantiales$comp_list,
#'                              sel_factor=as.factor(marchantiales$metadata$species),
#'                              sel_colors=marchantiales$metadata$color,
#'                              components=(nlevels(as.factor(marchantiales$metadata$species))-1))$`_selected_variables_`,
#' sample_colors=marchantiales$metadata$color, filename=NULL)
heatmap.selected_features <- function(feat_list, sel_feat, sel_names=NULL, sample_colors=NULL, filename=NULL, main="", scale="col", plot_width=6, plot_height=5, cex_col=0.5, cex_row=0.7) {
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
#' @importFrom grDevices colorRampPalette pdf dev.off rainbow
#' @importFrom graphics text
#' @importFrom utils tail
#' @importFrom circlize draw.sector
#' @importFrom plotrix arctext
#' @examples
#' sunBurstPlotFromSubstanceClasses(rownames(marchantiales$div_classes),
#' marchantiales$div_classes$frequency,
#' colorStart=0.0, colorAlpha=0.6)
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



# ---------- Make Classes at Level ----------
#' Function to extract information at a specific ontology level or below of classification information.
#'
#' @param level The level of the ontology to extract information to
#' @param div_samples Data frame with classes in rows and samples in columns.
#' @export
make_classes_at_level <- function(level=2, div_samples) {
	# Extract classes at specific level
	classes_samples_names <- NULL
	classes_samples_names <- c(classes_samples_names, lapply(X=strsplit(rownames(div_samples), '; '), FUN=function(x) { gsub(x=paste(x[1:level],sep='',collapse='; '),pattern='; NA',replacement='') }))

	# Count classes at specific level
	classes_samples <- data.frame()
	for (i in c(1:ncol(div_samples))) classes_samples <- rbind(classes_samples, rep(0, length(unique(classes_samples_names))))
	classes_samples <- t(classes_samples)
	colnames(classes_samples) <- colnames(div_samples)
	rownames(classes_samples) <- unique(classes_samples_names)
	for (i in rownames(classes_samples)) {
		for (j in c(1:ncol(div_samples))) {
			classes_samples[rownames(classes_samples)==i, j] <- sum(div_samples[grep(x=rownames(div_samples), pattern=i), j])
		}
	}

	# Create data frame
	classes <- classes_samples
	classes[is.na(classes)] <- 0
	classes <- apply(X=classes, MARGIN=1, FUN=function(x) { sum(x) })
	classes <- data.frame(row.names=names(classes), frequency=as.numeric(classes))

	# Imputation of NA with zeros
	classes[is.na(classes)] <- 0
	classes_samples[is.na(classes_samples)] <- 0

	# Classification list for statistics
	class_list <- as.data.frame(t(classes_samples))
	class_list[is.na(class_list)] <- 0

	# Return
	return(class_list)
}



# ---------- Gap Weighting based on Thiele (1993) ----------
#' Function to perform gap weighting on continuous values.
#'
#' This function converts continuous quantitative traits to distinct character states useful for calculating phylogenies.
#' Here, we apply the gap weighting algorithm by Thiele (1993), https://doi.org/10.1006/clad.1993.1020 .
#'
#' @param morph_char A vector with character values
#' @param states Number of states, defaults to 8
#' @export
#' @examples
#' gap_weighting(marchantiales$char_list$thallus.width, states=8)
gap_weighting <- function(morph_char, states=8) {
	suppressWarnings({ morph_num <- as.numeric(morph_char) })

	# Gap weighting
	gap_weighted <- as.factor(as.character(unlist(lapply(X=as.list(morph_char), FUN=function(x) {
		x
		if ((is.na(x)) | (x == "?") | (x == "-") | (x <= 0) | (x == "1")) {
			x
		} else {
			x = as.numeric(x)
			x = round(((x - as.numeric(min(morph_num[morph_num>0], na.rm=T))) / (as.numeric(max(morph_num, na.rm=T)) - as.numeric(min(morph_num[morph_num>0], na.rm=T))) * (states - 1)), digits=0)
		}
	}))))

	# Character states [ 0, 1, .., 9, A, .., Z]
	gap_weighted <- as.character(gap_weighted)
	for (i in c(10:36)) {
		gap_weighted[(gap_weighted==i)] <- LETTERS[(i-10+1)]
	}

	return(gap_weighted)
}



# ---------- Export peak list as MAF ----------
#' Function to export a peak table to MAF for use in MetaboLights.
#'
#' This function exports a peak table to MAF without any annotations. It should be run first.
#'
#' @param peak_list The peak table with MS1 features in rows and information and samples in columns
#' @param maf_filename The filename of the MAF file
#' @export
#' @importFrom utils write.table
export_maf <- function(peak_list, maf_filename) {
	# Preparations
	l <- nrow(peak_list)

	# These columns are defined by MetaboLights mzTab
	maf <- apply(X=data.frame(database_identifier=character(l),
							  chemical_formula=character(l),
							  smiles=character(l),
							  inchi=character(l),
							  metabolite_identification=character(l),
							  mass_to_charge=peak_list$mzmed,
							  fragmentation=character(l),
							  modifications=character(l),
							  charge=character(l),
							  retention_time=peak_list$rtmed,
							  taxid=character(l),
							  species=character(l),
							  database=character(l),
							  database_version=character(l),
							  reliability=character(l),
							  uri=character(l),
							  search_engine=character(l),
							  search_engine_score=character(l),
							  smallmolecule_abundance_sub=character(l),
							  smallmolecule_abundance_stdev_sub=character(l),
							  smallmolecule_abundance_std_error_sub=character(l),
							  xcms_identifier=rownames(peak_list),
							  peak_list,
							  stringsAsFactors=FALSE),
				 MARGIN=2, FUN=as.character)

	# Export MAF
	write.table(maf, file=maf_filename, row.names=FALSE, col.names=colnames(maf), quote=TRUE, sep="\t", na="\"\"")
}



# ---------- Update existing MAF with annotated compounds ----------
#' Function to update an unannotated MAF with classification information for use in MetaboLights.
#'
#' This function updates an unannotated MAF with classification information. It should be run second on an unannotated MAF.
#'
#' @param maf_input The unannotated MAF file
#' @param maf_output The filename of the resulting MAF
#' @export
#' @importFrom utils read.table write.table
#' @importFrom ontologyIndex get_ontology
annotate_maf_classes <- function(maf_input, maf_output) {
	# Read CHEMONT ontology
	obo <- get_ontology(file=url("https://raw.githubusercontent.com/stuchalk/cjld/master/ChemOnt_2_1.obo"), extract_tags="minimal")

	# Import MAF
	maf_out <- read.table(file=maf_input, quote="\"", sep="\t", na.strings="NA", header=TRUE, stringsAsFactors=TRUE)
	maf_out[is.na(maf_out)] <- as.character("")

	# Annotate classes
	for (i in c(1:length(maf_out$xcms_identifier))) {
		cl <- gsub(x=maf_out$primary_class[i], pattern='.*; ', replacement='')
		id <- as.character(obo$id[which(as.character(obo$name) %in% as.character(cl))])

		if (nchar(cl) > 1) {
			maf_out$database_identifier[i] <- paste(id, collapse='|')
			maf_out$metabolite_identification[i] <- paste(cl, collapse='|')
		}
	}

	# Export MAF
	write.table(maf_out, file=maf_output, row.names=FALSE, col.names=colnames(maf_out), quote=TRUE, sep="\t", na="\"\"")
}



# ---------- Update existing MAF with annotated compounds ----------
#' Function to annotate a MAF for use in MetaboLights.
#'
#' This function annotates a MAF with structure identification. It should be run second on an unannotated MAF or third on a MAF with classification information.
#'
#' @param maf_input The input MAF file
#' @param maf_output The filename of the resulting MAF
#' @param polarity The polarity used ("neg" or "pos")
#' @param xcms_id The id of XCMS peak detection, usually the FTxxxxx name
#' @param pol_mode The polarization mode for any feature (e.g. rep("neg",nrow(ms1_def_neg))
#' @param smiles Vector containing the SMILES information
#' @param names The UPAC names of the annotated compounds
#' @export
#' @importFrom utils read.table write.table
annotate_maf_compounds <- function(maf_input, maf_output, polarity, xcms_id, pol_mode, smiles, names) {
	# Import MAF
	maf_out <- read.table(file=maf_input, quote="\"", sep="\t", na.strings="NA", header=TRUE, stringsAsFactors=TRUE)
	maf_out[is.na(maf_out)] <- as.character("")

	maf_out$database_identifier <- as.character(maf_out$database_identifier)
	maf_out$metabolite_identification <- as.character(maf_out$metabolite_identification)

	# Annotate compounds
	for (i in c(1:length(xcms_id))) {
		if (pol_mode[i] == polarity) {
			j <- which(maf_out$xcms_identifier==xcms_id[i])
			if (length(j) > 0) {
				if (nchar(as.character(maf_out$database_identifier[j])) > 0) {
					maf_out$database_identifier[j] <- paste(c(smiles[i], maf_out$database_identifier[j]), collapse='|')
					maf_out$metabolite_identification[j] <- paste(c(names[i], maf_out$metabolite_identification[j]), collapse='|')
				} else {
					maf_out$database_identifier[j] <- paste(smiles[i], collapse='|')
					maf_out$metabolite_identification[j] <- paste(names[i], collapse='|')
				}
			}
		}
	}

	# Remove non-standard columns "ms_level", "primary_class"
	maf_out$ms_level <- NULL
	maf_out$primary_class <- NULL

	# Export MAF
	write.table(maf_out, file=maf_output, row.names=FALSE, col.names=colnames(maf_out), quote=TRUE, sep="\t", na="\"\"")
}



# ---------- Make high-level descriptors from CDK names ----------
#' Function to make abstracted high-level molecular descriptors from CDK names.
#'
#' @param descriptor_names A vector containing the names of the molecular descriptors
#' @export
make.high_level.descriptor <- function(descriptor_names) {
	descriptor_names <- gsub(x=descriptor_names, pattern="MolWeight", replacement="Molecular Mass")
	descriptor_names <- gsub(x=descriptor_names, pattern="nAtoms", replacement="Atom Count Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="numC", replacement="Number of C atoms")
	descriptor_names <- gsub(x=descriptor_names, pattern="numN", replacement="Number of N atoms")
	descriptor_names <- gsub(x=descriptor_names, pattern="numP", replacement="Number of P atoms")
	descriptor_names <- gsub(x=descriptor_names, pattern="numO", replacement="Number of O atoms")
	descriptor_names <- gsub(x=descriptor_names, pattern="numHydrogen", replacement="Number of Hydrogen atoms")
	descriptor_names <- gsub(x=descriptor_names, pattern="CNratio", replacement="Molecular C to N ratio")
	descriptor_names <- gsub(x=descriptor_names, pattern="XLogP", replacement="Partition Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="Fsp3", replacement="Molecular Non-flatness Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="MW", replacement="Molecular Weight")
	descriptor_names <- gsub(x=descriptor_names, pattern="LipinskiFailures", replacement="Drug Likeness")
	descriptor_names <- gsub(x=descriptor_names, pattern="nRotB", replacement="Molecular Bond Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="MLogP", replacement="Partition Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="nAtom.*", replacement="Atom Count Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="nB$", replacement="Molecular Bond Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="nBase", replacement="Basic Group Count Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="nAromBond.*", replacement="Aromaticity Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="naAromAtom.*", replacement="Aromaticity Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="(ALogP|ALogp2|AMR)", replacement="Partition Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="nAcid", replacement="Acidity Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="(nA|nR|nN|nD|nC|nF|nQ|nE|nG|nH|nI|nP|nL|nK|nM|nS|nT|nY|nV|nW)$", replacement="Amino Acid Count Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="nSmallRings", replacement="Molecular Ring Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="nAromRings", replacement="Molecular Ring Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="nRingBlocks", replacement="Molecular Ring Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="nAromBlocks", replacement="Molecular Ring Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="nRings\\d", replacement="Molecular Ring Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="tpsaEfficiency", replacement="Molecular Topological Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="Zagreb", replacement="Eccentric Connectivity Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="WPATH", replacement="Molecular Path Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="WPOL", replacement="Molecular Polarity Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="WTPT.*", replacement="Molecular Path Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="VAdjMat", replacement="Molecular Topological Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="VABC", replacement="Molecular Topological Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="TopoPSA", replacement="Molecular Topological Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="(topoShape|geomShape)", replacement="Molecular Topological Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="PetitjeanNumber", replacement="Molecular Topological Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="MDE(N|C|O)\\..*", replacement="Molecular Distance Edge Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="khs.*", replacement="Molecular Electron State Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="Kier\\d", replacement="Molecular Shape Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="HybRatio", replacement="Hybridization Ratio Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="fragC", replacement="Molecular Complexity Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="FMF", replacement="Molecular Complexity Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="ECCEN", replacement="Eccentric Connectivity Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="ATSc.*", replacement="Molecular Topological Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="ATSm.*", replacement="Molecular Topological Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="ATSp.*", replacement="Molecular Topological Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="(PPSA|PNSA|DPSA|FPSA|FNSA|WPSA|WNSA|RPC|RNC|THSA|TPSA|RHSA|RPSA).*", replacement="Molecular Topological Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="(SCH\\.|VCH\\.).*", replacement="Molecular Topological Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="(SC\\.|VC\\.).*", replacement="Molecular Topological Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="(SPC\\.|VPC\\.).*", replacement="Molecular Topological Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="(SP\\.|VP\\.).*", replacement="Molecular Topological Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="(C1SP1|C2SP1|C1SP2|C2SP2|C3SP2|C1SP3|C2SP3|C3SP3|C4SP3)", replacement="Carbon Connectivity Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="apol", replacement="Molecular Electronic Descriptor")
	descriptor_names <- gsub(x=descriptor_names, pattern="bpol", replacement="Molecular Electronic Descriptor")

	return(descriptor_names)
}


