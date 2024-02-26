## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
#devtools::install_github("https://github.com/ipb-halle/iESTIMATE")
library(iESTIMATE)

## ----fig.width=16, fig.height=16, out.width=800, out.height=800---------------
sunBurstPlotFromSubstanceClasses(rownames(marchantiales$div_classes), marchantiales$div_classes$frequency, colorStart=0.0, colorAlpha=0.6)

## ----fig.width=16, fig.height=16, out.width=800, out.height=800---------------
sunBurstPlotFromSubstanceClasses(rownames(marchantiales$div_npclasses), marchantiales$div_npclasses$frequency, colorStart=0.0, colorAlpha=0.6)

## -----------------------------------------------------------------------------
library(parallel)
library(e1071)
library(doMC)
nSlaves <- detectCores(all.tests=FALSE, logical=FALSE)
registerDoMC(nSlaves)

## ----message=FALSE, fig.width=8, fig.height=6, out.width=800, out.height=600----
suppressWarnings(
sel_pls_comp_list <- select_features_pls(feat_matrix=marchantiales$comp_list, sel_factor=as.factor(marchantiales$metadata$species), sel_colors=marchantiales$metadata$color, components=(nlevels(as.factor(marchantiales$metadata$species))-1))
)

## ----message=FALSE, fig.width=8, fig.height=6, out.width=800, out.height=600----
print(paste("Number of essential variables:", length(unique(unlist(sel_pls_comp_list$`_selected_variables_`)))))
print(sel_pls_comp_list$`_selected_variables_`)
print(sel_pls_comp_list$`_multiclass_metrics_`)
print(sel_pls_comp_list$`_model_r2_`)

## ----message=FALSE, fig.width=8, fig.height=6, out.width=800, out.height=600----
heatmap.selected_features(feat_list=marchantiales$comp_list, sel_feat=sel_pls_comp_list$`_selected_variables_`, sample_colors=marchantiales$colors, plot_width=10, plot_height=10, filename=NULL, main="PLS-DA")

## ----message=FALSE, fig.width=8, fig.height=6, out.width=800, out.height=600----
suppressWarnings(
sel_rf_comp_list <- select_features_random_forest(feat_matrix=marchantiales$comp_list, sel_factor=as.factor(marchantiales$metadata$species), sel_colors=marchantiales$metadata$color, tune_length=10, quantile_threshold=0.95, plot_roc_filename=NULL)
)

## ----message=FALSE, fig.width=8, fig.height=6, out.width=800, out.height=600----
print(paste("Number of essential variables:", length(unique(unlist(sel_rf_comp_list$`_selected_variables_`)))))
print(sel_rf_comp_list$`_selected_variables_`)
print(sel_rf_comp_list$`_multiclass_metrics_`)
print(sel_rf_comp_list$`_model_r2_`)

## ----message=FALSE, fig.width=8, fig.height=6, out.width=800, out.height=600----
heatmap.selected_features(feat_list=marchantiales$comp_list, sel_feat=sel_rf_comp_list$`_selected_variables_`, sample_colors=marchantiales$metadata$color, plot_width=8, plot_height=6, cex_col=0.1, cex_row=0.4, filename=NULL, main="Random Forest")

## ----message=FALSE, fig.width=8, fig.height=6, out.width=800, out.height=600----
library(heatmaply)
heatmaply(scale(marchantiales$comp_list[, which(colnames(marchantiales$comp_list) %in% sel_rf_comp_list$`_selected_variables_`)]), k_row=1, k_col=1, colors=colorRampPalette(c('darkblue','white','darkred'), alpha=0.1, bias=1)(256), file=NULL, selfcontained=TRUE, fontsize_row=6, fontsize_col=3)

## ----message=FALSE, fig.width=8, fig.height=6, out.width=800, out.height=600----
suppressWarnings(
sel_rf_class_list <- select_features_random_forest(feat_matrix=marchantiales$class_list, sel_factor=as.factor(marchantiales$metadata$species), sel_colors=marchantiales$metadata$color, tune_length=10, quantile_threshold=0.95, plot_roc_filename=NULL)
)

## ----message=FALSE, fig.width=8, fig.height=6, out.width=800, out.height=600----
print(paste("Number of essential variables:", length(unique(unlist(sel_rf_class_list$`_selected_variables_`)))))
print(sel_rf_class_list$`_selected_variables_`)
#print(sel_rf_class_list$`_multiclass_metrics_`)
print(sel_rf_class_list$`_model_r2_`)

## ----message=FALSE, fig.width=8, fig.height=6, out.width=800, out.height=600----
heatmaply(scale(marchantiales$class_list[, which(colnames(marchantiales$class_list) %in% sel_rf_class_list$`_selected_variables_`)]), k_row=1, k_col=1, colors=colorRampPalette(c('darkblue','white','darkred'), alpha=0.1, bias=1)(256), file=NULL, selfcontained=TRUE, fontsize_row=6, fontsize_col=3)

## ----message=FALSE, fig.width=8, fig.height=6, out.width=800, out.height=600----
suppressWarnings(
sel_rf_subclass_list <- select_features_random_forest(feat_matrix=marchantiales$subclass_list, sel_factor=as.factor(marchantiales$metadata$species), sel_colors=marchantiales$metadata$color, tune_length=10, quantile_threshold=0.95, plot_roc_filename=NULL)
)

## ----message=FALSE, fig.width=8, fig.height=6, out.width=800, out.height=600----
print(paste("Number of essential variables:", length(unique(unlist(sel_rf_subclass_list$`_selected_variables_`)))))
print(sel_rf_subclass_list$`_selected_variables_`)
#print(sel_rf_subclass_list$`_multiclass_metrics_`)
print(sel_rf_subclass_list$`_model_r2_`)

## ----message=FALSE, fig.width=8, fig.height=6, out.width=800, out.height=600----
heatmaply(scale(marchantiales$subclass_list[, which(colnames(marchantiales$subclass_list) %in% sel_rf_subclass_list$`_selected_variables_`)]), k_row=1, k_col=1, colors=colorRampPalette(c('darkblue','white','darkred'), alpha=0.1, bias=1)(256), file=NULL, selfcontained=TRUE, fontsize_row=6, fontsize_col=3)

## ----message=FALSE, fig.width=8, fig.height=6, out.width=800, out.height=600----
suppressWarnings(
sel_rf_superclass_list <- select_features_random_forest(feat_matrix=marchantiales$superclass_list, sel_factor=as.factor(marchantiales$metadata$species), sel_colors=marchantiales$metadata$color, tune_length=10, quantile_threshold=0.95, plot_roc_filename=NULL)
)

## ----message=FALSE, fig.width=8, fig.height=6, out.width=800, out.height=600----
print(paste("Number of essential variables:", length(unique(unlist(sel_rf_superclass_list$`_selected_variables_`)))))
print(sel_rf_superclass_list$`_selected_variables_`)
#print(sel_rf_superclass_list$`_multiclass_metrics_`)
print(sel_rf_superclass_list$`_model_r2_`)

## ----message=FALSE, fig.width=8, fig.height=6, out.width=800, out.height=600----
heatmaply(scale(marchantiales$superclass_list[, which(colnames(marchantiales$superclass_list) %in% sel_rf_superclass_list$`_selected_variables_`)]), k_row=1, k_col=1, colors=colorRampPalette(c('darkblue','white','darkred'), alpha=0.1, bias=1)(256), file=NULL, selfcontained=TRUE, fontsize_row=6, fontsize_col=3)

## ----message=FALSE, fig.width=8, fig.height=6, out.width=800, out.height=600----
suppressWarnings(
sel_rf_npclass_list <- select_features_random_forest(feat_matrix=marchantiales$npclass_list, sel_factor=as.factor(marchantiales$metadata$species), sel_colors=marchantiales$metadata$color, tune_length=10, quantile_threshold=0.95, plot_roc_filename=NULL)
)

## ----message=FALSE, fig.width=8, fig.height=6, out.width=800, out.height=600----
print(paste("Number of essential variables:", length(unique(unlist(sel_rf_npclass_list$`_selected_variables_`)))))
print(sel_rf_npclass_list$`_selected_variables_`)
#print(sel_rf_npclass_list$`_multiclass_metrics_`)
print(sel_rf_npclass_list$`_model_r2_`)

## ----message=FALSE, fig.width=8, fig.height=6, out.width=800, out.height=600----
heatmaply(scale(marchantiales$npclass_list[, which(colnames(marchantiales$npclass_list) %in% sel_rf_npclass_list$`_selected_variables_`)]), k_row=1, k_col=1, colors=colorRampPalette(c('darkblue','white','darkred'), alpha=0.1, bias=1)(256), file=NULL, selfcontained=TRUE, fontsize_row=6, fontsize_col=3)

## ----message=FALSE, fig.width=8, fig.height=6, out.width=800, out.height=600----
suppressWarnings(
sel_rf_nppathway_list <- select_features_random_forest(feat_matrix=marchantiales$nppathway_list, sel_factor=as.factor(marchantiales$metadata$species), sel_colors=marchantiales$metadata$color, tune_length=10, quantile_threshold=0.95, plot_roc_filename=NULL)
)

## ----message=FALSE, fig.width=8, fig.height=6, out.width=800, out.height=600----
print(paste("Number of essential variables:", length(unique(unlist(sel_rf_nppathway_list$`_selected_variables_`)))))
print(sel_rf_nppathway_list$`_selected_variables_`)
#print(sel_rf_nppathway_list$`_multiclass_metrics_`)
print(sel_rf_nppathway_list$`_model_r2_`)

## ----message=FALSE, fig.width=8, fig.height=6, out.width=800, out.height=600----
heatmaply(scale(marchantiales$nppathway_list[, which(colnames(marchantiales$nppathway_list) %in% sel_rf_nppathway_list$`_selected_variables_`)]), k_row=1, k_col=1, colors=colorRampPalette(c('darkblue','white','darkred'), alpha=0.1, bias=1)(256), file=NULL, selfcontained=TRUE, fontsize_row=6, fontsize_col=3)

## ----message=FALSE, fig.width=8, fig.height=6, out.width=800, out.height=600----
suppressWarnings(
sel_rf_mdes_list <- select_features_random_forest(feat_matrix=marchantiales$mdes_list, sel_factor=as.factor(marchantiales$metadata$species), sel_colors=marchantiales$metadata$color, tune_length=10, quantile_threshold=0.95, plot_roc_filename=NULL)
)

## ----message=FALSE, fig.width=8, fig.height=6, out.width=800, out.height=600----
print(paste("Number of essential variables:", length(unique(unlist(sel_rf_mdes_list$`_selected_variables_`)))))
print(sel_rf_mdes_list$`_selected_variables_`)
#print(sel_rf_mdes_list$`_multiclass_metrics_`)
print(sel_rf_mdes_list$`_model_r2_`)

## ----message=FALSE, fig.width=8, fig.height=6, out.width=800, out.height=600----
heatmaply(scale(marchantiales$mdes_list[, which(colnames(marchantiales$mdes_list) %in% sel_rf_mdes_list$`_selected_variables_`)]), k_row=1, k_col=1, colors=colorRampPalette(c('darkblue','white','darkred'), alpha=0.1, bias=1)(256), file=NULL, selfcontained=TRUE, fontsize_row=6, fontsize_col=3)

## -----------------------------------------------------------------------------
detach("package:doMC", unload=TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  f.export_maf(cbind(ms1_def_neg, t(feat_list_neg)), "data/metabolites_maf_neg.tsv")
#  f.annotate_maf_classes(maf_input="data/metabolites_maf_neg.tsv", maf_output="data/metabolites_maf_neg_classes.tsv")
#  f.annotate_maf_compounds(maf_input="data/metabolites_maf_neg_classes.tsv", maf_output="data/m_MTBLS2239_LC-MS_negative_reverse-phase_metabolite_profiling_v2_maf.tsv", polarity="neg", xcms_id=rownames(ms1_def_neg), pol_mode=rep("neg",nrow(ms1_def_neg)), smiles=ms1_def_neg$smiles, names=ms1_def_neg$name)

## ----eval=FALSE---------------------------------------------------------------
#  f.export_maf(cbind(ms1_def_pos, t(feat_list_pos)), "data/metabolites_maf_pos.tsv")
#  f.annotate_maf_classes(maf_input="data/metabolites_maf_pos.tsv", maf_output="data/metabolites_maf_pos_classes.tsv")
#  f.annotate_maf_compounds(maf_input="data/metabolites_maf_pos_classes.tsv", maf_output="data/m_MTBLS2239_LC-MS_positive_reverse-phase_metabolite_profiling_v2_maf.tsv", polarity="pos", xcms_id=rownames(ms1_def_pos), pol_mode=rep("pos",nrow(ms1_def_pos)), smiles=ms1_def_pos$smiles, names=ms1_def_pos$name)

