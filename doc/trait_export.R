## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
#devtools::install_github("https://github.com/ipb-halle/iESTIMATE")
library(iESTIMATE)

## -----------------------------------------------------------------------------
try_table <- data.frame(Species=marchantiales$metadata$id_species[seq(1, nrow(marchantiales$metadata), 3)])
try_table$`Sample ID` <- marchantiales$metadata$id_samples[seq(1, nrow(marchantiales$metadata), 3)]

## -----------------------------------------------------------------------------
try_table$Latitude <- marchantiales$metadata$latitude[seq(1, nrow(marchantiales$metadata), 3)]
try_table$Longitude <- marchantiales$metadata$longitude[seq(1, nrow(marchantiales$metadata), 3)]
try_table$Altitude <- ""
try_table$Exposition <- "natural environment"
try_table$`Voucher specimen` <- marchantiales$metadata$voucher_specimen[seq(1, nrow(marchantiales$metadata), 3)]
try_table$`Plant maturity` <- "mature gametophyte"

## -----------------------------------------------------------------------------
try_table$`Sample ID` <- unique(gsub(x=try_table$`Sample ID`, pattern="\\.\\d.*", replacement=""))
try_table$`Sample ID` <- gsub(x=try_table$`Sample ID`, pattern="GOT|OEL", replacement="SWE")
try_table$`Sample ID` <- gsub(x=try_table$`Sample ID`, pattern="OE", replacement="SWE")
try_table$`Sample ID` <- gsub(x=try_table$`Sample ID`, pattern="HAL", replacement="GER")

## -----------------------------------------------------------------------------
try_table_colnames_1st <- c("Species", "Sample ID", "Latitude", "Longitude", "Altitude", "Exposition", "Voucher specimen", "Plant maturity")
try_table_colnames_2nd <- c("", "", "", "", "", "", "", "")

## -----------------------------------------------------------------------------
try_table_char_list <- marchantiales$char_list
try_table_char_list <- try_table_char_list[, -grep("(\\.mean)", colnames(try_table_char_list))]
colnames(try_table_char_list) <- gsub(x=colnames(try_table_char_list), pattern="\\.", replacement=" ")
colnames(try_table_char_list) <- paste(toupper(substr(colnames(try_table_char_list), 1, 1)), substr(colnames(try_table_char_list), 2, nchar(colnames(try_table_char_list))), sep="")
colnames(try_table_char_list) <- gsub(x=colnames(try_table_char_list), pattern=" cross", replacement=" cross-section")

## -----------------------------------------------------------------------------
rownames(try_table_char_list) <- gsub(x=rownames(try_table_char_list), pattern="Asterella.gracilis", replacement="A.gracilis.SWE")
rownames(try_table_char_list) <- gsub(x=rownames(try_table_char_list), pattern="Athalamia.hyalina.var..suecica", replacement="A.hyalina.SWE")
rownames(try_table_char_list) <- gsub(x=rownames(try_table_char_list), pattern="Mannia.fragrans", replacement="M.fragrans.SWE")
rownames(try_table_char_list) <- gsub(x=rownames(try_table_char_list), pattern="Reboulia.hemisphaerica.subsp..hemisphaerica", replacement="R.hemisphaerica.SWE")
rownames(try_table_char_list) <- gsub(x=rownames(try_table_char_list), pattern="Riccia.beyrichiana", replacement="R.beyrichiana.SWE")
rownames(try_table_char_list) <- gsub(x=rownames(try_table_char_list), pattern="Riccia.bifurca", replacement="R.bifurca.SWE")
rownames(try_table_char_list) <- gsub(x=rownames(try_table_char_list), pattern="Riccia.canaliculata", replacement="R.canaliculata.SWE")
rownames(try_table_char_list) <- gsub(x=rownames(try_table_char_list), pattern="Riccia.cavernosa", replacement="R.cavernosa.SWE")
rownames(try_table_char_list) <- gsub(x=rownames(try_table_char_list), pattern="Riccia.ciliifera.HAL", replacement="R.ciliifera.GER")
rownames(try_table_char_list) <- gsub(x=rownames(try_table_char_list), pattern="Riccia.ciliifera.SWE", replacement="R.ciliifera.SWE")
rownames(try_table_char_list) <- gsub(x=rownames(try_table_char_list), pattern="Riccia.gothica", replacement="R.gothica.SWE")
rownames(try_table_char_list) <- gsub(x=rownames(try_table_char_list), pattern="Riccia.gougetiana.var..armatissima", replacement="R.gougetiana.GER")
rownames(try_table_char_list) <- gsub(x=rownames(try_table_char_list), pattern="Riccia.huebeneriana", replacement="R.huebeneriana.SWE")
rownames(try_table_char_list) <- gsub(x=rownames(try_table_char_list), pattern="Riccia.sorocarpa", replacement="R.sorocarpa.SWE")
rownames(try_table_char_list) <- gsub(x=rownames(try_table_char_list), pattern="Riccia.subbifurca.SWE1", replacement="R.subbifurca.SWE1")
rownames(try_table_char_list) <- gsub(x=rownames(try_table_char_list), pattern="Riccia.subbifurca.SWE2", replacement="R.subbifurca.SWE2")

## -----------------------------------------------------------------------------
try_table <- cbind(try_table, try_table_char_list[match(try_table$`Sample ID`, rownames(try_table_char_list)), ])
try_table_colnames_1st <- c(try_table_colnames_1st, colnames(try_table_char_list))
try_table_colnames_2nd <- c(try_table_colnames_2nd, rep("", times=length(colnames(try_table_char_list))))

## -----------------------------------------------------------------------------
try_table_chemodiv_full <- cbind(marchantiales$model_div$richness, marchantiales$model_div$shannon, marchantiales$model_div$pielou, marchantiales$model_div$hillfunc)
colnames(try_table_chemodiv_full) <- c("Chemical richness", "Chemical Shannon diversity", "Chemical Pielou\\'s evenness index", "Chemical functional Hill diversity")
rownames(try_table_chemodiv_full) <- rownames(marchantiales$comp_list)

## -----------------------------------------------------------------------------
try_table_chemodiv <- NULL
for (i in unique(gsub(x=rownames(try_table_chemodiv_full), pattern="\\.\\d.*", replacement=""))) try_table_chemodiv <- rbind(try_table_chemodiv, apply(X=try_table_chemodiv_full[gsub(x=rownames(try_table_chemodiv_full), pattern="\\.\\d.*", replacement="")==i, ], MARGIN=2, FUN=function(x) { median(x) } ))
rownames(try_table_chemodiv) <- unique(gsub(x=rownames(try_table_chemodiv_full), pattern="\\.\\d.*", replacement=""))
rownames(try_table_chemodiv) <- gsub(x=rownames(try_table_chemodiv), pattern="GOT|OEL", replacement="SWE")
rownames(try_table_chemodiv) <- gsub(x=rownames(try_table_chemodiv), pattern="OE", replacement="SWE")
rownames(try_table_chemodiv) <- gsub(x=rownames(try_table_chemodiv), pattern="HAL", replacement="GER")

## -----------------------------------------------------------------------------
try_table <- cbind(try_table, try_table_chemodiv)
try_table_colnames_1st <- c(try_table_colnames_1st, colnames(try_table_chemodiv))
try_table_colnames_2nd <- c(try_table_colnames_2nd, rep("", times=length(colnames(try_table_chemodiv))))

## -----------------------------------------------------------------------------
try_table_nppathway_list_pos_full <- make_classes_at_level(level=1, div_samples=marchantiales$div_npclasses_samples_pos)

## -----------------------------------------------------------------------------
try_table_nppathway_list_pos <- NULL
for (i in unique(gsub(x=rownames(try_table_nppathway_list_pos_full), pattern="\\.\\d.*", replacement=""))) try_table_nppathway_list_pos <- rbind(try_table_nppathway_list_pos, apply(X=try_table_nppathway_list_pos_full[gsub(x=rownames(try_table_nppathway_list_pos_full), pattern="\\.\\d.*", replacement="")==i, ], MARGIN=2, FUN=function(x) { median(x) } ))
rownames(try_table_nppathway_list_pos) <- unique(gsub(x=rownames(try_table_nppathway_list_pos_full), pattern="\\.\\d.*", replacement=""))
rownames(try_table_nppathway_list_pos) <- gsub(x=rownames(try_table_nppathway_list_pos), pattern="GOT|OEL", replacement="SWE")
rownames(try_table_nppathway_list_pos) <- gsub(x=rownames(try_table_nppathway_list_pos), pattern="OE", replacement="SWE")
rownames(try_table_nppathway_list_pos) <- gsub(x=rownames(try_table_nppathway_list_pos), pattern="HAL", replacement="GER")

## -----------------------------------------------------------------------------
try_table_nppathway_list_neg_full <- make_classes_at_level(level=1, div_samples=marchantiales$div_npclasses_samples_neg)

## -----------------------------------------------------------------------------
try_table_nppathway_list_neg <- NULL
for (i in unique(gsub(x=rownames(try_table_nppathway_list_neg_full), pattern="\\.\\d.*", replacement=""))) try_table_nppathway_list_neg <- rbind(try_table_nppathway_list_neg, apply(X=try_table_nppathway_list_neg_full[gsub(x=rownames(try_table_nppathway_list_neg_full), pattern="\\.\\d.*", replacement="")==i, ], MARGIN=2, FUN=function(x) { median(x) } ))
rownames(try_table_nppathway_list_neg) <- unique(gsub(x=rownames(try_table_nppathway_list_neg_full), pattern="\\.\\d.*", replacement=""))
rownames(try_table_nppathway_list_neg) <- gsub(x=rownames(try_table_nppathway_list_neg), pattern="GOT|OEL", replacement="SWE")
rownames(try_table_nppathway_list_neg) <- gsub(x=rownames(try_table_nppathway_list_neg), pattern="OE", replacement="SWE")
rownames(try_table_nppathway_list_neg) <- gsub(x=rownames(try_table_nppathway_list_neg), pattern="HAL", replacement="GER")

## -----------------------------------------------------------------------------
try_table_nppathway_list <- cbind(try_table_nppathway_list_pos, try_table_nppathway_list_neg)
try_table_nppathway_list <- cbind(sapply(unique(colnames(try_table_nppathway_list)[duplicated(colnames(try_table_nppathway_list))]), function(x) rowSums(try_table_nppathway_list[,grepl(paste(x, "$", sep=""), colnames(try_table_nppathway_list))])), try_table_nppathway_list[,!duplicated(colnames(try_table_nppathway_list)) & ! duplicated(colnames(try_table_nppathway_list), fromLast=TRUE)])

## -----------------------------------------------------------------------------
try_table_npsuperclass_list_pos_full <- make_classes_at_level(level=2, div_samples=marchantiales$div_npclasses_samples_pos)

## -----------------------------------------------------------------------------
try_table_npsuperclass_list_pos <- NULL
for (i in unique(gsub(x=rownames(try_table_npsuperclass_list_pos_full), pattern="\\.\\d.*", replacement=""))) try_table_npsuperclass_list_pos <- rbind(try_table_npsuperclass_list_pos, apply(X=try_table_npsuperclass_list_pos_full[gsub(x=rownames(try_table_npsuperclass_list_pos_full), pattern="\\.\\d.*", replacement="")==i, ], MARGIN=2, FUN=function(x) { median(x) } ))
rownames(try_table_npsuperclass_list_pos) <- unique(gsub(x=rownames(try_table_npsuperclass_list_pos_full), pattern="\\.\\d.*", replacement=""))
rownames(try_table_npsuperclass_list_pos) <- gsub(x=rownames(try_table_npsuperclass_list_pos), pattern="GOT|OEL", replacement="SWE")
rownames(try_table_npsuperclass_list_pos) <- gsub(x=rownames(try_table_npsuperclass_list_pos), pattern="OE", replacement="SWE")
rownames(try_table_npsuperclass_list_pos) <- gsub(x=rownames(try_table_npsuperclass_list_pos), pattern="HAL", replacement="GER")

## -----------------------------------------------------------------------------
try_table_npsuperclass_list_neg_full <- make_classes_at_level(level=2, div_samples=marchantiales$div_npclasses_samples_neg)

## -----------------------------------------------------------------------------
try_table_npsuperclass_list_neg <- NULL
for (i in unique(gsub(x=rownames(try_table_npsuperclass_list_neg_full), pattern="\\.\\d.*", replacement=""))) try_table_npsuperclass_list_neg <- rbind(try_table_npsuperclass_list_neg, apply(X=try_table_npsuperclass_list_neg_full[gsub(x=rownames(try_table_npsuperclass_list_neg_full), pattern="\\.\\d.*", replacement="")==i, ], MARGIN=2, FUN=function(x) { median(x) } ))
rownames(try_table_npsuperclass_list_neg) <- unique(gsub(x=rownames(try_table_npsuperclass_list_neg_full), pattern="\\.\\d.*", replacement=""))
rownames(try_table_npsuperclass_list_neg) <- gsub(x=rownames(try_table_npsuperclass_list_neg), pattern="GOT|OEL", replacement="SWE")
rownames(try_table_npsuperclass_list_neg) <- gsub(x=rownames(try_table_npsuperclass_list_neg), pattern="OE", replacement="SWE")
rownames(try_table_npsuperclass_list_neg) <- gsub(x=rownames(try_table_npsuperclass_list_neg), pattern="HAL", replacement="GER")

## -----------------------------------------------------------------------------
try_table_npsuperclass_list <- cbind(try_table_npsuperclass_list_pos, try_table_npsuperclass_list_neg)
try_table_npsuperclass_list <- cbind(sapply(unique(colnames(try_table_npsuperclass_list)[duplicated(colnames(try_table_npsuperclass_list))]), function(x) rowSums(try_table_npsuperclass_list[,grepl(paste(x, "$", sep=""), colnames(try_table_npsuperclass_list))])), try_table_npsuperclass_list[,!duplicated(colnames(try_table_npsuperclass_list)) & ! duplicated(colnames(try_table_npsuperclass_list), fromLast=TRUE)])

## -----------------------------------------------------------------------------
try_table$`Molecules in Biosynthetic pathways` <- rowSums(try_table_npsuperclass_list)
try_table_colnames_1st <- c(try_table_colnames_1st, "Molecules in Biosynthetic pathways")
try_table_colnames_2nd <- c(try_table_colnames_2nd, "")

## -----------------------------------------------------------------------------
try_table_selected_nppathway <- try_table_nppathway_list
try_table <- cbind(try_table, try_table_selected_nppathway)
try_table_colnames_1st <- c(try_table_colnames_1st, colnames(try_table_selected_nppathway))
try_table_colnames_2nd <- c(try_table_colnames_2nd, rep("", times=length(colnames(try_table_selected_nppathway))))

## -----------------------------------------------------------------------------
try_table_selected_npsuperclass <- try_table_npsuperclass_list[, grep(x=colnames(try_table_npsuperclass_list), pattern=".*; ", perl=TRUE)]
try_table <- cbind(try_table, try_table_selected_npsuperclass)
try_table_colnames_1st <- c(try_table_colnames_1st, gsub(x=colnames(try_table_selected_npsuperclass), pattern="; .*", replacement=""))
try_table_colnames_2nd <- c(try_table_colnames_2nd, gsub(x=colnames(try_table_selected_npsuperclass), pattern=".*; ", replacement=""))

## -----------------------------------------------------------------------------
try_table_mdes_list_full <- marchantiales$mdes_list

## -----------------------------------------------------------------------------
try_table_mdes_list <- NULL
for (i in unique(gsub(x=rownames(try_table_mdes_list_full), pattern="\\.\\d.*", replacement=""))) try_table_mdes_list <- rbind(try_table_mdes_list, apply(X=try_table_mdes_list_full[gsub(x=rownames(try_table_mdes_list_full), pattern="\\.\\d.*", replacement="")==i, ], MARGIN=2, FUN=function(x) { median(x) } ))
rownames(try_table_mdes_list) <- unique(gsub(x=rownames(try_table_mdes_list_full), pattern="\\.\\d.*", replacement=""))
rownames(try_table_mdes_list) <- gsub(x=rownames(try_table_mdes_list), pattern="GOT|OEL", replacement="SWE")
rownames(try_table_mdes_list) <- gsub(x=rownames(try_table_mdes_list), pattern="OE", replacement="SWE")
rownames(try_table_mdes_list) <- gsub(x=rownames(try_table_mdes_list), pattern="HAL", replacement="GER")

## -----------------------------------------------------------------------------
try_table_mdes_high_level <- make.high_level.descriptor(colnames(marchantiales$mdes_list))

## -----------------------------------------------------------------------------
try_table_mdes_low_level <- colnames(marchantiales$mdes_list)

## -----------------------------------------------------------------------------
try_table <- cbind(try_table, try_table_mdes_list)

try_table_colnames_1st <- c(try_table_colnames_1st, try_table_mdes_high_level)
try_table_colnames_2nd <- c(try_table_colnames_2nd, try_table_mdes_low_level)

## ----eval=FALSE---------------------------------------------------------------
#  write.table(x=rbind(try_table_colnames_1st, try_table_colnames_2nd, try_table), file="try_table_traits.txt", row.names=FALSE, col.names=FALSE, sep='\t')

