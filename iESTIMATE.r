


# ############################## iESTIMATE ##############################



# ---------- Preparations ----------
library(devtools)
library(roxygen2)
library(usethis)
library(testthat)



# ---------- Create new R package ----------
setwd("./")

# Create new package
#usethis::create_package("iESTIMATE")
#devtools::create(path="./")

# Add function from R files
usethis::use_r("functions.R")

# Build documentation
devtools::document()

# README
#usethis::use_readme_md()

# Unit tests
#usethis::use_test("shannon.diversity")
#devtools::test()

# Check requirements
devtools::check()

# Load functions
load_all()

# Add reference data
marchantiales <- list()
marchantiales[["comp_list"]] <- read.table(file="use-cases/marchantiales/comp_list.tsv", header=TRUE, sep="\t", quote="\"", comment.char="", fill=FALSE, dec=".", stringsAsFactors=FALSE)
marchantiales[["bina_list"]] <- read.table(file="use-cases/marchantiales/bina_list.tsv", header=TRUE, sep="\t", quote="\"", comment.char="", fill=FALSE, dec=".", stringsAsFactors=FALSE)
marchantiales[["uniq_list"]] <- read.table(file="use-cases/marchantiales/uniq_list.tsv", header=TRUE, sep="\t", quote="\"", comment.char="", fill=FALSE, dec=".", stringsAsFactors=FALSE)
marchantiales[["class_list"]] <- read.table(file="use-cases/marchantiales/class_list.tsv", header=TRUE, sep="\t", quote="\"", comment.char="", fill=FALSE, dec=".", stringsAsFactors=FALSE)
marchantiales[["subclass_list"]] <- read.table(file="use-cases/marchantiales/subclass_list.tsv", header=TRUE, sep="\t", quote="\"", comment.char="", fill=FALSE, dec=".", stringsAsFactors=FALSE)
marchantiales[["superclass_list"]] <- read.table(file="use-cases/marchantiales/superclass_list.tsv", header=TRUE, sep="\t", quote="\"", comment.char="", fill=FALSE, dec=".", stringsAsFactors=FALSE)
marchantiales[["mdes_list"]] <- read.table(file="use-cases/marchantiales/mdes_list.tsv", header=TRUE, sep="\t", quote="\"", comment.char="", fill=FALSE, dec=".", stringsAsFactors=FALSE)
marchantiales[["div_classes"]] <- read.table(file="use-cases/marchantiales/div_classes.tsv", header=TRUE, sep="\t", quote="\"", comment.char="", fill=FALSE, dec=".", stringsAsFactors=FALSE)
marchantiales[["div_subclasses"]] <- read.table(file="use-cases/marchantiales/div_subclasses.tsv", header=TRUE, sep="\t", quote="\"", comment.char="", fill=FALSE, dec=".", stringsAsFactors=FALSE)
marchantiales[["div_superclasses"]] <- read.table(file="use-cases/marchantiales/div_superclasses.tsv", header=TRUE, sep="\t", quote="\"", comment.char="", fill=FALSE, dec=".", stringsAsFactors=FALSE)
marchantiales[["div_npclasses"]] <- read.table(file="use-cases/marchantiales/div_npclasses.tsv", header=TRUE, sep="\t", quote="\"", comment.char="", fill=FALSE, dec=".", stringsAsFactors=FALSE)
marchantiales[["div_nppathway"]] <- read.table(file="use-cases/marchantiales/div_nppathway.tsv", header=TRUE, sep="\t", quote="\"", comment.char="", fill=FALSE, dec=".", stringsAsFactors=FALSE)
marchantiales[["div_npclasses_samples_pos"]] <- read.table(file="use-cases/marchantiales/div_npclasses_samples_pos.tsv", header=TRUE, sep="\t", quote="\"", comment.char="", fill=FALSE, dec=".", stringsAsFactors=FALSE)
marchantiales[["div_npclasses_samples_neg"]] <- read.table(file="use-cases/marchantiales/div_npclasses_samples_neg.tsv", header=TRUE, sep="\t", quote="\"", comment.char="", fill=FALSE, dec=".", stringsAsFactors=FALSE)
marchantiales[["model_div"]] <- read.table(file="use-cases/marchantiales/model_div.tsv", header=TRUE, sep="\t", quote="\"", comment.char="", fill=FALSE, dec=".", stringsAsFactors=FALSE)
marchantiales[["metadata"]] <- read.table(file="use-cases/marchantiales/metadata.tsv", header=TRUE, sep="\t", quote="\"", comment.char="", fill=FALSE, dec=".", stringsAsFactors=FALSE)
marchantiales[["phylo_trnLF"]] <- phangorn::read.phyDat(file="use-cases/marchantiales/phylo_trnLF.fasta", format="fasta")
marchantiales[["char_list"]] <- read.table(file="use-cases/marchantiales/char_list.tsv", header=TRUE, sep="\t", quote="\"", comment.char="", fill=FALSE, dec=".", stringsAsFactors=FALSE)
marchantiales[["trait_list"]] <- read.table(file="use-cases/marchantiales/trait_list.tsv", header=TRUE, sep="\t", quote="\"", comment.char="", fill=FALSE, dec=".", stringsAsFactors=FALSE)
usethis::use_data(marchantiales, overwrite=TRUE)
#use_r("marchantiales")

# Vignettes
#usethis::use_vignette("marchantiales")
devtools::build_vignettes()

# Commit using GIT
#usethis::use_git()

# GITHub
#usethis::use_git_remote("origin", url = NULL, overwrite = TRUE)
#usethis::use_github()

# Install package
devtools::install_github("ipb-halle/iESTIMATE")


