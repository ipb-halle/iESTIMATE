


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
tenerife <- list()
tenerife[["comp_list"]] <- read.table(file="use-cases/tenerife/comp_list.tsv", header=TRUE, sep="\t", quote="\"", comment.char="", fill=FALSE, dec=".", stringsAsFactors=FALSE)
tenerife[["bina_list"]] <- read.table(file="use-cases/tenerife/bina_list.tsv", header=TRUE, sep="\t", quote="\"", comment.char="", fill=FALSE, dec=".", stringsAsFactors=FALSE)
tenerife[["uniq_list"]] <- read.table(file="use-cases/tenerife/uniq_list.tsv", header=TRUE, sep="\t", quote="\"", comment.char="", fill=FALSE, dec=".", stringsAsFactors=FALSE)
tenerife[["class_list"]] <- read.table(file="use-cases/tenerife/class_list.tsv", header=TRUE, sep="\t", quote="\"", comment.char="", fill=FALSE, dec=".", stringsAsFactors=FALSE)
tenerife[["subclass_list"]] <- read.table(file="use-cases/tenerife/subclass_list.tsv", header=TRUE, sep="\t", quote="\"", comment.char="", fill=FALSE, dec=".", stringsAsFactors=FALSE)
tenerife[["superclass_list"]] <- read.table(file="use-cases/tenerife/superclass_list.tsv", header=TRUE, sep="\t", quote="\"", comment.char="", fill=FALSE, dec=".", stringsAsFactors=FALSE)
tenerife[["moldes_list"]] <- read.table(file="use-cases/tenerife/moldes_list.tsv", header=TRUE, sep="\t", quote="\"", comment.char="", fill=FALSE, dec=".", stringsAsFactors=FALSE)
usethis::use_data(tenerife, overwrite=TRUE)
#use_r("tenerife")

# Vignettes
#usethis::use_vignette("tenerife")
devtools::build_vignettes()

# Commit using GIT
#usethis::use_git()

# GITHub
#usethis::use_git_remote("origin", url = NULL, overwrite = TRUE)
#usethis::use_github()

# Install package
devtools::install_github("ipb-halle/iESTIMATE")


