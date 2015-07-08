library(knitr)
library(gdata)
library(annotate)
library(ggplot2)
library(xtable)
library(saps)
library(genefu)
library(hgu133plus2.db)

source(system.file("extdata", "patientselection.config", package="MetaGxBreast"))
source(system.file("extdata", "createEsetList.R", package="MetaGxBreast"))

## TODO: extract this function to not require MetaGx
source("~/repos/MetaGx/R/datasetMerging.R")
source("~/repos/MetaGx/R/getSubtype.R")
source("~/repos/MetaGx/R/setSubtype.R")
source("~/repos/MetaGx/R/subtypeClassification.R")
source("~/repos/MetaGx/R/stripWhiteSpace.R")

# TODO? order by publication date. Note that EXPO was set to the year it was public on GEO (2005).

esets <- lapply(esets, function(x) {
  factor.indices <- sapply(pData(x), is.factor)
  pData(x)[factor.indices] <- lapply(pData(x)[factor.indices], as.character)
  return(x)
})

# only keep patients with survival data
esets <- lapply(esets, function(eset) eset[,!is.na(eset$days_to_death) & !is.na(eset$vital_status)])

# The gene expression matrix of UCSF is over 8% NAs
esets <- esets[-which(names(esets) == "UCSF")]

# For TCGA, remove the 169 genes with NA values
esets$TCGA <- esets$TCGA[apply(exprs(esets$TCGA), 1, function(x) all(!is.na(x))),]

esets$NKI <- esets$NKI[apply(exprs(esets$NKI), 1, function(x) sum(is.na(x)) < 20),]
esets$NKI <- esets$NKI[,apply(exprs(esets$NKI), 2, function(x) sum(is.na(x))) < 5]
esets$NKI <- esets$NKI[apply(exprs(esets$NKI), 1, function(x) sum(is.na(x)) == 0),]

## Remove datasets that are empty
esets <- esets[sapply(esets, function(x) ncol(exprs(x)) > 0)]

esets <- lapply(esets, function(x) {
  x <- subtypeClassification(x, model = "scmod2")
  x$subtype <- experimentData(x)@other$class
  return(x)
})

pooled.eset.intersecting.genes <- datasetMerging(esets, method='intersect', nthread=parallel::detectCores())

save(pooled.eset.intersecting.genes, file="pooled.eset.intersecting.genes.RData")