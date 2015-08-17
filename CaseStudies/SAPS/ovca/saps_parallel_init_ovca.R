.libPaths(c("/mnt/work1/users/bhklab/Rlib", .libPaths()))
library(knitr)
library(gdata)
library(annotate)
library(ggplot2)
library(xtable)
library(saps)
library(genefu)
library(hgu133plus2.db)

source(system.file("extdata", "patientselection.config", package="MetaGxOvarian"))
source(system.file("extdata", "createEsetList.R", package="MetaGxOvarian"))

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

# Remove TCGA RNASeq
esets <- esets[names(esets) != "TCGA.RNASeqV2"]

# Remove GSE19829
esets <- esets[names(esets) != "GSE19829"]

#GSE51088
esets$GSE51088 <- esets$GSE51088[apply(exprs(esets$GSE51088), 1, function(x) all(!is.na(x))),]

#GSE8842
esets$GSE8842 <- esets$GSEGSE8842[apply(exprs(esets$GSE8842), 1, function(x) all(!is.na(x))),]

## Remove datasets that are empty
esets <- esets[sapply(esets, function(x) ncol(exprs(x)) > 0)]

pooled.ovca.eset.intersecting.genes <- datasetMerging(esets, method='intersect', nthread=parallel::detectCores())

save(pooled.ovca.eset.intersecting.genes, file="pooled.ovca.eset.intersecting.genes.RData")
