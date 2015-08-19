.libPaths(c("/mnt/work1/users/bhklab/Rlib", .libPaths()))
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

esets <- lapply(esets, function(x) {
  factor.indices <- sapply(pData(x), is.factor)
  pData(x)[factor.indices] <- lapply(pData(x)[factor.indices], as.character)
  return(x)
})

# only keep patients with rfs data
esets <- lapply(esets, function(eset) eset[,!is.na(eset$recurrence_status) & !is.na(eset$days_to_tumor_recurrence)  |  !is.na(eset$dmfs_status) & !is.na(eset$dmfs_days)])

esets <- lapply(esets, function(eset) eset[apply(exprs(eset), 1, function(x) all(!is.na(x))),])

## Remove datasets that are empty
esets <- esets[sapply(esets, function(x) ncol(exprs(x)) > 0)]

# Only keep datasets with at least 10000 genes
esets <- esets[sapply(esets, function(x) nrow(x) > 10000)]

esets <- lapply(esets, function(x) {
  x <- subtypeClassification(x, model = "scmod2")
  x$subtype <- experimentData(x)@other$class
  return(x)
})

pooled.eset.intersecting.genes <- datasetMerging(esets, method='intersect', nthread=parallel::detectCores())

# impute rfs with dmfs if rfs is not available
use.dmfs.logical <- is.na(pooled.eset.intersecting.genes$days_to_tumor_recurrence) & is.na(pooled.eset.intersecting.genes$recurrence_status) & !is.na(pooled.eset.intersecting.genes$dmfs_days) & !is.na(pooled.eset.intersecting.genes$dmfs_status)

pooled.eset.intersecting.genes$days_to_tumor_recurrence[use.dmfs.logical] <- pooled.eset.intersecting.genes$dmfs_days[use.dmfs.logical]
pooled.eset.intersecting.genes$recurrence_status[use.dmfs.logical] <- pooled.eset.intersecting.genes$dmfs_status[use.dmfs.logical]

save(pooled.eset.intersecting.genes, file="pooled.brca.eset.intersecting.genes.RData")
