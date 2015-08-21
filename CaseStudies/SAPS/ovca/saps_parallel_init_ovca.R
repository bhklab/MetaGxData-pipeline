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
source("../datasetMerging.R")
source("../stripWhiteSpace.R")

# TODO? order by publication date. Note that EXPO was set to the year it was public on GEO (2005).

esets <- lapply(esets, function(x) {
  factor.indices <- sapply(pData(x), is.factor)
  pData(x)[factor.indices] <- lapply(pData(x)[factor.indices], as.character)
  return(x)
})

# only keep patients with survival data
esets <- lapply(esets, function(eset) eset[,!is.na(eset$recurrence_status) & !is.na(eset$days_to_tumor_recurrence)])

# Remove TCGA RNASeq
esets <- esets[names(esets) != "TCGA.RNASeqV2"]

esets <- lapply(esets, function(eset) eset[apply(exprs(eset), 1, function(x) all(!is.na(x))),])

## Remove datasets that are empty
esets <- esets[sapply(esets, function(x) ncol(exprs(x)) > 0)]

# Only keep datasets with at least 10000 genes
esets <- esets[sapply(esets, function(x) nrow(x) > 10000)]

for(i in 1:length(esets)) {
  expression.matrix <- t(exprs(esets[[i]]))
  annot <- fData(esets[[i]])
  colnames(annot)[which(colnames(annot) == "EntrezGene.ID")] <- "entrezgene"
  angio <- genefu::ovcAngiogenic(data = expression.matrix, annot=annot, gmap="entrezgene", do.mapping = TRUE)
  esets[[i]]$Bentink.subtypes <- angio$subtype$subtype
}

pooled.ovca.eset.intersecting.genes <- datasetMerging(esets, method='intersect', nthread=parallel::detectCores())

save(pooled.ovca.eset.intersecting.genes, file="pooled.ovca.eset.intersecting.genes.RData")
