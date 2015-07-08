days.per.month <- 30.4368
days.per.year <- 365.242
par.original <- par()
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

source("GeneSetListGenerator.R")
  genes.in.pooled.intersect <- as.character(fData(pooled.eset.intersecting.genes)$EntrezGene.ID)
GO.gene.sets.intersect <- lapply(CompleteLists, function(x) intersect(x, genes.in.pooled.intersect))
GO.gene.sets.intersect <- GO.gene.sets.intersect[sapply(GO.gene.sets.intersect, function(x) length(x) >= 15)]
GO.gene.sets.intersect <- lapply(GO.gene.sets.intersect, function(x) gsub("^", "geneid.", x))
rowMax <- max(sapply(GO.gene.sets.intersect, length))
# now create output matrix by lengthening rows
GO.gene.sets.matrix <- do.call(rbind, lapply(GO.gene.sets.intersect, function(x){
  length(x) <- rowMax
  x
}))

i <- as.integer(Sys.getenv("SGE_TASK_ID"))

GO.gene.sets.matrix <- GO.gene.sets.matrix[i,,drop=FALSE]

out <- saps(GO.gene.sets.matrix, t(exprs(pooled.eset.intersecting.genes)), pooled.eset.intersecting.genes$days_to_death, as.integer(pooled.eset.intersecting.genes$vital_status == "deceased"))

save(out, file=paste0(rownames(GO.gene.sets.matrix), ".RData"))

