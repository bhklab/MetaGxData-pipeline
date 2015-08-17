.libPaths(c("/mnt/work1/users/bhklab/Rlib", .libPaths()))

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

subtype <- args[1]

out.dir <- 

gene.set.index <- as.integer(Sys.getenv("SGE_TASK_ID"))

days.per.month <- 30.4368
days.per.year <- 365.242

library(knitr)
library(gdata)
library(annotate)
library(ggplot2)
library(xtable)
library(saps)
library(genefu)
library(hgu133plus2.db)

load("pooled.ovca.eset.intersecting.genes.RData")

source("../GeneSetListGenerator.R", chdir=TRUE)
genes.in.pooled.intersect <- as.character(fData(pooled.ovca.eset.intersecting.genes)$EntrezGene.ID)
GO.gene.sets.intersect <- lapply(CompleteLists, function(x) intersect(x, genes.in.pooled.intersect))
GO.gene.sets.intersect <- GO.gene.sets.intersect[sapply(GO.gene.sets.intersect, function(x) length(x) >= 15)]
GO.gene.sets.intersect <- lapply(GO.gene.sets.intersect, function(x) gsub("^", "geneid.", x))
rowMax <- max(sapply(GO.gene.sets.intersect, length))
# now create output matrix by lengthening rows
GO.gene.sets.matrix <- do.call(rbind, lapply(GO.gene.sets.intersect, function(x){
  length(x) <- rowMax
  x
}))

## Subset by the subtype to consider
if(subtype != "All") {
  pooled.ovca.eset.intersecting.genes <- pooled.ovca.eset.intersecting.genes[,pooled.ovca.eset.intersecting.genes$Bentink.subtypes == subtype]
}

GO.gene.sets.matrix <- GO.gene.sets.matrix[gene.set.index,,drop=FALSE]

out <- saps(
  GO.gene.sets.matrix, 
  t(exprs(pooled.ovca.eset.intersecting.genes)), 
  pooled.ovca.eset.intersecting.genes$days_to_death, 
  as.integer(pooled.ovca.eset.intersecting.genes$vital_status == "deceased"),
  compute_qvalue = FALSE)

var.name <- paste0(rownames(GO.gene.sets.matrix), "_", subtype)

assign(var.name, out)

save(list=var.name, file=paste0("saps_output_ovca/", var.name, ".RData"))
