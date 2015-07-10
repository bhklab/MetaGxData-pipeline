.libPaths(c("/mnt/work1/users/bhklab/Rlib", .libPaths()))

task.id <- as.integer(Sys.getenv("SGE_TASK_ID"))
set.seed(6400 + task.id * 100)

gene.set.size <- task.id

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

subtype <- args[1]

out.dir <- "randomsigs_out_500genesets_100resamples_500genes"

library(knitr)
library(gdata)
library(annotate)
library(ggplot2)
library(xtable)
library(saps)
library(genefu)
library(hgu133plus2.db)

load("pooled.eset.over.5000.genes.RData")

# Print patients on first two principal components

#pc.out <- prcomp(t(exprs(pooled.eset.over.5000.genes)))
#scatterplot.data <- as.data.frame(pc.out$x[,c(1,2)])
#scatterplot.data$data.source <- pooled.eset.over.5000.genes$data.source
#pca.scatterplot <- ggplot(scatterplot.data, aes(x=PC1, y=PC2, colour=data.source)) + 
#  geom_point(shape=1) +
#  scale_colour_hue(l=50) +
#  ggtitle("PCA")
#ggsave(pca.scatterplot, filename = "pca.brca.datasets.png", width=13, height=13)
genes.in.pooled.intersect <- rownames(fData(pooled.eset.over.5000.genes))
genes.in.pooled.intersect <- sub("^", "geneid.", genes.in.pooled.intersect)

nki.sample.size <- 295

.getPVals <- function(eset, random=FALSE) {
  log.rank.pvals <- sapply(1:500, function(x) {
    print(x)
    current.eset <- eset
    if(random == TRUE) {
       current.eset <- current.eset[,sample(1:ncol(exprs(current.eset)), 1000)]
    }
    
    random.gene.indices <- sample(1:nrow(exprs(current.eset)), gene.set.size)
    resampled.pvals <- sapply(1:100, function(y) {
      samples.to.include <- sample(1:ncol(exprs(current.eset)), nki.sample.size, replace=FALSE)
      current.eset <- current.eset[random.gene.indices,samples.to.include]
      expression.matrix <- t(exprs(current.eset))
      pc.out <- prcomp(expression.matrix)
      pc1 <- pc.out$x[,1]
      quantiles <- cut(pc1, breaks=quantile(pc1, probs=seq(0,1,length=2+1)), include.lowest=TRUE)
      levels(quantiles) <- c("Low", "High")
      eset.surv <- Surv(time=current.eset$days_to_death, event=current.eset$vital_status == "deceased")
      coxph.out <- coxph(eset.surv ~ quantiles + strata(current.eset$data.source))
      pval <- summary(coxph.out)$sctest["pvalue"]
      return(pval)
    })
    return(mean(resampled.pvals))
  })
  return(log.rank.pvals)
}

if(subtype=="All") {
	pvals <- .getPVals(pooled.eset.over.5000.genes)
} else if (subtype=="Random") {
	pvals <- .getPVals(pooled.eset.over.5000.genes[,seq(1, ncol(exprs(pooled.eset.over.5000.genes)), 4)])
} else {
	pvals <- .getPVals(pooled.eset.over.5000.genes, random=TRUE)
}
 
var.name <- paste0("pvals.", subtype, ".", gene.set.size)

assign(var.name, pvals)

save(list=var.name, file=paste0(out.dir, "/", var.name, ".RData"))
