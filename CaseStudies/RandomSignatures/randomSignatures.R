.libPaths(c("/mnt/work1/users/bhklab/Rlib", .libPaths()))

task.id <- as.integer(Sys.getenv("SGE_TASK_ID"))
set.seed(6400 + task.id * 100)

gene.set.size <- task.id

out.dir <- "randomsigs_out"

library(knitr)
library(gdata)
library(annotate)
library(ggplot2)
library(xtable)
library(saps)
library(genefu)
library(hgu133plus2.db)

load("pooled.eset.intersecting.genes.RData")

# Print patients on first two principal components

#pc.out <- prcomp(t(exprs(pooled.eset.intersecting.genes)))
#scatterplot.data <- as.data.frame(pc.out$x[,c(1,2)])
#scatterplot.data$data.source <- pooled.eset.intersecting.genes$data.source
#pca.scatterplot <- ggplot(scatterplot.data, aes(x=PC1, y=PC2, colour=data.source)) + 
#  geom_point(shape=1) +
#  scale_colour_hue(l=50) +
#  ggtitle("PCA")
#ggsave(pca.scatterplot, filename = "pca.brca.datasets.png", width=13, height=13)
genes.in.pooled.intersect <- rownames(fData(pooled.eset.intersecting.genes))
genes.in.pooled.intersect <- sub("^", "geneid.", genes.in.pooled.intersect)

bootstrap.sample.size <- min(table(pooled.eset.intersecting.genes$subtype))

.getPVals <- function(eset) {
# now create output matrix by lengthening rows
  log.rank.pvals <- sapply(1:500, function(x) {
    random.gene.indices <- sample(1:nrow(exprs(eset)), gene.set.size)
    system.time(
    bootstrap.pvals <- sapply(1:500, function(x) {
      samples.to.include <- sample(1:ncol(exprs(eset)), bootstrap.sample.size, replace=TRUE)
      current.eset <- eset[random.gene.indices,samples.to.include]
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
    )
    return(mean(bootstrap.pvals))
  })
  return(pvals)
}

all.pvals <- .getPVals(pooled.eset.intersecting.genes)
basal.pvals <- .getPVals(pooled.eset.intersecting.genes[,pooled.eset.intersecting.genes$subtype=="Basal"])
her2.pvals <- .getPVals(pooled.eset.intersecting.genes[,pooled.eset.intersecting.genes$subtype=="Her2"])
lumA.pvals <- .getPVals(pooled.eset.intersecting.genes[,pooled.eset.intersecting.genes$subtype=="LumA"])
lumB.pvals <- .getPVals(pooled.eset.intersecting.genes[,pooled.eset.intersecting.genes$subtype=="LumB"])
  
pvals.list <- list(All=all.pvals.out,
                   Basal=basal.pvals,
                   Her2=her2.pvals,
                   LumA=lumA.pvals,
                   LumB=lumB.pvals)

var.name <- paste0("pvals.list.", gene.set.size)

assign(var.name, pvals.list)

save(list=var.name, file=paste0(out.dir, "/", var.name, ".RData"))