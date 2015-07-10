library(reshape)

for(filename in paste0("randomsigs_out_1000genesets_1000resamples/", grep(".RData", list.files("randomsigs_out_1000genesets_1000resamples/"), value = TRUE))) {
  load(filename)
}

gene.set.sizes <- seq(10, 80, 10)
#
#nki.sample.size <- 295
#for(gene.set.size in gene.set.sizes) {
#  random.gene.indices <- sample(1:nrow(exprs(eset)), gene.set.size)
# for(j in 1:9) { 
#  out <- lapply(c(pooled.eset.intersecting.genes, 
#           pooled.eset.intersecting.genes[,seq(1, ncol(exprs(pooled.eset.intersecting.genes)), 4)],
#           pooled.eset.intersecting.genes[,sample(1:ncol(exprs(pooled.eset.intersecting.genes)), 1000)]), function(eset) {
#    resampled.pvals <- sapply(1:100, function(x) {
#      num.to.sample <- round(295 * (as.numeric(table(eset$data.source)) / length(eset$data.source)))
#      #sampled.vals <- lapply(1:length(num.to.sample), function(i) {
#      #  current.num.to.sample <- num.to.sample[i]
#      #  current.indices <- 1:ncol(exprs(eset))
#      #  current.indices <- current.indices[eset$data.source==levels(eset$data.source)[i]]
#      #  sample(current.indices, current.num.to.sample, replace=FALSE)
#      #})
#      #samples.to.include <- do.call(c, sampled.vals)
#      samples.to.include <- sample(1:ncol(exprs(eset)), nki.sample.size, replace=FALSE)
#      
#      current.eset <- eset[random.gene.indices,samples.to.include]
#      expression.matrix <- t(exprs(current.eset))
#      pc.out <- prcomp(expression.matrix)
#      pc1 <- pc.out$x[,1]
#      quantiles <- cut(pc1, breaks=quantile(pc1, probs=seq(0,1,length=2+1)), include.lowest=TRUE)
#      levels(quantiles) <- c("Low", "High")
#      eset.surv <- Surv(time=current.eset$days_to_death, event=current.eset$vital_status == "deceased")
#      coxph.out <- coxph(eset.surv ~ quantiles + strata(current.eset$data.source))
#      pval <- summary(coxph.out)$sctest["pvalue"]
#      return(pval)
#    })
#    })
#  boxplot(out)
# }
#  return(mean(resampled.pvals))
#}



pval.list.by.size <- list()

for(gene.set.size in gene.set.sizes) {
  current.objnames <- grep(paste0("\\.", gene.set.size, "$"), ls(), value = TRUE)
  pval.list <- mget(sort(current.objnames))
  pval.list <- pval.list[order(names(pval.list))]
  if(length(pval.list) != 6) {
    stop(paste0("For size ", gene.set.size, ", could not find all RData objects"))
  }
  rm(list=current.objnames)
  names(pval.list) <- sub(paste0("\\.", gene.set.size, "$"), "", names(pval.list))
  names(pval.list) <- sub(paste0("^pvals\\."), "", names(pval.list))
  pval.list.by.size[[length(pval.list.by.size)+1]] <- pval.list
  names(pval.list.by.size)[length(pval.list.by.size)] <- paste0("GeneSet.size.", gene.set.size)
}

pval.list.by.size <- lapply(pval.list.by.size, as.data.frame)

pval.means.list <- lapply(pval.list.by.size, function(x) sapply(x, mean))
pval.mean.matrix <- as.matrix(do.call(rbind, pval.means.list))

proportion.below.cutoff.list <- lapply(pval.list.by.size, function(x) sapply(x, function(y) mean(y<0.20)))
proportion.below.cutoff.matrix <- as.matrix(do.call(rbind, proportion.below.cutoff.list))

.getHeatmap <- function(stat.matrix, stat.name, cluster=TRUE) {
  if(cluster==TRUE) {
    ord <- hclust(dist(stat.matrix, method="euclidean"))$order
  }
  stat.vals.m <- melt(stat.matrix)
  colnames(stat.vals.m) <- c("GeneSet", "Subtype", "value")
  
  if(cluster==TRUE) {
    #Set order of factor levels, which is reflected in the heatmap order for genesets
    stat.vals.m$GeneSet <- factor(stat.vals.m$GeneSet, levels = levels(stat.vals.m$GeneSet)[ord])
  } else {
    #Retain original order of rows
    stat.vals.m$GeneSet <- factor(stat.vals.m$GeneSet, levels = levels(stat.vals.m$GeneSet)[match(rownames(stat.matrix), levels(stat.vals.m$GeneSet))])
  }
  
  p <- ggplot(stat.vals.m, aes(Subtype, GeneSet)) + 
    geom_tile(aes(fill = value), colour = "white") + 
    scale_fill_gradient(name=stat.name, low="#58A4DE", high="black") + 
    ggtitle("Significance Analysis of Prognostic Signatures") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
  
  return(p)
}

.getHeatmap(proportion.below.cutoff.matrix, "Proportion with p < 0.10", cluster=FALSE) + theme(axis.text.y = element_text(size=20)) + geom_text(label=as.character(proportion.below.cutoff.matrix), colour="white")
.getHeatmap(pval.mean.matrix, "Mean p-value", cluster=FALSE) + theme(axis.text.y = element_text(size=20))  + geom_text(label=sprintf("%.3f", pval.mean.matrix), colour="white")

par(mfrow=c(2,4))
lapply(names(pval.list.by.size), function(pval.list.name) {
  boxplot(value ~ variable, melt(pval.list.by.size[pval.list.name]), main=pval.list.name)
  })
