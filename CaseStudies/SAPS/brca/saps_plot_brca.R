library(reshape)
library(ggplot2)

for(filename in paste0("saps_output_brca_rfs/", grep(".RData", list.files("saps_output_brca_rfs/"), value = TRUE))) {
  load(filename)
}

all.subtypes.saps.output <- mget(sort(grep("_All", ls(), value = TRUE)))
rm(list=grep("_All$", ls(), value = TRUE))
names(all.subtypes.saps.output) <- sub("_All$", "", names(all.subtypes.saps.output))
basal.saps.output <- mget(sort(grep("_Basal", ls(), value = TRUE)))
rm(list=grep("_Basal$", ls(), value = TRUE))
names(basal.saps.output) <- sub("_Basal$", "", names(all.subtypes.saps.output))
her2.saps.output <- mget(sort(grep("_Her2", ls(), value = TRUE)))
rm(list=grep("_Her2$", ls(), value = TRUE))
names(her2.saps.output) <- sub("_Her2$", "", names(all.subtypes.saps.output))
lumB.saps.output <- mget(sort(grep("_LumB", ls(), value = TRUE)))
rm(list=grep("_LumB$", ls(), value = TRUE))
names(lumB.saps.output) <- sub("_LumB$", "", names(all.subtypes.saps.output))
lumA.saps.output <- mget(sort(grep("_LumA", ls(), value = TRUE)))
rm(list=grep("_LumA$", ls(), value = TRUE))
names(lumA.saps.output) <- sub("_LumA$", "", names(all.subtypes.saps.output))

all.out <- list(All=all.subtypes.saps.output,
                Basal=basal.saps.output,
                Her2=her2.saps.output,
                LumA=lumA.saps.output,
                LumB=lumB.saps.output
                )
# order alphabetically
all.out <- lapply(all.out, function(x) x[order(names(x))])

#intersecting.gene.set.names <- scan("../intersecting.gene.sets.txt", what=character(0))
#all.out <- lapply(all.out, function(x) x[intersecting.gene.set.names])

# Make sure names are all equal
my.names <- lapply(all.out, names)
if(!(all.equal(my.names[[1]], my.names[[2]])
    && all.equal(my.names[[1]], my.names[[3]])
    && all.equal(my.names[[1]], my.names[[4]])
    && all.equal(my.names[[1]], my.names[[5]]))) {
  stop("Names of saps output are not equal - are all RData files present?")
}

saps.scores <- lapply(all.out, function(x) lapply(x, function(y) y$genesets[[1]]$saps_unadjusted["saps_score"]))
saps.scores <- sapply(saps.scores, unlist)
saps.scores <- abs(saps.scores)
rownames(saps.scores) <- sub(".saps_score$", "", rownames(saps.scores))

p.random <- lapply(all.out, function(x) lapply(x, function(y) y$genesets[[1]]$saps_unadjusted["p_random"]))
p.random <- sapply(p.random, unlist)
rownames(p.random) <- sub(".p_random$", "", rownames(p.random))

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
    scale_fill_gradient(name=stat.name, low="#58A4DE", high="#000000") + 
    ggtitle("Significance Analysis of Prognostic Signatures") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
  
  return(p)
}
.getHeatmapCategorical <- function(stat.matrix, stat.name, cluster=TRUE, breaks=c(0,0.05,0.25,1)) {
  if(cluster==TRUE) {
    ord <- hclust(dist(stat.matrix, method="euclidean"))$order
  }
  original.rownames <- rownames(stat.matrix)
  stat.matrix <- apply(stat.matrix, 2, cut, breaks=breaks)
  rownames(stat.matrix) <- original.rownames
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
    scale_fill_discrete(name=stat.name) + 
    ggtitle("Significance Analysis of Prognostic Signatures") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
  
  return(p)
}

p.random.fdr <- apply(p.random, 2, p.adjust, method="fdr")
p.random.fdr.cutoffs <- apply(p.random.fdr, 2, cut, breaks=c(0,0.05, 0.25,1))
rownames(p.random.fdr.cutoffs) <- rownames(p.random.fdr)

p.random.plot <- .getHeatmap(-log10(p.random), "-log(p-random)")

p.random.plot.fdr.cutoffs <- .getHeatmapCategorical(p.random.fdr, "p-value fdr")

#ggsave(saps.plot, filename = "saps.png", width=13, height=15)
ggsave(p.random.plot, filename = "p-random_brca_fullgeneset.png", width=13, height=15)
