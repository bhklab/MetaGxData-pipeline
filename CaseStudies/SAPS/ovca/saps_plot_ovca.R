library(reshape)

for(filename in paste0("saps_output_ovca/", grep(".RData", list.files("saps_output_ovca/"), value = TRUE))) {
  load(filename)
}

all.subtypes.saps.output <- mget(sort(grep("_All", ls(), value = TRUE)))
rm(list=grep("_All$", ls(), value = TRUE))
names(all.subtypes.saps.output) <- sub("_All$", "", names(all.subtypes.saps.output))

angiogenic.saps.output <- mget(sort(grep("_Angiogenic", ls(), value = TRUE)))
rm(list=grep("_Angiogenic$", ls(), value = TRUE))
names(angiogenic.saps.output) <- sub("_Angiogenic$", "", names(angiogenic.saps.output))

nonAngiogenic.saps.output <- mget(sort(grep("_nonAngiogenic", ls(), value = TRUE)))
rm(list=grep("_nonAngiogenic$", ls(), value = TRUE))
names(nonAngiogenic.saps.output) <- sub("_nonAngiogenic$", "", names(nonAngiogenic.saps.output))

all.out <- list(All=all.subtypes.saps.output,
                Angiogenic=angiogenic.saps.output,
                nonAngiogenic=nonAngiogenic.saps.output
                )
# order alphabetically
all.out <- lapply(all.out, function(x) x[order(names(x))])

intersecting.gene.set.names <- scan("../intersecting.gene.sets.txt", what=character(0))
all.out <- lapply(all.out, function(x) x[intersecting.gene.set.names])

# Make sure names are all equal
my.names <- lapply(all.out, names)
if(!( all.equal(my.names[[1]], my.names[[2]])
    && all.equal(my.names[[1]], my.names[[3]]))) {
  stop("Names of saps output are not equal - are all RData files present?")
}

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

p.random.plot <- .getHeatmap(-log10(p.random), "-log(p-random)")

ggsave(p.random.plot, filename = "p-random_ovca_intersect.png", width=13, height=15)
