library(doppelgangR)
library(Biobase)
filenames <- paste("./esets/mapped_esets2/", list.files("./esets/mapped_esets2"),sep="")
lapply(filenames, load, .GlobalEnv)

datasets <- read.csv("datasets.csv")
dataset.names <- datasets$Dataset
esets <- lapply(as.character(dataset.names), get)
names(esets) <- dataset.names
# save(esets, file="./esets/esets.rda)
# load("./esets/esets.rda")
esets.mapped <- list()
for(eset in esets){
  Biobase::exprs(eset) <- exprs(eset)[fData(eset)$best_probe,]
  Biobase::fData(eset) <- fData(eset)[fData(eset)$best_probe,]
  rownames(fData(eset)) <- rownames(exprs(eset)) <- paste("geneid.", fData(eset)$ENTREZID, sep="")
  esets.mapped <- c(esets.mapped, eset)
}
names(esets.mapped) <- names(esets)
save(esets.mapped, file="./esets/esetsMapped.rda")


## bonf.prob=2.0 to get all correlations of all pairs, then set cut off according to distribution
#duplicates <- doppelgangR(esets.mapped, outlierFinder.expr.args=list(bonf.prob=2.0, transFun=atanh, tail="upper"), phenoFinder.args = NULL, cache.dir=NULL)
#save(duplicates, file ="./esets/duplicatesMETABRICandTCGA.rda")
