## load the expression sets

filenames <- paste("./esets/mapped_esets2/", list.files("./esets/mapped_esets2"),sep="")
lapply(filenames, load, .GlobalEnv)

datasets <- read.csv("datasets.csv")
dataset.names <- datasets$Dataset
esets <- lapply(as.character(dataset.names), get)
names(esets) <- dataset.names

genes <- list()
for ( dataset in dataset.names){
	eset <- get(dataset)
	ds.genes <- unique(fData(eset)$SYMBOL)
	if(length(grep("///", ds.genes))>0){
	ds.genes <- ds.genes[-grep("///",ds.genes)]
	}

	genes[[dataset]] <- ds.genes
}

SYMBOLS <- lapply(genes, function(i){as.character(i)})
SYMBOLS <- unique(unlist(SYMBOLS))

gene.map <- NULL
for (dataset in dataset.names){
    eset <- get(dataset)
    tmp <- rep(0, each = length(SYMBOLS))
    tmp[is.element(SYMBOLS, fData(eset)$SYMBOL)] <- 1
    gene.map <- cbind(gene.map, tmp)
}
colnames(gene.map) <- dataset.names
rownames(gene.map) <- SYMBOLS
save(gene.map, file="./esets/genemap.rda")