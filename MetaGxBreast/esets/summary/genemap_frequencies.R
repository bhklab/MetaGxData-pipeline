## load the expression sets



source("./patientselection.config")
min.number.of.genes <- 0
rm(remove.duplicates)
source("./R/createEsetList2.R")


genes <- list()
for(i in 1:length(esets)){
	eset <- esets[[i]]
	ds.genes <- unique(fData(eset)$EntrezGene.ID)
	if(length(grep("///", ds.genes))>0){
		ds.genes <- ds.genes[-grep("///",ds.genes)]
	}
	genes[[names(esets)[i]]] <- ds.genes
}
SYMBOLS <- lapply(genes, function(i){as.character(i)})
SYMBOLS <- unique(unlist(SYMBOLS))



gene.map <- NULL
for (i in 1:length(esets)){
    eset <-esets[[i]]
    tmp <- rep(0, each = length(SYMBOLS))
    tmp[is.element(SYMBOLS, fData(eset)$EntrezGene.ID)] <- 1
    gene.map <- cbind(gene.map, tmp)
}
colnames(gene.map) <- names(esets)
rownames(gene.map) <- paste("geneid.", SYMBOLS, sep="")

frequency <- NULL
for(i in 1:nrow(gene.map)){
	tmp <- sum(gene.map[i,])
	frequency <- c(frequency, tmp)
}
gene.map <- cbind(gene.map, frequency)
colnames(gene.map)[ncol(gene.map)] <- "frequency"

save(gene.map, file="./esets/summary/breast_genemap.rda")

