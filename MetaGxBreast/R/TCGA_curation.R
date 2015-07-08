library(org.Hs.eg.db)
load("./preprocessing/TCGAcurations.RData")
expr <- TCGA_exprs

## if symbol is duplicated, rename to be __.2 for expr
duplicates.names <- colnames(expr)[duplicated(colnames(expr))]

## save symbols for fData annotations
gene <- colnames(expr)

tmp <- colnames(expr)
for( i in 1:length(duplicates.names)){
	indices <- which(duplicates.names[i]== colnames(expr))
	for(n in 2:length(indices)){
		tmp[indices[n]] <- paste(tmp[indices[n]], ".", as.character(n), sep="")
	}
}

colnames(expr) <- tmp
rm(tmp)

##extract rows of expr of samples in available pData
pData <- read.table("./curation/breast/curated/TCGA_curated.txt", header=T)

expr <- expr[is.element(rownames(expr),pData$sample_name),]
pData <- pData[is.element(pData$sample_name, rownames(expr)),]

##create featureData annotations

feature <- data.frame(probeset = colnames(expr), gene=gene, row.names=colnames(expr))
x <- org.Hs.egSYMBOL2EG
x <- toTable(x)
feature <- feature[is.element(feature$gene, x$symbol),]
EntrezGene.ID <- list()
for(i in 1:nrow(feature)){
  gene.id <- x[x$symbol == feature$gene[i], "gene_id"]
  if(length(gene.id) > 1){gene.id <- paste(gene.id[1], gene.id[2], sep="///")}
  EntrezGene.ID <- c(EntrezGene.ID, gene.id)
}
EntrezGene.ID <- unlist(EntrezGene.ID)
feature<- cbind(feature, EntrezGene.ID)
rownames(feature) <- feature$probeset 
expr <- t(expr)
write.table(expr, file="./data/TCGA/PROCESSED/DEFAULT/TCGA_default_curated_exprs.txt", sep="\t")
write.table(pData, file= "./curation/breast/curated/TCGA_curated.txt", row.names=FALSE, sep="\t")
save(feature, file ="./annotations/TCGA.rda")
