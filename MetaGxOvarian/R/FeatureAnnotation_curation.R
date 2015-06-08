#### save soft files for annotations ###
datasets <- read.csv("./datasets.csv")
platforms <- grep("GPL", datasets$Platform,value=T)

if(!file.exists("./annotations")){
  dir.create("./annotations")
}
if(!file.exists("./platforms")){
	dir.create("./platforms")
}



library(org.Hs.eg.db)
library(GEOquery)

##--------
## GPL570
##--------
gpl <- getGEO("GPL570", destdir ="./platforms")
feature <- Table(gpl)
feature <-  data.frame(probeset = feature$ID, gene = feature[,"Gene Symbol"], EntrezGene.ID = feature[,"ENTREZ_GENE_ID"], row.names = feature$ID)
save(feature, file="./annotations/GPL570.rda")

##--------
## GPL96
##--------
gpl <- getGEO("GPL96", destdir ="./platforms")
feature <- Table(gpl)
feature <-  data.frame(probeset = feature$ID, gene = feature[,"Gene Symbol"], EntrezGene.ID = feature[,"ENTREZ_GENE_ID"], row.names = feature$ID)
save(feature, file="./annotations/GPL96.rda")

##--------
## GPL887
##--------
gpl <- getGEO("GPL887", destdir="./platforms")
feature <- Table(gpl)
feature <- data.frame(probeset = feature$ID, gene = feature$GENE_SYMBOL, row.names = feature$ID)
x <- org.Hs.egSYMBOL2EG
x <- toTable(x)
feature <- feature[is.element(feature$gene, x$symbol),]
EntrezGene.ID <- list()
for(i in 1:nrow(feature)){
  gene.id <- x[x$symbol == feature$gene[i], "gene_id"]
  if(length(gene.id) > 1){gene.id <- paste(gene.id, collapse="///")}
  EntrezGene.ID <- c(EntrezGene.ID, gene.id)
}
EntrezGene.ID <- unlist(EntrezGene.ID)
feature<- cbind(feature, EntrezGene.ID)
rownames(feature) <- feature$probeset 
save(feature, file="./annotations/GPL887.rda")

##--------
## GPL6104
##--------
gpl <- getGEO("GPL6104", destdir="./platforms")
feature <- Table(gpl)
feature <- data.frame(probeset=feature$ID, gene=feature$Symbol, EntrezGene.ID=feature$Entrez_Gene_ID, row.names = feature$ID)
save(feature, file="./annotations/GPL6104.rda")

##--------
## GPL5886
##--------
gpl <- getGEO("GPL5886", destdir="./platforms")
feature <- Table(gpl)
feature <- data.frame(probeset=feature$ID, gene=feature$OligoSet_geneSymbol, row.names = feature$ID)
x <- org.Hs.egSYMBOL2EG
x <- toTable(x)
feature <- feature[is.element(feature$gene, x$symbol),]
EntrezGene.ID <- list()
for(i in 1:nrow(feature)){
  gene.id <- x[x$symbol == feature$gene[i], "gene_id"]
  if(length(gene.id) > 1){gene.id <- paste(gene.id, collapse="///")}
  EntrezGene.ID <- c(EntrezGene.ID, gene.id)
}
EntrezGene.ID <- unlist(EntrezGene.ID)
feature<- cbind(feature, EntrezGene.ID)
save(feature, file="./annotations/GPL5886.rda")

##--------
## GPL7759
##--------
gpl <- getGEO("GPL7759", destdir="./platforms")
feature <- Table(gpl)
feature <- data.frame(probeset=feature$ID, gene=feature[,"Gene Symbol"], row.names = feature$ID)
x <- org.Hs.egSYMBOL2EG
x <- toTable(x)
feature <- feature[is.element(feature$gene, x$symbol),]
EntrezGene.ID <- list()
for(i in 1:nrow(feature)){
  gene.id <- x[x$symbol == feature$gene[i], "gene_id"]
  if(length(gene.id) > 1){gene.id <- paste(gene.id, collapse="///")}
  EntrezGene.ID <- c(EntrezGene.ID, gene.id)
}
EntrezGene.ID <- unlist(EntrezGene.ID)
feature<- cbind(feature, EntrezGene.ID)
save(feature, file="./annotations/GPL7759.rda")


##--------
## GPL6480
##--------
gpl <- getGEO("GPL6480", destdir="./platforms")
feature <- Table(gpl)
feature <- data.frame(probeset=feature$ID, gene=feature$GENE_SYMBOL, EntrezGene.ID=feature$GENE)
feature <- feature[-which(is.na(feature$probeset)),]
rownames(feature) <- feature$probeset
save(feature, file="./annotations/GPL6480.rda")

##--------
## GPL8300
##--------
gpl <- getGEO("GPL8300", destdir="./platforms")
feature <- Table(gpl)
feature <- data.frame(probeset=feature$ID, gene=feature[,"Gene Symbol"], EntrezGene.ID=feature[, "ENTREZ_GENE_ID"], row.names = feature$ID)
save(feature, file="./annotations/GPL8300.rda")


# ##--------
# ## GPL2005
# ##--------
# gpl <- gpl <- getGEO("GPL2005", destdir="./platforms")
# feature <- Table(gpl)
# feature <- data.frame(probeset=feature$ID, gene=feature[,"Gene Symbol"], EntrezGene.ID=feature[, "ENTREZ_GENE_ID"])
# save(feature, file="./annotations/GPL8300.rda")



# ##--------
# ## GPL6801
# ##--------

# gpl <- getGEO("GPL6801", destdir="./platforms")
# feature <- Table(gpl)
# feature <- data.frame(probeset=feature$ID, gene=feature[,"Gene Symbol"], EntrezGene.ID=feature[, "ENTREZ_GENE_ID"])
# save(feature, file="./annotations/GPL6801.rda")


##--------
## GPL13728
##--------

gpl <- getGEO("GPL13728", destdir="./platforms")
feature <- Table(gpl)
feature <- data.frame(probeset=feature$ID, gene=feature[,"ORF_LIST"], row.names = feature$ID)
x <- org.Hs.egSYMBOL2EG
x <- toTable(x)
feature <- feature[is.element(feature$gene, x$symbol),]
EntrezGene.ID <- list()
for(i in 1:nrow(feature)){
  gene.id <- x[x$symbol == feature$gene[i], "gene_id"]
  if(length(gene.id) > 1){gene.id <- paste(gene.id, collapse="///")}
  EntrezGene.ID <- c(EntrezGene.ID, gene.id)
}
EntrezGene.ID <- unlist(EntrezGene.ID)
feature<- cbind(feature, EntrezGene.ID)
save(feature, file="./annotations/GPL13728.rda")


##--------
## GPL80
##--------

gpl <- getGEO("GPL80", destdir="./platforms")
feature <- Table(gpl)
feature <- data.frame(probeset=feature$ID, gene=feature[,"Gene Symbol"], EntrezGene.ID=feature[, "ENTREZ_GENE_ID"], row.names = feature$ID)
save(feature, file="./annotations/GPL80.rda")


##--------
## GPL4685
##--------

gpl <- getGEO("GPL4685", destdir="./platforms")
feature <- Table(gpl)
feature <- data.frame(probeset=feature$ID, gene=feature[,"Gene Symbol"], EntrezGene.ID=feature[, "Entrez Gene"], row.names = feature$ID)
save(feature, file="./annotations/GPL4685.rda")

##--------
## GPL9530
##--------

gpl <- getGEO("GPL9530", destdir="./platforms")
feature <- Table(gpl)
feature <- data.frame(probeset=feature$ID, gene=feature[,"GENE NAME"], row.names = feature$ID)
x <- org.Hs.egSYMBOL2EG
x <- toTable(x)
feature <- feature[is.element(feature$gene, x$symbol),]
EntrezGene.ID <- list()
for(i in 1:nrow(feature)){
  gene.id <- x[x$symbol == feature$gene[i], "gene_id"]
  if(length(gene.id) > 1){gene.id <- paste(gene.id, collapse="///")}
  EntrezGene.ID <- c(EntrezGene.ID, gene.id)
}
EntrezGene.ID <- unlist(EntrezGene.ID)
feature<- cbind(feature, EntrezGene.ID)
save(feature, file="./annotations/GPL9530.rda")



##################### Extra GPLs for extension (other potential esets)


##--------
## GPL7264
##--------

gpl <- getGEO("GPL7264", destdir ="./platforms")
feature <- Table(gpl)
feature <- data.frame(probeset=feature$ID, gene=feature$GENE_SYMBOL, EntrezGene.ID=feature$GENE)
feature <- feature[-which(is.na(feature$probeset)),]
rownames(feature) <- feature$probeset
save(feature, file="./annotations/GPL7264.rda")


##--------
## GPL2986
##--------

gpl <- getGEO("GPL2986", destdir ="./platforms")
feature <- Table(gpl)
feature <- data.frame(probeset=feature$ID, gene=feature[,"Gene Symbol"], EntrezGene.ID =feature$GENE)
rownames(feature) <- feature$probeset
save(feature, file="./annotations/GPL2986.rda")

## -------
## GPL5689
##---------

gpl <- getGEO("GPL5689", destdir="./platforms")
feature <- Table(gpl)
feature <- data.frame(probeset = feature$ID, GB_LIST=feature$GB_ACC)
x<- toTable(org.Hs.egACCNUM2EG)
common <- unlist(strsplit(x=as.character(feature$GB_LIST), split=","))
x <-x[is.element(x$accession,common),]
EntrezGene.ID <- list()
EntrezGene.ID <- rep(NA, nrow(feature))

for(i in 1:nrow(x)){
  EntrezGene.ID[grep(x=feature$GB_LIST, pattern=as.character(x$accession[i]))] <- x$gene_id[i]  
}
feature <- cbind(feature, EntrezGene.ID)
feature <- na.omit(feature)

x <- toTable(org.Hs.egSYMBOL)
x <- x[is.element(x$gene_id,feature$EntrezGene.ID),]

gene <- list()
for(i in 1:nrow(feature)){
  gene <- c(gene, x[x$gene_id == feature$EntrezGene.ID[i], "symbol"])
}
gene <- unlist(gene)
feature <- data.frame(probeset = feature$probeset, gene = gene, EntrezGene.ID = feature$EntrezGene.ID)
save(feature, file="./annotations/GPL5689.rda")

## -------
## Illumina
## -------

library(org.Hs.eg.db)

load("./FULLVdata/TCGA.RNASeqV2_eset.rda")
annot <- fData(TCGA.RNASeqV2_eset)
save(annot, file="./annotations/TCGA.RNASeqV2_annot.rda")
EntrezGene.ID <- gsub(pattern = ".+[[:punct:]]", replacement="", x=annot$probeset)
x <- toTable(org.Hs.egSYMBOL2EG)
x <- x[is.element(x$gene_id, EntrezGene.ID),]
probeset <- NULL
for(i in 1:nrow(x)){
  probeset <- c(probeset, annot$probeset[which(EntrezGene.ID == x$gene_id[i])])
}
feature <- data.frame(probeset=probeset, gene=x$symbol, EntrezGene.ID=x$gene_id, row.names=probeset)
save(feature, file="./annotations/Illumina.rda")




