library(GEOquery)
library(org.Hs.eg.db)

## create directories
if(!file.exists("./platforms")){
  dir.create("./platforms")
}


# myfn <- file.path("platforms", paste(platforms, ".soft", sep=""))
# if (!file.exists(myfn)) {
#   gpl <- getGEO(platforms, destdir="./platforms")
# } else {
#   ## the platform has already been downloaded
#   gpl <- getGEO(filename=paste("./platforms/", platforms, ".soft", sep=""))
# }

## -------
## saving GPL files
## -------
# feature <-  data.frame(probeset = feature$ID, gene = feature[,"GeneSymbol"], EntrezGene.ID = feature[,"EntrezGene"], row.names = feature$ID)

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
## GPL97
##--------
gpl <- getGEO("GPL97", destdir ="./platforms")
feature <- Table(gpl)
feature <-  data.frame(probeset = feature$ID, gene = feature[,"Gene Symbol"], EntrezGene.ID = feature[,"ENTREZ_GENE_ID"], row.names = feature$ID)
save(feature, file="./annotations/GPL97.rda")

##--------
## GPL1352
##--------
gpl <- getGEO("GPL1352", destdir="./platforms")
feature <- Table(gpl)
feature <- data.frame(probeset = feature$ID, gene = feature[,"Gene Symbol"], EntrezGene.ID = feature[,"ENTREZ_GENE_ID"], row.names= feature$ID)
save(feature, file="./annotations/GPL1352.rda")


##--------
## GPL1390
##--------
gpl <- getGEO("GPL1390", destdir="./platforms")
feature <- Table(gpl)
feature <- data.frame(probeset = feature$ID, gene = feature$GENE_NAME, row.names = feature$ID)
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
rownames(feature) <- feature$probeset <- paste("probe.", feature$probeset, sep="")
save(feature, file="./annotations/GPL1390.rda")

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
  if(length(gene.id) > 1){gene.id <- paste(gene.id[1], gene.id[2], sep="///")}
  EntrezGene.ID <- c(EntrezGene.ID, gene.id)
}
EntrezGene.ID <- unlist(EntrezGene.ID)
feature<- cbind(feature, EntrezGene.ID)
rownames(feature) <- feature$probeset <- paste("probe.", feature$probeset, sep="")
save(feature, file="./annotations/GPL887.rda")

##--------
## GPL885
##--------
gpl <- getGEO("GPL885", destdir="./platforms")
feature <- Table(gpl)
feature <- data.frame(probeset = feature$ID, gene = feature$GENE_SYMBOL, row.names = feature$ID)
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
rownames(feature) <- feature$probeset <- paste("probe.", feature$probeset, sep="")
save(feature, file="./annotations/GPL885.rda")

##--------
## GPL2777
##--------

gpl <- getGEO("GPL2777", destdir="./platforms")
feature <- Table(gpl)
feature <- data.frame(probeset = feature$SPOT_ID, GB_LIST = feature$GB_LIST)

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
save(feature, file="./annotations/GPL2777.rda")


##--------
## GPL2776
##--------

gpl <- getGEO("GPL2776", destdir="./platforms")
feature <- Table(gpl)
feature <- data.frame(probeset = feature$SPOT_ID, GB_LIST = feature$GB_LIST)

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
save(feature, file="./annotations/GPL2776.rda")

##--------
## GPL180
##--------

gpl <- getGEO("GPL180", destdir="./platforms")
feature <- Table(gpl)
feature <- data.frame(probeset = feature$CLONE_ID, gene = feature$GENE_SYM)
feature<- unique(feature)
rownames(feature) <- feature$CLONE_ID
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
save(feature, file="./annotations/GPL180.rda")


##--------
## GPL1261
##--------

gpl <- getGEO("GPL1261", destdir="./platforms")
feature <- Table(gpl)
feature <- data.frame(probeset = feature$ID, gene = feature[,"Gene Symbol"], EntrezGene.ID = feature$ENTREZ_GENE_ID, row.names=feature$ID)
feature$EntrezGene.ID <- gsub(pattern=" /// ", replacement="///", x=feature$EntrezGene.ID)
feature$gene <- gsub(pattern=" /// ", replacement="///", x=feature$gene)
save(feature, file="./annotations/GPL1261.rda")

##--------
## GPL14374
##--------

feature <- getGEO(filename="./platforms/GPL14374.soft")
feature <- Table(feature)
feature <- data.frame(probeset = feature$ID, gene = feature[,"OligoSet_geneSymbol"],  row.names=feature$ID)
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
save(feature, file="./annotations/GPL14374.rda")



##--------
## GPL8300
##--------

feature <- getGEO(filename="./platforms/GPL8300.soft")
feature <- Table(feature)
feature <- data.frame(probeset = feature$ID, gene = feature[,"Gene Symbol"], EntrezGene.ID = feature$ENTREZ_GENE_ID, row.names=feature$ID)
feature$EntrezGene.ID <- gsub(pattern=" /// ", replacement="///", x=feature$EntrezGene.ID)
feature$gene <- gsub(pattern=" /// ", replacement="///", x=feature$gene)

save(feature, file="./annotations/GPL8300.rda")

##--------
## GPL6486
##--------
feature <- getGEO(filename="./platforms/GPL6486.soft")
feature <- Table(feature)
feature <- data.frame(probeset = feature$ID, gene=feature$gene_symbol, row.names=feature$ID)
feature$gene <- gsub(pattern="-", replacement="",x=feature$gene)
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
save(feature, file="./annotations/GPL6486.rda")


##--------
## GPL6106
##--------
feature <- getGEO(filename="./platforms/GPL6106.soft")
feature <- Table(feature)
feature <- data.frame(probeset = feature$ID, gene = feature[,"Symbol"], row.names = feature$ID)
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

feature$probeset <- paste( "probe_",feature$probeset, sep="")
rownames(feature)<- feature$probeset

save(feature, file="./annotations/GPL6106.rda")



##--------
## GPL3883
##--------

feature <- getGEO(filename="./platforms/GPL3883.soft")
feature <- Table(feature)
feature <-  data.frame(probeset = feature$ID, gene = feature[,"GeneSymbol"], EntrezGene.ID = feature[,"EntrezGene"], row.names = feature$ID)
feature$EntrezGene.ID <- gsub(pattern="N/A", replacement="", x=feature$EntrezGene.ID)
feature$gene <- gsub(pattern="N/A", replacement="", x=feature$gene)

feature$EntrezGene.ID <- gsub(pattern="; ", replacement="///", x=feature$EntrezGene.ID)
feature$gene <- gsub(pattern="; ", replacement="///", x=feature$gene)

feature$EntrezGene.ID<- gsub(pattern="\\|", replacement="///", x=feature$EntrezGene.ID)
feature$gene <- gsub(pattern="\\|", replacement="///", x=feature$gene)

save(feature, file="./annotations/GPL3883.rda")




##--------
## GPL4819
##--------

feature <- getGEO(filename="./platforms/GPL4819.soft")
feature <- Table(feature)
feature <- data.frame(probeset = feature$ID, gene= feature[, "Gene Symbol"], row.names = feature$ID)
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
save(feature, file="./annotations/GPL4819.rda")



##################
## Affymetrix HGU
##################

strEsets <- c("CAL","MDA4", "SUPERTAM_HGU133A", "SUPERTAM_HGU133PLUS2", "UNT", "VDX")

load(paste("./annotations/", strEsets[1],"_annot.rda", sep=""))
tmp <- data.frame(probeset = annot$probe, gene = annot$Gene.symbol, EntrezGene.ID = annot$Gene.ID)

for(i in 2:length(strEsets)){
  load(paste("./annotations/", strEsets[i],"_annot.rda", sep=""))
  tmp <- rbind(tmp,data.frame(probeset = annot$probe, gene = annot$Gene.symbol, EntrezGene.ID = annot$Gene.ID))
}
tmp <- unique(tmp)
rownames(tmp) <- tmp$probeset
feature <- tmp
save(feature, file="./annotations/Affymetrix HGU.rda")


##################
## Affymetrix HGU95
##################
StrEsets <- "KOO"
load(paste("./annotations/", strEsets[1],"_annot.rda", sep=""))
feature <- data.frame(probeset = annot$probe, gene = annot$Gene.symbol, EntrezGene.ID = annot$Gene.ID)
rownames(feature) <- feature$probeset
save(feature, file = "./annotations/Affymetrix HGU95.rda")

##################
## Illumina
##################

strEsets <- "HLP"
load(paste("./annotations/", strEsets,"_annot.rda", sep=""))
feature <- data.frame(probeset= annot$probe, gene =annot$symbol, EntrezGene.ID =annot$EntrezGene.ID)
rownames(feature) <- feature$probeset
x <- toTable(org.Hs.egSYMBOL)
x <- x[is.element(x$gene_id,feature$EntrezGene.ID),]
feature <- feature[is.element(feature$EntrezGene.ID, x$gene_id),]

gene <- list()
for(i in 1:nrow(feature)){
    print(i)
   gene <- c(gene, x[x$gene_id == feature$EntrezGene.ID[i], "symbol"])
}
gene <- unlist(gene)
feature$gene <- gene
save(feature,file="./annotations/Illumina.rda")

##################
## Agilent
##################
strEsets <- "NKI"
load(paste("./annotations/", strEsets,"_annot.rda", sep=""))
feature <- data.frame(probeset = annot$probe, gene=annot$NCBI.gene.symbol, EntrezGene.ID = annot$EntrezGene.ID)
rownames(feature) <- feature$probeset
save(feature,file="./annotations/Agilent.rda")


##################
## In-house cDNA
##################
strEsets <- c("NCI", "UCSF")
load(paste("./annotations/", strEsets[1],"_annot.rda", sep=""))
feature <- data.frame(probeset = annot$probe, gene=annot$NCBI.gene.symbol, EntrezGene.ID = annot$EntrezGene.ID)
load(paste("./annotations/", strEsets[2],"_annot.rda", sep=""))
feature <- rbind(feature, data.frame(probeset = annot$probe, gene=annot$Name, EntrezGene.ID = annot$EntrezGene.ID))

rownames(feature) <- feature$probeset
save(feature, file="./annotations/In-house cDNA.rda")


##################
## METABRIC
##################
load(paste("./annotations/", "METABRIC","_annot.rda", sep=""))
annot <- data.frame(annot)
feature <- data.frame(probeset = annot$Probe_Id, gene = annot$Symbol, EntrezGene.ID = annot$EntrezID)
save(feature, file="./annotations/METABRIC.rda")
