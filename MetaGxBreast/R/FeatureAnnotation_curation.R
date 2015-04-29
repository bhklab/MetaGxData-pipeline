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
# feature <-  data.frame(PROBE = feature$ID, SYMBOL = feature[,"GeneSymbol"], ENTREZID = feature[,"EntrezGene"], row.names = feature$ID)

##--------
## GPL1390
##--------
gpl <- getGEO("GPL1390", destdir="./platforms")
feature <- Table(gpl)
feature <- data.frame(PROBE = feature$ID, SYMBOL = feature$GENE_NAME, row.names = feature$ID)
x <- org.Hs.egSYMBOL2EG
x <- toTable(x)
feature <- feature[is.element(feature$SYMBOL, x$symbol),]
ENTREZID <- list()
for(i in 1:nrow(feature)){
  gene.id <- x[x$symbol == feature$SYMBOL[i], "gene_id"]
  if(length(gene.id) > 1){gene.id <- paste(gene.id[1], gene.id[2], sep="///")}
  ENTREZID <- c(ENTREZID, gene.id)
}
ENTREZID <- unlist(ENTREZID)
feature<- cbind(feature, ENTREZID)
rownames(feature) <- feature$PROBE <- paste("probe.", feature$PROBE, sep="")
save(feature, file="./annotations/GPL1390.rda")

##--------
## GPL887
##--------
gpl <- getGEO("GPL887", destdir="./platforms")
feature <- Table(gpl)
feature <- data.frame(PROBE = feature$ID, SYMBOL = feature$GENE_SYMBOL, row.names = feature$ID)
x <- org.Hs.egSYMBOL2EG
x <- toTable(x)
feature <- feature[is.element(feature$SYMBOL, x$symbol),]
ENTREZID <- list()
for(i in 1:nrow(feature)){
  gene.id <- x[x$symbol == feature$SYMBOL[i], "gene_id"]
  if(length(gene.id) > 1){gene.id <- paste(gene.id[1], gene.id[2], sep="///")}
  ENTREZID <- c(ENTREZID, gene.id)
}
ENTREZID <- unlist(ENTREZID)
feature<- cbind(feature, ENTREZID)
rownames(feature) <- feature$PROBE <- paste("probe.", feature$PROBE, sep="")
save(feature, file="./annotations/GPL887.rda")

##--------
## GPL885
##--------
gpl <- getGEO("GPL885", destdir="./platforms")
feature <- Table(gpl)
feature <- data.frame(PROBE = feature$ID, SYMBOL = feature$GENE_SYMBOL, row.names = feature$ID)
x <- org.Hs.egSYMBOL2EG
x <- toTable(x)
feature <- feature[is.element(feature$SYMBOL, x$symbol),]
ENTREZID <- list()
for(i in 1:nrow(feature)){
  gene.id <- x[x$symbol == feature$SYMBOL[i], "gene_id"]
  if(length(gene.id) > 1){gene.id <- paste(gene.id[1], gene.id[2], sep="///")}
  ENTREZID <- c(ENTREZID, gene.id)
}
ENTREZID <- unlist(ENTREZID)
feature<- cbind(feature, ENTREZID)
rownames(feature) <- feature$PROBE <- paste("probe.", feature$PROBE, sep="")
save(feature, file="./annotations/GPL885.rda")

##--------
## GPL2777
##--------

gpl <- getGEO("GPL2777", destdir="./platforms")
feature <- Table(gpl)
feature <- data.frame(PROBE = feature$SPOT_ID, GB_LIST = feature$GB_LIST)

x<- toTable(org.Hs.egACCNUM2EG)
common <- unlist(strsplit(x=as.character(feature$GB_LIST), split=","))
x <-x[is.element(x$accession,common),]
ENTREZID <- list()
ENTREZID <- rep(NA, nrow(feature))
for(i in 1:nrow(x)){
  ENTREZID[grep(x=feature$GB_LIST, pattern=as.character(x$accession[i]))] <- x$gene_id[i]  
}
feature <- cbind(feature, ENTREZID)
feature <- na.omit(feature)

x <- toTable(org.Hs.egSYMBOL)
x <- x[is.element(x$gene_id,feature$ENTREZID),]

SYMBOL <- list()
for(i in 1:nrow(feature)){
  SYMBOL <- c(SYMBOL, x[x$gene_id == feature$ENTREZID[i], "symbol"])
}
SYMBOL <- unlist(SYMBOL)
feature <- data.frame(PROBE = feature$PROBE, SYMBOL = SYMBOL, ENTREZID = feature$ENTREZID)
save(feature, file="./annotations/GPL2777.rda")


##--------
## GPL2776
##--------

gpl <- getGEO("GPL2776", destdir="./platforms")
feature <- Table(gpl)
feature <- data.frame(PROBE = feature$SPOT_ID, GB_LIST = feature$GB_LIST)

x<- toTable(org.Hs.egACCNUM2EG)
common <- unlist(strsplit(x=as.character(feature$GB_LIST), split=","))
x <-x[is.element(x$accession,common),]
ENTREZID <- list()
ENTREZID <- rep(NA, nrow(feature))
for(i in 1:nrow(x)){
  ENTREZID[grep(x=feature$GB_LIST, pattern=as.character(x$accession[i]))] <- x$gene_id[i]  
}
feature <- cbind(feature, ENTREZID)
feature <- na.omit(feature)

x <- toTable(org.Hs.egSYMBOL)
x <- x[is.element(x$gene_id,feature$ENTREZID),]

SYMBOL <- list()
for(i in 1:nrow(feature)){
  SYMBOL <- c(SYMBOL, x[x$gene_id == feature$ENTREZID[i], "symbol"])
}
SYMBOL <- unlist(SYMBOL)
feature <- data.frame(PROBE = feature$PROBE, SYMBOL = SYMBOL, ENTREZID = feature$ENTREZID)
save(feature, file="./annotations/GPL2776.rda")

##--------
## GPL180
##--------

gpl <- getGEO("GPL180", destdir="./platforms")
feature <- Table(gpl)
feature <- data.frame(PROBE = feature$CLONE_ID, SYMBOL = feature$GENE_SYM)
feature<- unique(feature)
rownames(feature) <- feature$CLONE_ID
x <- org.Hs.egSYMBOL2EG
x <- toTable(x)
feature <- feature[is.element(feature$SYMBOL, x$symbol),]
ENTREZID <- list()
for(i in 1:nrow(feature)){
  gene.id <- x[x$symbol == feature$SYMBOL[i], "gene_id"]
  if(length(gene.id) > 1){gene.id <- paste(gene.id[1], gene.id[2], sep="///")}
  ENTREZID <- c(ENTREZID, gene.id)
}
ENTREZID <- unlist(ENTREZID)
feature<- cbind(feature, ENTREZID)
rownames(feature) <- feature$PROBE
save(feature, file="./annotations/GPL180.rda")


##--------
## GPL1261
##--------

gpl <- getGEO("GPL1261", destdir="./platforms")
feature <- Table(gpl)
feature <- data.frame(PROBE = feature$ID, SYMBOL = feature[,"Gene Symbol"], ENTREZID = feature$ENTREZ_GENE_ID, row.names=feature$ID)
feature$ENTREZID <- gsub(pattern=" /// ", replacement="///", x=feature$ENTREZID)
feature$SYMBOL <- gsub(pattern=" /// ", replacement="///", x=feature$SYMBOL)
save(feature, file="./annotations/GPL1261.rda")

##--------
## GPL14374
##--------

feature <- getGEO(filename="./platforms/GPL14374.soft")
feature <- Table(feature)
feature <- data.frame(PROBE = feature$ID, SYMBOL = feature[,"OligoSet_geneSymbol"],  row.names=feature$ID)
x <- org.Hs.egSYMBOL2EG
x <- toTable(x)
feature <- feature[is.element(feature$SYMBOL, x$symbol),]
ENTREZID <- list()
for(i in 1:nrow(feature)){
  gene.id <- x[x$symbol == feature$SYMBOL[i], "gene_id"]
  if(length(gene.id) > 1){gene.id <- paste(gene.id[1], gene.id[2], sep="///")}
  ENTREZID <- c(ENTREZID, gene.id)
}
ENTREZID <- unlist(ENTREZID)
feature<- cbind(feature, ENTREZID)
save(feature, file="./annotations/GPL14374.rda")



##--------
## GPL8300
##--------

feature <- getGEO(filename="./platforms/GPL8300.soft")
feature <- Table(feature)
feature <- data.frame(PROBE = feature$ID, SYMBOL = feature[,"Gene Symbol"], ENTREZID = feature$ENTREZ_GENE_ID, row.names=feature$ID)
feature$ENTREZID <- gsub(pattern=" /// ", replacement="///", x=feature$ENTREZID)
feature$SYMBOL <- gsub(pattern=" /// ", replacement="///", x=feature$SYMBOL)

save(feature, file="./annotations/GPL8300.rda")

##--------
## GPL6486
##--------
feature <- getGEO(filename="./platforms/GPL6486.soft")
feature <- Table(feature)
feature <- data.frame(PROBE = feature$ID, SYMBOL=feature$gene_symbol, row.names=feature$ID)
feature$SYMBOL <- gsub(pattern="-", replacement="",x=feature$SYMBOL)
x <- org.Hs.egSYMBOL2EG
x <- toTable(x)
feature <- feature[is.element(feature$SYMBOL, x$symbol),]
ENTREZID <- list()
for(i in 1:nrow(feature)){
  gene.id <- x[x$symbol == feature$SYMBOL[i], "gene_id"]
  if(length(gene.id) > 1){gene.id <- paste(gene.id[1], gene.id[2], sep="///")}
  ENTREZID <- c(ENTREZID, gene.id)
}
ENTREZID <- unlist(ENTREZID)
feature<- cbind(feature, ENTREZID)
save(feature, file="./annotations/GPL6486.rda")


##--------
## GPL6106
##--------
feature <- getGEO(filename="./platforms/GPL6106.soft")
feature <- Table(feature)
feature <- data.frame(PROBE = feature$ID, SYMBOL = feature[,"Symbol"], row.names = feature$ID)
x <- org.Hs.egSYMBOL2EG
x <- toTable(x)
feature <- feature[is.element(feature$SYMBOL, x$symbol),]
ENTREZID <- list()
for(i in 1:nrow(feature)){
  gene.id <- x[x$symbol == feature$SYMBOL[i], "gene_id"]
  if(length(gene.id) > 1){gene.id <- paste(gene.id[1], gene.id[2], sep="///")}
  ENTREZID <- c(ENTREZID, gene.id)
}
ENTREZID <- unlist(ENTREZID)
feature<- cbind(feature, ENTREZID)

feature$PROBE <- paste( "probe_",feature$PROBE, sep="")
rownames(feature)<- feature$PROBE

save(feature, file="./annotations/GPL6106.rda")



##--------
## GPL3883
##--------

feature <- getGEO(filename="./platforms/GPL3883.soft")
feature <- Table(feature)
feature <-  data.frame(PROBE = feature$ID, SYMBOL = feature[,"GeneSymbol"], ENTREZID = feature[,"EntrezGene"], row.names = feature$ID)
feature$ENTREZID <- gsub(pattern="N/A", replacement="", x=feature$ENTREZID)
feature$SYMBOL <- gsub(pattern="N/A", replacement="", x=feature$SYMBOL)

feature$ENTREZID <- gsub(pattern="; ", replacement="///", x=feature$ENTREZID)
feature$SYMBOL <- gsub(pattern="; ", replacement="///", x=feature$SYMBOL)

feature$ENTREZID<- gsub(pattern="\\|", replacement="///", x=feature$ENTREZID)
feature$SYMBOL <- gsub(pattern="\\|", replacement="///", x=feature$SYMBOL)

save(feature, file="./annotations/GPL3883.rda")




##--------
## GPL4819
##--------

feature <- getGEO(filename="./platforms/GPL4819.soft")
feature <- Table(feature)
feature <- data.frame(PROBE = feature$ID, SYMBOL= feature[, "Gene Symbol"], row.names = feature$ID)
x <- org.Hs.egSYMBOL2EG
x <- toTable(x)
feature <- feature[is.element(feature$SYMBOL, x$symbol),]
ENTREZID <- list()
for(i in 1:nrow(feature)){
  gene.id <- x[x$symbol == feature$SYMBOL[i], "gene_id"]
  if(length(gene.id) > 1){gene.id <- paste(gene.id[1], gene.id[2], sep="///")}
  ENTREZID <- c(ENTREZID, gene.id)
}
ENTREZID <- unlist(ENTREZID)
feature<- cbind(feature, ENTREZID)
save(feature, file="./annotations/GPL4819.rda")



##################
## Affymetrix HGU
##################

strEsets <- c("CAL","MDA4", "SUPERTAM_HGU133A", "SUPERTAM_HGU133PLUS2", "UNT", "VDX")

load(paste("./annotations/", strEsets[1],"_annot.rda", sep=""))
tmp <- data.frame(PROBE = annot$probe, SYMBOL = annot$Gene.symbol, ENTREZID = annot$Gene.ID)

for(i in 2:length(strEsets)){
  load(paste("./annotations/", strEsets[i],"_annot.rda", sep=""))
  tmp <- rbind(tmp,data.frame(PROBE = annot$probe, SYMBOL = annot$Gene.symbol, ENTREZID = annot$Gene.ID))
}
tmp <- unique(tmp)
rownames(tmp) <- tmp$PROBE
feature <- tmp
save(feature, file="./annotations/Affymetrix HGU.rda")


##################
## Affymetrix HGU95
##################
StrEsets <- "KOO"
load(paste("./annotations/", strEsets[1],"_annot.rda", sep=""))
feature <- data.frame(PROBE = annot$probe, SYMBOL = annot$Gene.symbol, ENTREZID = annot$Gene.ID)
rownames(feature) <- feature$PROBE
save(feature, file = "./annotations/Affymetrix HGU95.rda")

##################
## Illumina
##################

strEsets <- "HLP"
load(paste("./annotations/", strEsets,"_annot.rda", sep=""))
feature <- data.frame(PROBE= annot$probe, SYMBOL =annot$symbol, ENTREZID =annot$EntrezGene.ID)
rownames(feature) <- feature$PROBE
x <- toTable(org.Hs.egSYMBOL)
x <- x[is.element(x$gene_id,feature$ENTREZID),]
feature <- feature[is.element(feature$ENTREZID, x$gene_id),]

SYMBOL <- list()
for(i in 1:nrow(feature)){
    print(i)
   SYMBOL <- c(SYMBOL, x[x$gene_id == feature$ENTREZID[i], "symbol"])
}
SYMBOL <- unlist(SYMBOL)
feature$SYMBOL <- SYMBOL
save(feature,file="./annotations/Illumina.rda")

##################
## Agilent
##################
strEsets <- "NKI"
load(paste("./annotations/", strEsets,"_annot.rda", sep=""))
feature <- data.frame(PROBE = annot$probe, SYMBOL=annot$NCBI.gene.symbol, ENTREZID = annot$EntrezGene.ID)
rownames(feature) <- feature$PROBE
save(feature,file="./annotations/Agilent.rda")


##################
## In-house cDNA
##################
strEsets <- c("NCI", "UCSF")
load(paste("./annotations/", strEsets[1],"_annot.rda", sep=""))
feature <- data.frame(PROBE = annot$probe, SYMBOL=annot$NCBI.gene.symbol, ENTREZID = annot$EntrezGene.ID)
load(paste("./annotations/", strEsets[2],"_annot.rda", sep=""))
feature <- rbind(feature, data.frame(PROBE = annot$probe, SYMBOL=annot$Name, ENTREZID = annot$EntrezGene.ID))

rownames(feature) <- feature$PROBE
save(feature, file="./annotations/In-house cDNA.rda")


##################
## METABRIC
##################
load(paste("./annotations/", "METABRIC","_annot.rda", sep=""))
annot <- data.frame(annot)
feature <- data.frame(PROBE = annot$Probe_Id, SYMBOL = annot$Symbol, ENTREZID = annot$EntrezID)
save(feature, file="./annotations/METABRIC.rda")
