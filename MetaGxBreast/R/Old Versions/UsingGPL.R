###################################
## Natchar Ratanasirigulchai
## April 13, 2015
###################################


## ----------------------------------------
## Create expression sets from each dataset
## ----------------------------------------

library(Biobase)
library(GEOquery)  
library(WGCNA) 
datasets <- read.csv("datasets.csv")
dataset.names <- datasets$Dataset
eset <- NULL

## create directories
if(!file.exists("./platforms")){
  dir.create("./platforms")
}

if(!file.exists("./esets")){
  dir.create("./esets")
  dir.create("./esets/mapped_esets")
}

## map each eset using GPL or existing annotation

for(i in (1:length(dataset.names))){
  message(i)
  dataset.name <- dataset.names[i]
  expr <- read.table(paste("./data/", dataset.name, "/PROCESSED/DEFAULT/", dataset.name, "_default_curated_exprs.txt", sep=""))
  if(dataset.name == "METABRIC"){
    expr <- t(expr)
  }
  platforms <- as.character(datasets[datasets$Dataset == dataset.name,"Platform"])
  platforms <- strsplit(platforms, split="///", fixed=TRUE)[[1]]
  
  feature <- NULL
  
  if (length(platforms) < 1 | length(platforms) >1){
    ## No specific platform identified, get from annot
    load(paste("./annotations/", dataset.name,"_annot.rda", sep=""))
    
    ## following conditionals because difference in labels
    if(dataset.name == "HLP"| dataset.name =="UNC4"){
      feature <- annot[, c("probe", "symbol", "EntrezGene.ID")]
    } else if(dataset.name == "NCI" | dataset.name == "NKI" | dataset.name=="STNO2") {
      feature <- annot[, c("probe", "HUGO.gene.symbol", "EntrezGene.ID")]
    } else if(dataset.name == "UCSF"){
      feature <- annot[,c("probe", "Name2", "EntrezGene.ID")]
    } else if(dataset.name =="METABRIC"){
      feature <- annot[, c("Probe_Id", "Symbol", "EntrezID")]
      feature <- data.frame(feature)
    } else {
      feature <- annot[, c("probe", "Gene.symbol", "Gene.ID")]
    }
  } else {
    myfn <- file.path("platforms", paste(platforms, ".soft", sep=""))
    if (!file.exists(myfn)) {
      gpl <- getGEO(platforms, destdir="./platforms")
    } else {
      ## the platform has already been downloaded
      gpl <- getGEO(filename=paste("./platforms/", platforms, ".soft", sep=""))
    }
    if(!is.element("ENTREZ_GENE_ID", colnames(Table(gpl)))){ 
    ## GPL does not have entrez genes at all, look at annot from journal
      load(paste("./annotations/", dataset.name,"_annot.rda", sep=""))
      
      ## following conditionals because differences in labels
      if(platforms == "GPL4819" | platforms == "GPL14374"){
        feature <- annot[, c("probe", "HUGO.gene.symbol", "EntrezGene.ID")]
        row.names(feature) <- feature$probe
      } else if(platforms == "GPL6106"){
        feature <- annot[, c("probe", "Symbol", "EntrezGene.ID")]
        row.names(feature) <- feature$probe
      } else if(platforms == "GPL6486"){
        feature <- annot[,c("probe", "gene_symbol", "EntrezGene.ID")]
        row.names(feature) <- feature$probe
      } else{
        feature <- annot[, c("probe", "Gene.symbol", "Gene.ID")]
        row.names(feature) <- feature$probe
      }
    } else {
      feature <- Table(gpl)[Table(gpl)$ID == intersect(rownames(expr),Table(gpl)$ID),c("ID", "Gene Symbol", "ENTREZ_GENE_ID")]
      row.names(feature) <- feature$ID
    }
  }  
  message("got feature of ", dataset.name)
  colnames(feature) <- c("PROBE", "SYMBOL", "ENTREZID")
  feature[feature==""|feature==" "] <- NA
  feature <- na.omit(feature)
#   rownames(feature) <- feature$PROBE
  ## x is entrez ID's with exprs(eset)
  expr <- expr[intersect(row.names(expr),rownames(feature)),]
  x <- cbind(feature[,"ENTREZID"],expr)
  
 
  ## probe gene mapping options : MaxMean, maxRowVariance, Average
  Rprof("profile_mm_wgcna.txt")
  collapsed.wgcna.mrv <- collapseRows(x[,-1], rowID=rownames(x), rowGroup=x[,1], method="maxRowVariance")
  Rprof(NULL)
  
  best_probe <- collapsed.wgcna.mrv[[3]]
  feature <- cbind(feature, best_probe)
  message("collapsed rows of ", dataset.name)
  ## set best_probe for ambiguous entrez to false
  feature[grep("///", feature$ENTREZID), "best_probe"] <- FALSE

  ### create the expressionsets
  feature.df <- AnnotatedDataFrame(feature)
  pData <- read.table(paste("./curation/breast/curated/",dataset.names[i], "_curated.txt", sep="" ), header=TRUE)
  pData <- AnnotatedDataFrame(pData)
  rownames(pData) <- colnames(expr)

  #!!! expression set MUST BE MATRIX
  # !!! Colnames of exprs must match rownames of pData (no transpose of exprs!)
  eset <- ExpressionSet(assayData=(as.matrix(expr)), featureData=feature.df, phenoData=pData)
  assign(as.character(dataset.name), eset)
  save(list=as.character(dataset.name), file =paste("esets/mapped_esets/",dataset.names[i], "_eset.rda", sep=""))
  
}

