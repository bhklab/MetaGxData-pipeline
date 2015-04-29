###################################
## Natchar Ratanasirigulchai
## April 14, 2015
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
  dir.create("./esets/mapped_esets2")
}
nbcore <- 8
availcore <- parallel::detectCores()
if (is.null(nbcore) || nbcore > availcore) { nbcore <- availcore }
options("mc.cores"=nbcore)
## map each eset using GPL or existing annotation

parallel::mclapply((1:length(dataset.names)), function(i){
  message(i)
  dataset.name <- dataset.names[i]
  expr <- read.table(paste("./data/", dataset.name, "/PROCESSED/DEFAULT/", dataset.name, "_default_curated_exprs.txt", sep=""))
  if(dataset.name == "METABRIC"){
    expr <- t(expr)
  }
  
  platforms <- as.character(datasets[datasets$Dataset == dataset.name,"Platform"])
  platforms <- unlist(strsplit(platforms, split="///", fixed=TRUE))
  if(length(platforms) < 1){
    warning("There is no platform listed. Please provide the platform annotations in the annotations directory and update the dataset.csv file.")
    next
  }
  
  load(paste("./annotations/", platforms[1], ".rda", sep=""))
  if(length(platforms)>1){
    for(n in 1:length(platforms)){
      tmp <- feature
      load(paste("./annotations/", platforms[n], ".rda", sep=""))
      feature <- rbind(tmp, feature)
      feature <- unique(feature)
      rownames(feature) <- feature$PROBE
    }
  }
  
  message("got feature of ", dataset.name)
  #   colnames(feature) <- c("PROBE", "SYMBOL", "ENTREZID")
  feature[feature==""|feature==" "] <- NA
  feature <- na.omit(feature)
  #   rownames(feature) <- feature$PROBE
  ## x is entrez ID's with exprs(eset)
  expr <- expr[is.element(row.names(expr),rownames(feature)),]
  feature <- feature[is.element(rownames(feature),rownames(expr)),]
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
  feature <- feature[match(rownames(expr),rownames(feature)),]
  ### create the expressionsets
  feature.df <- AnnotatedDataFrame(feature)
  pData <- read.table(paste("./curation/breast/curated/",dataset.names[i], "_curated.txt", sep="" ), header=TRUE)
  pData <- AnnotatedDataFrame(pData)
  rownames(pData) <- colnames(expr)
  
  #!!! expression set MUST BE MATRIX
  # !!! Colnames of exprs must match rownames of pData (no transpose of exprs!)
  eset <- ExpressionSet(assayData=(as.matrix(expr)), featureData=feature.df, phenoData=pData)
  assign(as.character(dataset.name), eset)
  save(list=as.character(dataset.name), file =paste("esets/mapped_esets2/",dataset.names[i], "_eset.rda", sep=""))
  
})
