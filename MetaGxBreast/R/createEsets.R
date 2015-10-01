###################################
## Natchar Ratanasirigulchai
## April 14, 2015
###################################
setwd("/mnt/work1/users/bhklab/Data/MetaGxData/MetaGxBreast") #Ensure working directory is correct!!

mapping.method <- "maxRowVariance"
#mapping.method <- "maxMean"
mapping.group <- "EntrezGene.ID"
#mapping.group <- "gene"

## ----------------------------------------
## Create expression sets from each dataset
## ----------------------------------------
library(Biobase)
library(GEOquery)  
library(WGCNA) 
datasets <- read.csv("datasets.csv", stringsAsFactors = FALSE)
dataset.names <- datasets$Dataset
eset <- NULL
eData <- read.csv("./etc/Breast Cancer Database Literature Review.csv", stringsAsFactors=FALSE)
if(!file.exists("./esets")){
  dir.create("./esets")
  dir.create("./esets/mapped_esets2")
}

## map each eset using GPL or existing annotation
if(Sys.getenv("SGE_TASK_ID") == "") {
	vals.to.run <- 1:length(dataset.names)
} else {
	vals.to.run <- as.integer(Sys.getenv("SGE_TASK_ID"))
}

for(i in vals.to.run) {
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
#       rownames(feature) <- feature$probe
    }
  }

  message("got feature of ", dataset.name)

#   colnames(feature) <- c("probe", "gene", "EntrezGene.ID")
  feature[feature==""|feature==" "] <- NA
  feature <- na.omit(feature)
  #   rownames(feature) <- feature$probe
  ## x is entrez ID's with exprs(eset)
  expr <- expr[is.element(row.names(expr),rownames(feature)),]
  feature <- feature[is.element(rownames(feature),rownames(expr)),]
  x <- cbind(feature[,mapping.group],expr)
  
  
  ## probe gene mapping options : MaxMean, maxRowVariance, Average
  Rprof("profile_mm_wgcna.txt")
  collapsed.wgcna.mrv <- collapseRows(x[,-1], rowID=rownames(x), rowGroup=x[,1], method=mapping.method)
  Rprof(NULL)
  
  best_probe <- collapsed.wgcna.mrv[[3]]
  feature <- cbind(feature, best_probe)
  message("collapsed rows of ", dataset.name)
  ## set best_probe for ambiguous entrez to false
  feature[grep("///", feature[, mapping.group]), "best_probe"] <- FALSE
  feature <- feature[match(rownames(expr),rownames(feature)),]
  ### create the expressionsets
  feature.df <- AnnotatedDataFrame(feature)
  pData <- read.table(paste("./curation/breast/curated/",dataset.names[i], "_curated.txt", sep=""), header=TRUE, stringsAsFactors=FALSE)
  if(dataset.name == "TCGA"){
    pData <- pData[is.element(gsub("-", "\\.",pData[,"sample_name"]), colnames(expr)), ]
  }
  duplicates <- rep(NA, each=nrow(pData))
  pData <- cbind(pData, duplicates)
  pData <- AnnotatedDataFrame(pData)
  rownames(attributes(pData)$data) <- colnames(expr)
  
  #!!! expression set MUST BE MATRIX
  # !!! Colnames of exprs must match rownames of pData (no transpose of exprs!)
  eset <- ExpressionSet(assayData=(as.matrix(expr)), featureData=feature.df, phenoData=pData)
  
  if(is.element(dataset.name, c("EXPO", "TRANSBIG") )){
    sampleNames(eset) <- paste(dataset.name, "_", sampleNames(eset), sep="")
  }

  experimentData(eset) <- MIAME(pubMedIds=as.character(eData[eData$Dataset == dataset.name, "PMID"]), 
                                contact=as.character(eData[eData$Dataset == dataset.name, "Paper"]),
                                url=as.character(eData[eData$Dataset == dataset.name, "Link.to.original.dataset"]),
                                abstract=as.character(eData[eData$Dataset == dataset.name, "Objectives"]),
                                other=list(summary = as.character(eData[eData$Dataset == dataset.name, "Short.Summary.of.Results"]),
        
                                           mapping.method=mapping.method, 
                                           mapping.group=mapping.group,
                                           preprocessing= "As published by original author."))
  assign(as.character(dataset.name), eset)

  save(list=as.character(dataset.name), file =paste("esets/mapped_esets2/",dataset.names[i], "_eset.rda", sep=""))
  
}

#If run in parallel, this code section needs to be run separately after the above finishes
      if(Sys.getenv("SGE_TASK_ID") == "") {
        ## identify duplicates and annotate pData
        source("./R/benDuplicateFinder.R")
        message("Duplicates Identified. Annotating pData.")
        for(i in 1:length(remove)){
          replicates <- remove[[i]]
          for(n in 1:length(replicates)){
        #     dataset.name <- as.character(gsub("\\..+", "", replicates)[n])
            dataset.name <- as.character(unlist(strsplit(x=replicates[n], split="\\."))[1])
            eset <- get(dataset.name)
            pData(eset)[paste(as.character(unlist(strsplit(x=replicates[n], split="\\."))[-1]), collapse="."), "duplicates"] <- paste(replicates[-n], collapse="///")
            assign(as.character(dataset.name), eset)
          }
        }
        for(e in 1:length(dataset.names)){
          save(list=as.character(dataset.names[e]), file =paste("esets/mapped_esets2/",dataset.names[e], "_eset.rda", sep=""))
        }
      }

