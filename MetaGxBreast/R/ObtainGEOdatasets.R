# Deena M.A. Gendoo
# September 25, 2016

# Download Breast Cancer Datasets from GEO
# Use GEOquery to extract exprs, pData, fData
# NB: pData has to be externally formatted (heavy formatting) before re-reading into R
  #- Add tissue, dataset, series columns
  #- Rename some columns
# Accordingly, there are txt files with the externally formated pData
# These files are then modified before creating the RData object!

library(GEOquery)
library(data.table)

#/mnt/work1/users/bhklab/Data/MetaGxData/MetaGxBreast/GEOdata

#Read in list of GSE datasets
if(file.exists("datafilesGEO.txt"))
  readdata.files2 <- read.table("datafilesGEO.txt", as.is=TRUE)[, 1]
#list of GSEs to be processed
sep=".RData"
GSEs <- unlist(lapply(readdata.files2, function(x) strsplit(x, sep)[[1]][1]))

#############
# GSE25066 # 
#############

  #Get Expression Set for the Series Matrix
  gse <- getGEO(GEO = "GSE25066",GSEMatrix = TRUE)
  gse<-gse[[1]]
  
  #Get featureData
  annot<-fData(gse)
  #rename some columns and add some new ones
  #colnames(annot)[c("ID","Gene Title","Gene Symbol","ENTREZ_GENE_ID")]
  colnames(annot)[c(1,10,11,12)]<-c("probe","Gene.Title","Gene.symbol","EntrezGene.ID")
  annot$Gene.ID<-annot$EntrezGene.ID
  
  # Get Gene Expression data
  data<-t(exprs(gse))
  
  #Get phenoData.
  # NB: phenoData [pData(gse)] has been cleaned externally (some heavy formatting required!) before further modification in R
  demo<-read.delim("GSE25066_pData.txt")
  # Renamed drfs_1_event_0_censorted to e.dmfs
  # Renamed drfs_even_time_years to t.dmfs
  
  rownames(demo)<-demo$X
  demo<-demo[,-1]
  sep=": "
  #age in years
  demo$age_at_initial_pathologic_diagnosis<-as.numeric(unlist(lapply(as.character(demo$age_at_initial_pathologic_diagnosis), function(x) strsplit(x, sep)[[1]][2])))
  demo$er<-unlist(lapply(as.character(demo$er), function(y) strsplit(y, sep)[[1]][2]))
  demo$pgr<-unlist(lapply(as.character(demo$pgr), function(y) strsplit(y, sep)[[1]][2]))
  demo$her2<-unlist(lapply(as.character(demo$her2), function(y) strsplit(y, sep)[[1]][2]))
  demo$T<-unlist(lapply(as.character(demo$T), function(y) strsplit(y, sep)[[1]][2]))
  demo$N<-unlist(lapply(as.character(demo$N), function(y) strsplit(y, sep)[[1]][2]))
  demo$grade<-unlist(lapply(as.character(demo$grade), function(y) strsplit(y, sep)[[1]][2]))
  demo$grade<-unlist(lapply(as.character(demo$grade), function(y) strsplit(y, "=")[[1]][1]))
  demo$e.dmfs<-as.numeric(unlist(lapply(as.character(demo$e.dmfs), function(x) strsplit(x, sep)[[1]][2])))
  demo$esr1_status<-unlist(lapply(as.character(demo$esr1_status), function(y) strsplit(y, sep)[[1]][2]))
  demo$erbb2_status<-unlist(lapply(as.character(demo$erbb2_status), function(y) strsplit(y, sep)[[1]][2]))
  demo$chemosensitivity_prediction<-unlist(lapply(as.character(demo$chemosensitivity_prediction), function(y) strsplit(y, sep)[[1]][2]))
  demo$GGI_prediction<-unlist(lapply(as.character(demo$GGI_prediction), function(y) strsplit(y, sep)[[1]][2]))
  demo$PAM50_prediction<-unlist(lapply(as.character(demo$PAM50_prediction), function(y) strsplit(y, sep)[[1]][2]))
  demo$dlda30_prediction<-unlist(lapply(as.character(demo$dlda30_prediction), function(y) strsplit(y, sep)[[1]][2]))
  demo$RCB_prediction<-unlist(lapply(as.character(demo$RCB_prediction), function(y) strsplit(y, sep)[[1]][2]))
  
  demo$t.dmfs<-as.numeric(unlist(lapply(as.character(demo$t.dmfs), function(x) strsplit(x, sep)[[1]][2])))
  days.per.year <- 365.242 
  sep=": "
  demo$t.dmfs<-demo$t.dmfs*days.per.year
  
  # Assume drfs (distant relapse free survival) equates with dmfs (distant metastasis free survival)
  demo$er[demo$er == "P"]<-1
  demo$er[demo$er == "N"]<-0
  
  demo$pgr[demo$pgr == "P"]<-1
  demo$pgr[demo$pgr == "N"]<-0
  
  demo$her2[demo$her2 == "P"]<-1
  demo$her2[demo$her2 == "N"]<-0
  
  demo$N<-unlist(lapply(as.character(demo$N), function(y) strsplit(y, "N")[[1]][2]))
  
  demo$treatment<-99 #NA for treatment
  
  #Save all into an R file
  save(data,annot,demo,file="GSE25066.RData")

#############
# GSE32646 # 
#############
  
  #Get Expression Set for the Series Matrix
  gse <- getGEO(GEO = "GSE32646",GSEMatrix = TRUE)
  gse<-gse[[1]]
  
  #Get featureData
  annot<-fData(gse)
  #rename some columns and add some new ones
  #colnames(annot)[c("ID","Gene Title","Gene Symbol","ENTREZ_GENE_ID")]
  colnames(annot)[c(1,10,11,12)]<-c("probe","Gene.Title","Gene.symbol","EntrezGene.ID")
  annot$Gene.ID<-annot$EntrezGene.ID
  
  # Get Gene Expression data
  data<-t(exprs(gse))
  
  #Get phenoData.
  # NB: phenoData [demo<-pData(gse)] has been cleaned externally (some heavy formatting required!) before further modification in R
  demo<-read.delim("GSE32646_pData.txt")
  rownames(demo)<-demo$X
  demo<-demo[,-1]
  sep=": "
  demo$age_at_initial_pathologic_diagnosis<-as.numeric(unlist(lapply(as.character(demo$age_at_initial_pathologic_diagnosis), function(x) strsplit(x, sep)[[1]][2])))
  demo$er<-unlist(lapply(as.character(demo$er), function(y) strsplit(y, sep)[[1]][2]))
  demo$pgr<-unlist(lapply(as.character(demo$pgr), function(y) strsplit(y, sep)[[1]][2]))
  demo$her2<-unlist(lapply(as.character(demo$her2), function(y) strsplit(y, sep)[[1]][2]))
  demo$T<-unlist(lapply(as.character(demo$T), function(y) strsplit(y, sep)[[1]][2]))
  demo$N<-unlist(lapply(as.character(demo$N), function(y) strsplit(y, sep)[[1]][2]))
  demo$grade<-unlist(lapply(as.character(demo$grade), function(y) strsplit(y, sep)[[1]][2]))

  demo$N[demo$N == "positive"]<-1
  demo$N[demo$N == "negative"]<-0
  
  demo$treatment<-99 #NA for treatment
  
  #Save all into an R file
  save(data,annot,demo,file="GSE32646.RData")
  
  
#############
# GSE58644 # 
#############
  
  #Get Expression Set for the Series Matrix
  gse <- getGEO(GEO = "GSE58644",GSEMatrix = TRUE)
  gse<-gse[[1]]
  
  #Get featureData
  annot<-fData(gse)
  annot<-annot[,-c(2:10)]
  
  #rename some columns and add some new ones
  #colnames(annot)[c("ID","Gene Title","Gene Symbol","ENTREZ_GENE_ID")]
  colnames(annot)[c(1)]<-c("probe")
  
  # Some more fixes
  gpl6244<-getGEO('GPL6244',destdir = ".")
  colnames(Table(gpl6244))
  gpltable<-data.frame(Table(gpl6244))
  gpltable$genedescription<-unlist(lapply(as.character(gpltable$gene_assignment), function(y) strsplit(y, "//")[[1]][3]))
  gpltable$entrez<-as.numeric(unlist(lapply(as.character(gpltable$gene_assignment), function(y) strsplit(y, "//")[[1]][5])))
  gpltable$gene_assignment<-unlist(lapply(as.character(gpltable$gene_assignment), function(y) strsplit(y, "//")[[1]][2]))
  
  annot$Gene.Symbol<-gpltable$gene_assignment[match(annot$probe,gpltable$ID)]
  annot$Gene.Title<-gpltable$genedescription[match(annot$probe,gpltable$ID)]
  annot$EntrezGene.ID<-gpltable$entrez[match(annot$probe,gpltable$ID)]
  annot$Gene.ID<-annot$EntrezGene.ID
  
  # Get Gene Expression data
  data<-t(exprs(gse))
  
  #Get phenoData.
  # NB: phenoData [pData(gse)] has been cleaned externally (some heavy formatting required!) before further modification in R
  demo<-read.delim("GSE58644_pData.txt")
  rownames(demo)<-demo$X
  demo<-demo[,-1]
  sep=": "
  demo$age_at_initial_pathologic_diagnosis<-as.numeric(unlist(lapply(as.character(demo$age_at_initial_pathologic_diagnosis), function(x) strsplit(x, sep)[[1]][2])))
  demo$er<-unlist(lapply(as.character(demo$er), function(y) strsplit(y, sep)[[1]][2]))
  demo$her2<-unlist(lapply(as.character(demo$her2), function(y) strsplit(y, sep)[[1]][2]))
  demo$T<-unlist(lapply(as.character(demo$T), function(y) strsplit(y, sep)[[1]][2]))
  demo$N<-unlist(lapply(as.character(demo$N), function(y) strsplit(y, sep)[[1]][2]))
  demo$grade<-unlist(lapply(as.character(demo$grade), function(y) strsplit(y, sep)[[1]][2]))
  demo$e.dmfs<-unlist(lapply(as.character(demo$vital_status), function(y) strsplit(y, sep)[[1]][2]))
  demo$tumor_size<-unlist(lapply(as.character(demo$tumor_size), function(y) strsplit(y, sep)[[1]][2]))
  
  demo$chemo<-unlist(lapply(as.character(demo$chemo), function(y) strsplit(y, sep)[[1]][2]))
  demo$tamoxifen<-unlist(lapply(as.character(demo$tamoxifen), function(y) strsplit(y, sep)[[1]][2]))
  demo$herceptin<-unlist(lapply(as.character(demo$herceptin), function(y) strsplit(y, sep)[[1]][2]))
  
  demo$treatment<-99
  demo$treatment[(demo$chemo==NA & demo$tamoxifen==NA & demo$herceptin==NA)]<-99 #NA for treatment
  demo$treatment[(demo$chemo==0 & demo$tamoxifen==0 & demo$herceptin==0)]<-0 #no treatment
  demo$treatment[(demo$chemo==1 & demo$tamoxifen==0)]<-1 #chemo only
  demo$treatment[demo$tamoxifen==1]<-2 #hormonal only
  demo$treatment[(demo$chemo==1 & demo$tamoxifen==1)]<-6 #chemo + hormono
  demo$treatment[(demo$chemo==0 & demo$tamoxifen==1 & demo$herceptin==0)]<-2 #hormonal only
  
  #Assume TIME variable represents days_to_death
  days.per.year <- 365.242 
  sep=": "
  demo$t.dmfs<-demo$t.dmfs*days.per.year
  
  #HOW TO DEAL WITH ANNOT?!
  
  #Save all into an R file
  save(data,annot,demo,file="GSE58644.RData")

  
#############
# GSE48091 # 
#############
  
  #Get Expression Set for the Series Matrix
  gse <- getGEO(GEO = "GSE48091",GSEMatrix = TRUE)
  gse<-gse[[1]]
  
  #Get featureData
  annot<-fData(gse)
  #rename some columns and add some new ones
  #colnames(annot)[c("ID","Gene Title","Gene Symbol","ENTREZ_GENE_ID")]
  colnames(annot)[c(1,2,3)]<-c("probe","EntrezGene.ID","Gene.symbol")
  annot$Gene.ID<-annot$EntrezGene.ID
  
  # Get Gene Expression data
  data<-t(exprs(gse))
  
  #Get phenoData.
  # NB: phenoData [pData(gse)] has been cleaned externally (some heavy formatting required!) before further modification in R
  demo<-read.delim("GSE48091_pData.txt")
  rownames(demo)<-demo$X
  demo<-demo[,-1]

  #Save all into an R file
  save(data,annot,demo,file="GSE48091.RData")
  
  
