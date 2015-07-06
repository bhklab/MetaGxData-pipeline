rm(list=ls(all=TRUE))

## create directory
saveres <- "saveres"
if(!file.exists(saveres)) { dir.create(saveres, showWarnings=FALSE) }

nbcore <- 16
availcore <- parallel::detectCores()
if (is.null(nbcore) || nbcore > availcore) { ncore <- availcore }
options("mc.cores"=nbcore)

## all human gene symbols
library(org.Hs.eg.db)
gs <- toTable(org.Hs.egSYMBOL)
## get all gene symbols
gs2 <- as.character(unique(gs[ , "symbol"]))

## retrieve the breast cancer data
library(cgdsr)
# Create CGDS object
mycgds <- CGDS("http://www.cbioportal.org/public-portal/")
## test the connection
test(mycgds)

# Get list of cancer studies at server
## get the breast cancer data
mycancerstudy <- getCancerStudies(mycgds)[15, 1] #Study is brca_tcga (PROVISIONAL!)
## get all patient ids
mycaselist <- getCaseLists(mycgds,mycancerstudy)[9, 1] #Study is brca_tcga_rna_seq_v2_mrna, 1098 samples

## Get the list of gene expression data
mygeneticprofile <- getGeneticProfiles(mycgds,mycancerstudy) 
print(mygeneticprofile)
## 4: Log2 copy-number values for each gene (from Affymetrix SNP6).
## 6: Methylation (HM450) beta-values for genes in 616 cases. For genes with multiple methylation probes, the probe least correlated with expression
## 8: Expression levels for 20532 genes in 914 brca cases (RNA Seq V2 RSEM)


## let's start with gene expression
mygeneticprofile <- getGeneticProfiles(mycgds,mycancerstudy)[3, 1] #brca_tcga_pub_mrna is SELECTED
## parallel version
splitix <- parallel::splitIndices(nx=length(gs2), ncl=ceiling(length(gs2) / 100))
mcres <- parallel::mclapply(splitix, function(x, mycgds, gs2, mygeneticprofile, mycaselist) {
	dd <- getProfileData(x=mycgds, genes=gs2[x], geneticProfiles=mygeneticprofile, caseList=mycaselist)
  cix <- apply(dd, 2, function(x) { return(all(is.na(x))) })
	dd <- dd[ , !cix, drop=FALSE]
	return(dd)
}, mycgds=mycgds, gs2=gs2, mygeneticprofile=mygeneticprofile, mycaselist=mycaselist, mc.cores=nbcore)

# Check if any NAs, get rid of those
#   test <- NULL
#   for(i in 1:length(mcres)){test <- c(test, ncol(mcres[[i]]))}
#   mcres2 <- mcres[-which(test==0)]
#   lapply(mcres2,FUN=dim) #check if any expr sets are empty!
#   #data.ge<-lapply(mcres, function(n){ if (ncol(n)>0 && nrow(n)==526) {do.call(cbind,mcres)}})

data.ge <- do.call(cbind, mcres) #Get 1098 samples x 20498 genes
rix <- apply(data.ge, 1, function(x) { return(all(is.na(x))) })
data.ge <- data.ge[!rix, , drop=FALSE]
data.ge <- log2(data.ge + 1)
data.ge <- data.matrix(data.ge)
#     
#     ## retrieval of copy number variations
#     mygeneticprofile <- getGeneticProfiles(mycgds,mycancerstudy)[4, 1]
#     ## parallel version
#     splitix <- parallel::splitIndices(nx=length(gs2), ncl=ceiling(length(gs2) / 100))
#     mcres <- parallel::mclapply(splitix, function(x, mycgds, gs2, mygeneticprofile, mycaselist) {
#     	dd <- getProfileData(x=mycgds, genes=gs2[x], geneticProfiles=mygeneticprofile, caseList=mycaselist)
#       cix <- apply(dd, 2, function(x) { return(all(is.na(x))) })
#     	dd <- dd[ , !cix, drop=FALSE]
#     	return(dd)
#     }, mycgds=mycgds, gs2=gs2, mygeneticprofile=mygeneticprofile, mycaselist=mycaselist, mc.cores=nbcore)
#     data.cnv <- do.call(cbind, mcres)
#     rix <- apply(data.cnv, 1, function(x) { return(all(is.na(x))) })
#     data.cnv <- data.cnv[!rix, , drop=FALSE]
#     
#     ## retrieval of methylation
#     mygeneticprofile <- getGeneticProfiles(mycgds,mycancerstudy)[6, 1]
#     ## parallel version
#     splitix <- parallel::splitIndices(nx=length(gs2), ncl=ceiling(length(gs2) / 100))
#     mcres <- parallel::mclapply(splitix, function(x, mycgds, gs2, mygeneticprofile, mycaselist) {
#     	dd <- getProfileData(x=mycgds, genes=gs2[x], geneticProfiles=mygeneticprofile, caseList=mycaselist)
#       cix <- apply(dd, 2, function(x) { return(all(is.na(x))) })
#     	dd <- dd[ , !cix, drop=FALSE]
#     	return(dd)
#     }, mycgds=mycgds, gs2=gs2, mygeneticprofile=mygeneticprofile, mycaselist=mycaselist, mc.cores=nbcore)
#     data.methyl <- do.call(cbind, mcres)
#     rix <- apply(data.methyl, 1, function(x) { return(all(is.na(x))) })
#     data.methyl <- data.methyl[!rix, , drop=FALSE]
#     
#     ## retrieval of mutations
#     mygeneticprofile <- getGeneticProfiles(mycgds,mycancerstudy)[7, 1]
#     ## parallel version
#     splitix <- parallel::splitIndices(nx=length(gs2), ncl=ceiling(length(gs2) / 100))
#     mcres <- parallel::mclapply(splitix, function(x, mycgds, gs2, mygeneticprofile, mycaselist) {
#     	dd <- getProfileData(x=mycgds, genes=gs2[x], geneticProfiles=mygeneticprofile, caseList=mycaselist)
#       cix <- apply(dd, 2, function(x) { return(all(is.na(x))) })
#     	dd <- dd[ , !cix, drop=FALSE]
#     	return(dd)
#     }, mycgds=mycgds, gs2=gs2, mygeneticprofile=mygeneticprofile, mycaselist=mycaselist, mc.cores=nbcore)
#     data.mut <- do.call(cbind, mcres)
#     rix <- apply(data.mut, 1, function(x) { return(all(is.na(x))) })
#     data.mut <- data.mut[!rix, , drop=FALSE]


## get the clinical data
data.clin <- getClinicalData(mycgds, mycaselist)

TCGA_exprs<-data.ge
TCGA_annot<-data.clin
save(TCGA_exprs,TCGA_annot,file="TCGAeset.RData")

# Collect TCGA Clinical Information from the Cancer Genome Atlas Portal:
#   https://wiki.nci.nih.gov/display/TCGA/Downloading+Clinical+Data+Using+the+Data+Matrix
#   https://tcga-data.nci.nih.gov/tcga/dataAccessFileProcessing.htm
#   https://tcga-data.nci.nih.gov/tcga/dataAccessMatrix.htm
#   http://cancergenome.nih.gov/publications/publicationguidelines

#Read in the Data from the Cancer Genome Atlas Data Portal (Clinical Info on Patients, most up-to-date)
#Select from the pData what you need to match the curatedBreastDatasets
ClinicalTCGA<-read.delim("nationwidechildrens.org_clinical_patient_brca.txt")
ClinicalTCGA_sub<-ClinicalTCGA[,c(1,2,12,14,15,16,17,21,44,50,56,93,109)]

    #TCGA_annot_sub<-TCGA_annot[,c(1,4,8,10,11,31,32,43)]
ClinicalTCGA_sub <- ClinicalTCGA_sub[-c(1,2),]
write.csv(ClinicalTCGA_sub, "../curation/breast/uncurated/TCGA.csv")
#Final Cleanups
test<-TCGA_exprs
rownames(test)<-gsub(rownames(test),pattern=".01",replacement="")
rownames(test)<-gsub(rownames(test),pattern="\\.",replacement="-")

TCGA_exprs<-test
rownames(ClinicalTCGA_sub)<-ClinicalTCGA_sub$bcr_patient_barcode
save.image("TCGAcurations.RData")

#     library(RTCGAToolbox)
#     stddata = getFirehoseRunningDates(last=3)
#     stddata
#     brcaData = getFirehoseData (dataset="BRCA", runDate="20141206", Clinic=TRUE, RNAseq_Gene=TRUE, mRNA_Array=TRUE)
