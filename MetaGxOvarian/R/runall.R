## runall
## excludes TCGA.mirna.8x15kv2_eset

## 1) first download all esets from FULLVcuratedOvarianData in FULLVdata
if(!file.exists("./FULLVdata")){
	stop("Please download original expression sets from FULLVcuratedOvarianData and store in FULLVdata.")
}

source("./R/getDatafromEsets.R")
source("./R/FeatureAnnotation_curation.R")
source("./R/createEsets.R") #Note: calls duplicate finder internally
createEsets()