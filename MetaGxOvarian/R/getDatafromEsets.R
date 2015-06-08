## getDatafromEsets.R
### assuming esets from FULLVdata are in folder FULLVdata
### reads names of esets from datalist file in FULLVdata

`datasetMerging` <- 
  function (esets, method=c("union", "intersect"), nthread=1) {

    method <- match.arg(method)

    
    ## all unique Entrez gene ids
    ## gene ids
    ugid <- lapply(esets, function(x) { return(Biobase::fData(x)) })
    ugid <- do.call(rbind, ugid)
    ugid <- ugid[!is.na(ugid[ , "probeset"]) & !duplicated(ugid[ , "probeset"]), , drop=FALSE]
    rownames(ugid) <- gsub(sprintf("(%s).", paste(names(esets), collapse="|")), "", rownames(ugid))
    switch (method,
            "union" = {
              feature.merged <- ugid
            },
            "intersect" = {
              feature.merged <- lapply(esets, function(x) { return(trimWhiteSpace(as.character(Biobase::fData(x)[ , "probeset"]))) })
              feature.merged <- table(unlist(feature.merged))
              feature.merged <- names(feature.merged)[feature.merged == length(esets)]
              feature.merged <- ugid[match(feature.merged, trimWhiteSpace(as.character(ugid[ , "probeset"]))), , drop=FALSE]
            },
{
  stop("Unknown method")
}
    )
## expression data
exprs.merged <- lapply(esets, function (x, y) {
  ee <- Biobase::exprs(x)[is.element(rownames(exprs(x)),rownames(feature.merged)),]
  print(dim(ee))
  eem <- matrix(NA, nrow=length(y), ncol=ncol(ee), dimnames=list(y, colnames(ee)))
  print(dim(eem))
  print(length(intersect(rownames(ee),rownames(eem))))
  eem[rownames(ee), colnames(ee)] <- ee
  return (eem)
}, y=rownames(feature.merged))
exprs.merged <- do.call(cbind, exprs.merged)
## clinical info
ucid <- lapply(esets, function(x) { return(colnames(Biobase::pData(x))) })
ucid <- table(unlist(ucid))
ucid <- names(ucid)[ucid == length(esets)]
clinicinfo.merged <- lapply(esets, function (x , y) {
  ee <- Biobase::pData(x)[ , y, drop=FALSE]
}, y=ucid)
clinicinfo.merged <- do.call(rbind, clinicinfo.merged)
rownames(clinicinfo.merged) <- colnames(exprs.merged)

#   ## create a merged expressionSet object
eset.merged <- ExpressionSet(assayData=exprs.merged, phenoData=AnnotatedDataFrame(data=clinicinfo.merged), featureData=AnnotatedDataFrame(data=feature.merged))
experimentData(eset.merged)@preprocessing <- list("normalization"="mixed", package="unspecified", version="0")
annotation(eset.merged) <- "mixed"
## standardization

 ## robust scaling followed by quantile normalization
 ee <- exprs(eset.merged)
 
 # ee <- apply(ee, 2, genefu::rescale)
 splitix <- parallel::splitIndices(nx=ncol(ee), ncl=nthread)
 mcres <- parallel::mclapply(splitix, function(x, data) {
   res <- apply(data[ , x, drop=FALSE], 2, function (dx) {
     return ((genefu::rescale(dx, q=0.05, na.rm=TRUE) - 0.5) * 2)
   })
   return (res)
 }, data=ee, mc.cores=nthread)
 ee <- do.call(cbind, mcres)

 ## quantile normalization
 ee <- limma::normalizeBetweenArrays(object=ee, method="quantile")
 exprs(eset.merged) <- ee
 
return (eset.merged)
}


nthread<-1
availcore <- parallel::detectCores()
if (nthread > availcore) { nthread <- availcore }
options("mc.cores"=nthread)

if(!file.exists("./FULLVdata/clinicalinfo")){
	dir.create("./FULLVdata/clinicalinfo", recursive=T)
}
if(!file.exists("./FULLVdata/expr")){
	dir.create("./FULLVdata/expr", recursive=T)
}

datasets <- read.table("./FULLVdata/datalist")[,1]
n <- which(datasets == "TCGA.mirna.8x15kv2_eset")
datasets <- datasets[-n]

## load all the datasets
for(i in 1:length(datasets)){
	dataset <- as.character(datasets[i])
	load(paste("./FULLVdata/", dataset, ".rda", sep=""))
}

for(i in 1:length(datasets)){

	dataset <- as.character(datasets[i])
  dn <- gsub(pattern="_eset", replacement="", x=dataset)
	if(dataset == "GSE19829.GPL570_eset" | dataset== "GSE19829.GPL8300_eset"){
		## these will be merged and dealt with separately
		next
	}
	load(paste("./FULLVdata/", dataset, ".rda", sep=""))
	## get expression data
	expr <- exprs(get(dataset))
    colnames(expr) <- paste(dn, "_", colnames(expr), sep="")
    
	if(dataset == "GSE32062.GPL6480_eset"){
		write.table(expr, file="./FULLVdata/expr/GSE32062_eset_exprs.txt", sep="\t")
		PData <- pData(get(dataset))
		write.table(PData, file="./FULLVdata/clinicalinfo/GSE32062_eset_pdata.txt", sep="\t")
	} else {
		write.table(expr, paste("./FULLVdata/expr/",dataset,"_exprs.txt", sep=""), sep="\t")
		## get clinical annotations (already curated)
		PData <- pData(get(dataset))
		write.table(PData, paste("./FULLVdata/clinicalinfo/", dataset,"_pdata.txt", sep=""), sep="\t")
	}


}

## merge these specific datasets since split between platforms
esetsToMerge <- c(get("GSE19829.GPL570_eset"),get("GSE19829.GPL8300_eset"))
names(esetsToMerge) <- c("GSE19829.GPL570_eset", "GSE19829.GPL8300_eset")
GSE19829_eset <- datasetMerging(esetsToMerge, method="union", nthread =nthread)
dataset <- "GSE19829_eset"
expr <- exprs(get(dataset))
colnames(expr) <- paste("GSE19829_", colnames(expr), sep="")
write.table(expr, paste("./FULLVdata/expr/GSE19829_eset_exprs.txt", sep=""), sep="\t")
PData <- pData(get(dataset))
write.table(PData, paste("./FULLVdata/clinicalinfo/GSE19829_eset_pdata.txt", sep=""), sep="\t")
experimentData(GSE19829_eset) <- experimentData(GSE19829.GPL8300_eset)
save(GSE19829_eset, file="./FULLVdata/GSE19829_eset.rda")


