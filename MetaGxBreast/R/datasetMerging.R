rescale <- function (x, na.rm = FALSE, q = 0) 
{
    if (q == 0) {
        ma <- max(x, na.rm = na.rm)
        mi <- min(x, na.rm = na.rm)
    }
    else {
        ma <- quantile(x, probs = 1 - (q/2), na.rm = na.rm)
        mi <- quantile(x, probs = q/2, na.rm = na.rm)
    }
    xx <- (x - mi)/(ma - mi)
    attributes(xx) <- list(names = names(x), q1 = mi, q2 = ma)
    return(xx)
}


`datasetMerging` <- 
  function (esets, method=c("union", "intersect"), nthread=1) {

    method <- match.arg(method)

    
    ## all unique Entrez gene ids
    ## gene ids
    ugid <- lapply(esets, function(x) { return(Biobase::fData(x)[, is.element(colnames(fData(x)), c("probeset", "EntrezGene.ID"))]) })
    ugid <- do.call(rbind, ugid)
    ugid <- ugid[!is.na(ugid[ , "EntrezGene.ID"]) & !duplicated(ugid[ , "EntrezGene.ID"]), , drop=FALSE]
    rownames(ugid) <- gsub(sprintf("(%s).", paste(names(esets), collapse="|")), "", rownames(ugid))
    switch (method,
            "union" = {
              feature.merged <- ugid
            },
            "intersect" = {
              feature.merged <- lapply(esets, function(x) { return(trimWhiteSpace(as.character(Biobase::fData(x)[ , "EntrezGene.ID"]))) })
              feature.merged <- table(unlist(feature.merged))
              feature.merged <- names(feature.merged)[feature.merged == length(esets)]
              feature.merged <- ugid[match(feature.merged, trimWhiteSpace(as.character(ugid[ , "EntrezGene.ID"]))), , drop=FALSE]
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

 # DG-22/09/2015: Want only normalization by patient, avoid scaling & normalization by gene
 # Mask genefu normalization for now!
 # ee <- apply(ee, 2, rescale)
#      splitix <- parallel::splitIndices(nx=ncol(ee), ncl=nthread)
#      mcres <- parallel::mclapply(splitix, function(x, data) {
#        res <- apply(data[ , x, drop=FALSE], 2, function (dx) {
#          return ((rescale(dx, q=0.05, na.rm=TRUE) - 0.5) * 2)
#        })
#        return (res)
#      }, data=ee, mc.cores=nthread)
#      ee <- do.call(cbind, mcres)


# impute.knn.args <- list(k = 10, rowmax = 0.5, colmax = 1, maxp = 1500, rng.seed=362436069)
# impute.knn.args$data <- ee
# ee <- do.call(impute::impute.knn, args=impute.knn.args)

 ## quantile normalization
 ee <- limma::normalizeBetweenArrays(object=ee, method="quantile")
 exprs(eset.merged) <- ee
 
return (eset.merged)
}

