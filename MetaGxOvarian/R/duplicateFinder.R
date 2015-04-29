`duplicateFinder`<- function (eset, topvar.genes=1000, dupl.cor=0.99, method=c("pearson", "spearman", "kendall"), cor.matrix=FALSE, nthread=1) {
  
  method <- match.arg(method)
  if (topvar.genes < 3) { topvar.genes <- length(featureNames(eset)) }
  ## select the most variant genes
  ## at least in 80% of the datasets
  iix <- apply(exprs(eset), 1, function (x, y) {
    return ((sum(is.na(x)) / length(x)) < ( 1 - y))
  }, y=0.8)
  varg <- Biobase::featureNames(eset)[iix][order(apply(exprs(eset)[iix, , drop=FALSE], 1, var, na.rm=TRUE), decreasing=TRUE)[1:topvar.genes]]
  
  
  ## using mRMRe
  nn <- mRMRe::get.thread.count()
  mRMRe::set.thread.count(nthread)
  expr <- mRMRe::mRMR.data(data=data.frame(Biobase::exprs(eset)[varg, , drop=FALSE]))
  cor.samples <- mRMRe::mim(object=expr, continuous_estimator=method, method="cor")
  mRMRe::set.thread.count(nn)
  
  if (cor.matrix) { return (cor.samples) }
  diag(cor.samples) <- NA
  ## create list of duplicates for each sample
  duplix <- apply(cor.samples, 1, function (x, y) {
    res <- names(x)[!is.na(x) & x > y]
    return (res)
  }, y=dupl.cor)
  duplix <- duplix[sapply(duplix, length) > 0]
  return (duplix)
}
