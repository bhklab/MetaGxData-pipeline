## Some code for a random forest gene-pairs classifier

library(MetaGxBreast)

source(system.file("extdata", "patientselection.config", package="MetaGxBreast"))
source(system.file("extdata", "createEsetList.R", package="MetaGxBreast"))

esets <- lapply(esets, function(x) {
  factor.indices <- sapply(pData(x), is.factor)
  pData(x)[factor.indices] <- lapply(pData(x)[factor.indices], as.character)
  return(x)
  })
esets.subtype.scmod2 <- lapply(esets, function(x) {
  x <- subtypeClassification(x, model = "scmod2")
  return(x)
  })
esets.subtype.pam50 <- lapply(esets, function(x) {
  x <- subtypeClassification(x, model = "pam50")
  return(x)
})

esets.subtype.genepairs <- lapply(esets, function(x) {
  x <- getGenePairBreastCancerSubtypes(x)
  return(x)
})

eset.names <- names(esets)
esets <- lapply(names(esets), function(x) {
  eset.toreturn <- esets[[x]]
  eset.toreturn$data.source <- x
  return(eset.toreturn)
  })
names(esets) <- eset.names

getGenePairBreastCancerSubtypes <- function(eset) {
  
  ## TODO: rescale both metabric.eset and eset per gene (e.g. exprs(eset) <- t(scale(t(exprs(eset)))))
  
  metabric.eset <- esets.subtype.scmod2$METABRIC
  
  # Remove genes with NA values
  eset <- eset[complete.cases(exprs(eset)),]
  metabric.eset <- metabric.eset[complete.cases(exprs(metabric.eset)),]
  
  # Get top 20 entrez gene ids. Node that scmod2.robust$mod matrices are ordered in decreasing absolute value of coefficient.
  top.20.entrez.gene.ids <- as.numeric(sapply(genefu::scmod2.robust$mod, function(mod) mod[1:20,2]))
  entrez.gene.ids <- as.numeric(top.20.entrez.gene.ids)
  
  train.labels <- metabric.eset$subtype
  
  intersecting.entrez.ids <- intersect(as.character(fData(metabric.eset)$EntrezGene.ID), as.character(fData(eset)$EntrezGene.ID))
  intersecting.entrez.ids <- intersect(intersecting.entrez.ids, entrez.gene.ids)
  
  print("Training Random Forest classifier...")
  
  train.expression.matrix <- t(exprs(metabric.eset)[match(intersecting.entrez.ids, fData(metabric.eset)$EntrezGene.ID),])
  
  train.pairwise.matrix <-
    apply(combn(1:length(intersecting.entrez.ids),2), 2, function(pair) train.expression.matrix[,pair[1]] > train.expression.matrix[,pair[2]])
  train.pairwise.vals <- as.data.frame(train.pairwise.matrix)
  train.pairwise.vals[] <- lapply(train.pairwise.vals, function(x) factor(x, levels=c("FALSE", "TRUE")))
  
  rf.model <- randomForest(x=train.pairwise.vals, y=train.labels)
  
  test.expression.matrix <- t(exprs(eset)[match(intersecting.entrez.ids, fData(eset)$EntrezGene.ID),])
  
  test.pairwise.matrix <-
    apply(combn(1:length(intersecting.entrez.ids),2), 2, function(pair) test.expression.matrix[,pair[1]] > test.expression.matrix[,pair[2]])
  test.pairwise.vals <- as.data.frame(test.pairwise.matrix)
  test.pairwise.vals[] <- lapply(test.pairwise.vals, function(x) factor(x, levels=c("FALSE", "TRUE")))
  
  
  my.predictions <- predict(rf.model, newdata = test.pairwise.vals)
  #my.predictions.probs <- predict(nb.model, newdata = test.pairwise.matrix, type = 'prob')
  
  eset$subtypes.genepairs.rf.2 <- my.predictions
  return(eset)
}