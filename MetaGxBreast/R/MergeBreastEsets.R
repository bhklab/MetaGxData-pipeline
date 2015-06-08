
# Merge breast datasets
library(MetaGxBreast)
source(system.file("extdata", "patientselection.config",  package="MetaGxBreast"))
source(system.file("extdata", "createEsetList.R", package="MetaGxBreast"))

source("./R/datasetMerging.R")

merged.breast <- datasetMerging(esets)

save(merged.breast, file="./esets/summary/breast/meargedBreast.rda")