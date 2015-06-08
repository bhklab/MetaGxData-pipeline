# merge ovarian datasets
library(MetaGxOvarian)
source(system.file("extdata", "patientselection.config",  package="MetaGxOvarian"))
source(system.file("extdata", "createEsetList.R", package="MetaGxOvarian"))

source("./R/datasetMerging.R")

merged.ovarian <- datasetMerging(esets)

save(merged.ovarian, file="./esets/summary/ovarian/mergedOvarian.rda")