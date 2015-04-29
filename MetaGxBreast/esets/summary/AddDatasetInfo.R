#### Add breast dataset info ###

datasets <- read.csv("./datasets.csv")
load("./esets/summary/breast/SampleNumberSummaries_breast.rda")
patients <- SampleNumberSummary[datasets$Dataset, "NumberOfSamples"]
datasets <- cbind(datasets, patients)

source("./patientselection.config")
min.number.of.genes <- 0
rm(remove.duplicates)
rm(probe.gene.mapping)
source("./R/createEsetList2.R")

probe <- NULL
genes <- NULL
for (eset in esets){
	probe <- c(probe, length(featureNames(eset)))
	genes <- c(genes, sum(fData(eset)$best_probe))
}

tmp <- data.frame(probe=probe, genes=genes, row.names=names(esets))
genes <- tmp[datasets$Dataset, "genes"]
probes <- tmp[datasets$Dataset, "probe"]

datasets <- cbind(datasets, genes, probes)

write.csv(datasets, file="./esets/summary/breast/datasets.csv")