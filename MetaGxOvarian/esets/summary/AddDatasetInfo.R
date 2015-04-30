#### Add ovarian dataset info ###

datasets <- read.csv("./datasets.csv")
load("./esets/summary/ovarian/SampleNumberSummaries_ovarian.rda")
patients <- rep(NA, each=nrow(datasets))
for(i in 1:nrow(datasets)){
	if(is.element(datasets$Dataset[i], rownames(SampleNumberSummaryAll))){
		patients[i] <- SampleNumberSummaryAll[as.character(datasets$Dataset[i]), "NumberOfSamples"]
	}
}
datasets <- cbind(datasets, patients)


source("./patientselection.config")
min.number.of.genes <- 0
rm(remove.duplicates)
rm(probe.gene.mapping)
source("./R/createEsetList.R")

probe <- NULL
genes <- NULL
for (eset in esets){
	probe <- c(probe, length(featureNames(eset)))
	genes <- c(genes, sum(fData(eset)$best_probe))
}

tmp <- data.frame(probe=probe, genes=genes, row.names=names(esets))
genes <- tmp[as.character(datasets$Dataset), "genes"]
probes <- tmp[as.character(datasets$Dataset), "probe"]

datasets <- cbind(datasets, genes, probes)

write.csv(datasets, file="./esets/summary/ovarian/datasets.csv")