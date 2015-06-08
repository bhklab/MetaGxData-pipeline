#### Add ovarian dataset info ###

#datasets <- read.csv("./datasets.csv")
#load("./esets/summary/ovarian/SampleNumberSummaries_ovarian.rda")
#patients <- rep(NA, each=nrow(datasets))
#for(i in 1:nrow(datasets)){
#	if(is.element(datasets$Dataset[i], rownames(SampleNumberSummaryAll))){
#		patients[i] <- SampleNumberSummaryAll[as.character(datasets$Dataset[i]), "NumberOfSamples"]
#	}
#}
#datasets$patients <- patients
datasets <- read.csv("./datasets.csv")
library(MetaGxOvarian)
source(system.file("extdata", "patientselection.config", package="MetaGxOvarian"))
min.number.of.genes <- 0
rm(remove.duplicates)
rm(probe.gene.mapping)
rm(remove.subsets)
source(system.file("extdata", "createEsetList.R", package="MetaGxOvarian"))

probe <- rep(NA, each=nrow(datasets))
genes <- rep(NA, each=nrow(datasets))
patients <- rep(NA, each=nrow(datasets))
datasets$probes <- probe
datasets$genes <- genes
datasets$patients <- patients

for (i in 1:length(esets)){
    if(is.element(names(esets)[i], datasets$Dataset))
	datasets$probes[which(datasets$Dataset==names(esets)[i])] <- length(featureNames(esets[[i]]))
	datasets$genes[which(datasets$Dataset==names(esets)[i])] <- sum(fData(esets[[i]])$best_probe)
	datasets$patients[which(datasets$Dataset==names(esets)[[i]])] <- length(sampleNames(esets[[i]]))
}


write.csv(datasets, file="./esets/summary/ovarian/datasets.csv")
write.csv(datasets, file="./datasets.csv")