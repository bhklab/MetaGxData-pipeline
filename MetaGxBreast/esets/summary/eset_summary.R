## BREAST ##

# ##---------------------------
# ## load and do probe gene mapping for esets
# ##---------------------------

# filenames <- paste("./esets/mapped_esets2/", list.files("./esets/mapped_esets2"),sep="")
# lapply(filenames, load, .GlobalEnv)

# datasets <- read.csv("datasets.csv")
# dataset.names <- datasets$Dataset
# sampleNames(EXPO) <- paste("EXPO_", sampleNames(EXPO), sep="")
# sampleNames(TRANSBIG) <- paste("TRANSBIG_", sampleNames(TRANSBIG), sep="")
# sampleNames(SUPERTAM_HGU133A) <- paste("SUPERTAM_HGU133A", sampleNames(SUPERTAM_HGU133A), sep="")
# sampleNames(SUPERTAM_HGU133PLUS2) <- paste("SUPERTAM_HGU133PLUS2", sampleNames(SUPERTAM_HGU133PLUS2), sep="")
# esets <- lapply(as.character(dataset.names), get)
# names(esets) <- dataset.names
# # save(esets, file="./esets/esets.rda)
# # load("./esets/esets.rda")
# esets.mapped <- list()
# for(eset in esets){
#   Biobase::exprs(eset) <- exprs(eset)[fData(eset)$best_probe,]
#   Biobase::fData(eset) <- fData(eset)[fData(eset)$best_probe,]
#   rownames(fData(eset)) <- rownames(exprs(eset)) <- fData(eset)$ENTREZID
#   esets.mapped <- c(esets.mapped, eset)
# }
# names(esets.mapped) <- names(esets)

# setwd("~/Documents/curatedBreastData/")

if(!file.exists("./esets/summary/breast")){
	dir.create("./esets/summary/breast", recursive=T)
}

library(MetaGxBreast)
source(system.file("extdata", "patientselection.config", package="MetaGxBreast"))
min.number.of.genes <- 0
rm(remove.duplicates)
source(system.file("extdata", "createEsetList.R", package="MetaGxBreast"))
esets.mapped <- esets


##---------------------------
## Number of samples
##---------------------------


## generate dataframe for summary of number of samples
numSamples <- NULL
for(i in 1:length(esets.mapped)){
	numSamples <- c(numSamples, length(sampleNames(esets.mapped[[i]])))
}

SampleNumberSummary <- data.frame(NumberOfSamples = numSamples, row.names = names(esets.mapped))
SampleNumberSummary <- rbind(SampleNumberSummary, sum(SampleNumberSummary[,"NumberOfSamples"]))

rownames(SampleNumberSummary)[nrow(SampleNumberSummary)] <- "Total"
SampleNumberSummaryPercent <- data.frame(Percent = SampleNumberSummary$NumberOfSamples/SampleNumberSummary$NumberOfSamples[nrow(SampleNumberSummary)])
rownames(SampleNumberSummaryPercent) <- rownames(SampleNumberSummary)
# save(SampleNumberSummary, file="./esets/summary/SampleNumberSummary.rda")
# save(SampleNumberSummaryPercent, file="./esets/summary/SampleNumberSummaryPercent.rda")
save(list = c("SampleNumberSummary", "SampleNumberSummaryPercent"), file = "./esets/summary/breast/SampleNumberSummaries_breast.rda")

##---------------------------
## Overall pData available
##---------------------------

pDataID <- c("er","pgr", "her2", "age_at_initial_pathologic_diagnosis", "grade", "dmfs_days", "dmfs_status", "days_to_tumor_recurrence", "recurrence_status", "days_to_death", "vital_status", "sample_type", "treatment")

pDataPercentSummaryTable <- NULL
pDataSummaryNumbersTable <- NULL
for(e in 1:length(esets.mapped)){
	eset <- esets.mapped[[e]]
	pDataPercentSummary <- NULL
	pDataSummaryNumbers <- NULL
	for(p in 1:length(pDataID)){
		pDataSummaryNumbers <- c(pDataSummaryNumbers, sum(!is.na(pData(eset)[,pDataID[p]])))
		pDataPercentSummary <- c(pDataPercentSummary, (sum(!is.na(pData(eset)[,pDataID[p]]))/nrow(pData(eset)))*100)

	}
	if(e == 1){
		pDataSummaryNumbersTable <- data.frame(test = pDataSummaryNumbers)
		pDataPercentSummaryTable <- data.frame(test = pDataPercentSummary)
	} else {
		pDataPercentSummaryTable <- cbind(pDataPercentSummaryTable,pDataPercentSummary)
		pDataSummaryNumbersTable <- cbind(pDataSummaryNumbersTable, pDataSummaryNumbers)
	}
}
rownames(pDataSummaryNumbersTable) <- pDataID
rownames(pDataPercentSummaryTable) <- pDataID
colnames(pDataSummaryNumbersTable) <- names(esets.mapped)
colnames(pDataPercentSummaryTable) <- names(esets.mapped)

pDataSummaryNumbersTable <- rbind(pDataSummaryNumbersTable, SampleNumberSummary["Total",])
rownames(pDataSummaryNumbersTable)[nrow(pDataSummaryNumbersTable)] <- "Total"

# save(pDataSummaryNumbersTable, file="./esets/summary/pDataSummaryNumbersTable.rda")
# save(pDataPercentSummaryTable, file="./esets/summary/pDataPercentSummaryTable.rda")
save(list=c("pDataSummaryNumbersTable", "pDataPercentSummaryTable"), file = "./esets/summary/breast/OverallpDataSummary_breast.rda")


##---------------------------
## Specific pData available
##---------------------------


## er summary
ER_numbers <- NULL
ER_percent <- NULL

for(i in 1:length(esets.mapped)){
	positive <- length(which(as.character(pData(esets.mapped[[i]])$er)=="positive"))
	negative <- length(which(as.character(pData(esets.mapped[[i]])$er)=="negative"))
	missing <- sum(is.na(pData(esets.mapped[[i]])$er))
	total <- length(sampleNames(esets.mapped[[i]]))
	tmp <- c(positive, negative, missing, total)
	tmpp <- tmp/tmp[4]*100
	if(i == 1){
		ER_numbers <- tmp
		ER_percent <- tmpp
	} else {
		ER_numbers <- rbind(ER_numbers, tmp)
		ER_percent <- rbind(ER_percent, tmpp)
	}
}
colnames(ER_numbers) <- colnames(ER_percent) <- c("positive", "negative", "missing", "total")
rownames(ER_numbers) <- rownames(ER_percent) <- names(esets.mapped)

ER_numbers <- t(ER_numbers)
ER_percent <- t(ER_percent)
total <- c(sum(ER_numbers["positive",])/sum(ER_numbers["total",]),
	sum(ER_numbers["negative",])/sum(ER_numbers["total",]), 
	sum(ER_numbers["missing",])/sum(ER_numbers["total",]), 
	sum(ER_numbers["total",])/sum(ER_numbers["total",]) )
total <- total*100
ER_percent <- cbind(ER_percent, total)
#save(ER_percent,file="./esets/summary/ER_percent.rda")
# save(ER_numbers,file="./esets/summary/ER_numbers.rda")



## pgr summary
PGR_numbers <- NULL
PGR_percent <- NULL

for(i in 1:length(esets.mapped)){
	positive <- length(which(as.character(pData(esets.mapped[[i]])$pgr)=="positive"))
	negative <- length(which(as.character(pData(esets.mapped[[i]])$pgr)=="negative"))
	missing <- sum(is.na(pData(esets.mapped[[i]])$pgr))
	total <- length(sampleNames(esets.mapped[[i]]))
	tmp <- c(positive, negative, missing, total)
	tmpp <- tmp/tmp[4]*100
	if(i == 1){
		PGR_numbers <- tmp
		PGR_percent <- tmpp

	} else {
		PGR_numbers <- rbind(PGR_numbers, tmp)
		PGR_percent <- rbind(PGR_percent, tmpp)
	}
}
colnames(PGR_numbers) <- colnames(PGR_percent) <- c("positive", "negative", "missing", "total")
rownames(PGR_numbers) <- rownames(PGR_percent) <- names(esets.mapped)
PGR_numbers <- t(PGR_numbers)
PGR_percent <- t(PGR_percent)
total <- c(sum(PGR_numbers["positive",])/sum(PGR_numbers["total",]),
	sum(PGR_numbers["negative",])/sum(PGR_numbers["total",]), 
	sum(PGR_numbers["missing",])/sum(PGR_numbers["total",]), 
	sum(PGR_numbers["total",])/sum(PGR_numbers["total",]) )
total <- total*100
PGR_percent <- cbind(PGR_percent, total)
#save(PGR_percent,file="./esets/summary/PGR_percent.rda")
# save(PGR_numbers,file="./esets/summary/PGR_numbers.rda")




## her2 summary

HER2_numbers <- NULL
HER2_percent <- NULL

for(i in 1:length(esets.mapped)){
	positive <- length(which(as.character(pData(esets.mapped[[i]])$her2)=="positive"))
	negative <- length(which(as.character(pData(esets.mapped[[i]])$her2)=="negative"))
	missing <- sum(is.na(pData(esets.mapped[[i]])$her2))
	total <- length(sampleNames(esets.mapped[[i]]))
	tmp <- c(positive, negative, missing, total)
	tmpp <- tmp/tmp[4]*100
	if(i == 1){
		HER2_numbers <- tmp
		HER2_percent <- tmpp

	} else {
		HER2_numbers <- rbind(HER2_numbers, tmp)
		HER2_percent <- rbind(HER2_percent, tmpp)
	}
}
colnames(HER2_numbers) <- colnames(HER2_percent) <- c("positive", "negative", "missing", "total")
rownames(HER2_numbers) <- rownames(HER2_percent) <- names(esets.mapped)
HER2_numbers <- t(HER2_numbers)
HER2_percent <- t(HER2_percent)
total <- c(sum(HER2_numbers["positive",])/sum(HER2_numbers["total",]),
	sum(HER2_numbers["negative",])/sum(HER2_numbers["total",]), 
	sum(HER2_numbers["missing",])/sum(HER2_numbers["total",]), 
	sum(HER2_numbers["total",])/sum(HER2_numbers["total",]) )
total <- total*100
HER2_percent <- cbind(HER2_percent, total)
#save(HER2_percent,file="./esets/summary/HER2_percent.rda")
# save(HER2_numbers,file="./esets/summary/HER2_numbers.rda")



## N summary

N_numbers <- NULL
N_percent <- NULL

for(i in 1:length(esets.mapped)){
	one <- length(which(pData(esets.mapped[[i]])$N==1))
	zero <- length(which(pData(esets.mapped[[i]])$N==0))
	missing <- sum(is.na(pData(esets.mapped[[i]])$N))
	total <- length(sampleNames(esets.mapped[[i]]))
	tmp <- c(zero, one, missing, total)
	tmpp <- tmp/tmp[4]*100
	if(i == 1){
		N_numbers <- tmp
		N_percent <- tmpp

	} else {
		N_numbers <- rbind(N_numbers, tmp)
		N_percent <- rbind(N_percent, tmpp)
	}
}
colnames(N_numbers) <- colnames(N_percent) <- c("zero", "one", "missing", "total")
rownames(N_numbers) <- rownames(N_percent) <- names(esets.mapped)
N_numbers <- t(N_numbers)
N_percent <- t(N_percent)
total <- c(sum(N_numbers["zero",])/sum(N_numbers["total",]),
	sum(N_numbers["one",])/sum(N_numbers["total",]), 
	sum(N_numbers["missing",])/sum(N_numbers["total",]), 
	sum(N_numbers["total",])/sum(N_numbers["total",]) )
total <- total*100
N_percent <- cbind(N_percent, total)
#save(N_percent,file="./esets/summary/N_percent.rda")
# save(N_numbers,file="./esets/summary/N_numbers.rda")

## grade summary

GRADE_numbers <- NULL
GRADE_percent <- NULL

for(i in 1:length(esets.mapped)){
	one <- length(which(pData(esets.mapped[[i]])$grade==1))
	two <- length(which(pData(esets.mapped[[i]])$grade==2))
	three <- length(which(pData(esets.mapped[[i]])$grade==3))
	zero <- length(which(pData(esets.mapped[[i]])$grade==0))
	missing <- sum(is.na(pData(esets.mapped[[i]])$grade))
	total <- length(sampleNames(esets.mapped[[i]]))
	tmp <- c(zero, one, two, three, missing, total)
	tmpp <- tmp/tmp[6]*100
	if(i == 1){
		GRADE_numbers <- tmp
		GRADE_percent <- tmpp

	} else {
		GRADE_numbers <- rbind(GRADE_numbers, tmp)
		GRADE_percent <- rbind(GRADE_percent, tmpp)
	}
}
colnames(GRADE_numbers) <- colnames(GRADE_percent) <- c("zero", "one", "two", "three", "missing", "total")
rownames(GRADE_numbers) <- rownames(GRADE_percent) <- names(esets.mapped)
GRADE_numbers <- t(GRADE_numbers)
GRADE_percent <- t(GRADE_percent)
total <- c(sum(GRADE_numbers["zero",])/sum(GRADE_numbers["total",]),
	sum(GRADE_numbers["one",])/sum(GRADE_numbers["total",]), 
	sum(GRADE_numbers["two",])/sum(GRADE_numbers["total",]), 
	sum(GRADE_numbers["three",])/sum(GRADE_numbers["total",]),  
	sum(GRADE_numbers["missing",])/sum(GRADE_numbers["total",]), 
	sum(GRADE_numbers["total",])/sum(GRADE_numbers["total",]) )
total <- total*100
GRADE_percent <- cbind(GRADE_percent, total)
#save(GRADE_percent,file="./esets/summary/GRADE_percent.rda")
# save(GRADE_numbers,file="./esets/summary/GRADE_numbers.rda")



## dmfs_status summary

dmfs_status_numbers <- NULL
dmfs_status_percent <- NULL

for(i in 1:length(esets.mapped)){
	deceased_or_recurrence <- length(which(as.character(pData(esets.mapped[[i]])$dmfs_status)=="deceased_or_recurrence"))
	living_norecurrence <- length(which(as.character(pData(esets.mapped[[i]])$dmfs_status)=="living_norecurrence"))
	missing <- sum(is.na(pData(esets.mapped[[i]])$dmfs_status))
	total <- length(sampleNames(esets.mapped[[i]]))
	tmp <- c(deceased_or_recurrence, living_norecurrence, missing, total)
	tmpp <- tmp/tmp[4]*100
	if(i == 1){
		dmfs_status_numbers <- tmp
		dmfs_status_percent <- tmpp

	} else {
		dmfs_status_numbers <- rbind(dmfs_status_numbers, tmp)
		dmfs_status_percent <- rbind(dmfs_status_percent, tmpp)
	}
}
colnames(dmfs_status_numbers) <- colnames(dmfs_status_percent) <- c("deceased_or_recurrence", "living_norecurrence", "missing", "total")
rownames(dmfs_status_numbers) <- rownames(dmfs_status_percent) <- names(esets.mapped)
dmfs_status_numbers <- t(dmfs_status_numbers)
dmfs_status_percent <- t(dmfs_status_percent)
total <- c(sum(dmfs_status_numbers["deceased_or_recurrence",])/sum(dmfs_status_numbers["total",]),
	sum(dmfs_status_numbers["living_norecurrence",])/sum(dmfs_status_numbers["total",]),
	 sum(dmfs_status_numbers["missing",])/sum(dmfs_status_numbers["total",]),
	  sum(dmfs_status_numbers["total",])/sum(dmfs_status_numbers["total",]) )
total <- total*100
dmfs_status_percent <- cbind(dmfs_status_percent, total)
#save(dmfs_status_percent,file="./esets/summary/dmfs_status_percent.rda")
# save(dmfs_status_numbers,file="./esets/summary/dmfs_status_numbers.rda")

## recurrence_status summary
recurrence_status_numbers <- NULL
recurrence_status_percent <- NULL

for(i in 1:length(esets.mapped)){
	deceased_or_recurrence <- length(which(as.character(pData(esets.mapped[[i]])$recurrence_status)=="deceased_or_recurrence"))
	living_norecurrence <- length(which(as.character(pData(esets.mapped[[i]])$recurrence_status)=="living_norecurrence"))
	missing <- sum(is.na(pData(esets.mapped[[i]])$recurrence_status))
	total <- length(sampleNames(esets.mapped[[i]]))
	tmp <- c(deceased_or_recurrence, living_norecurrence, missing, total)
	tmpp <- tmp/tmp[4]*100
	if(i == 1){
		recurrence_status_numbers <- tmp
		recurrence_status_percent <- tmpp

	} else {
		recurrence_status_numbers <- rbind(recurrence_status_numbers, tmp)
		recurrence_status_percent <- rbind(recurrence_status_percent, tmpp)
	}
}
colnames(recurrence_status_numbers) <- colnames(recurrence_status_percent) <- c("deceased_or_recurrence", "living_norecurrence", "missing", "total")
rownames(recurrence_status_numbers) <- rownames(recurrence_status_percent) <- names(esets.mapped)
recurrence_status_numbers <- t(recurrence_status_numbers)
recurrence_status_percent <- t(recurrence_status_percent)
total <- c(sum(recurrence_status_numbers["deceased_or_recurrence",])/sum(recurrence_status_numbers["total",]),
	sum(recurrence_status_numbers["living_norecurrence",])/sum(recurrence_status_numbers["total",]),
	 sum(recurrence_status_numbers["missing",])/sum(recurrence_status_numbers["total",]),
	  sum(recurrence_status_numbers["total",])/sum(recurrence_status_numbers["total",]) )
total <- total*100
recurrence_status_percent <- cbind(recurrence_status_percent, total)
#save(recurrence_status_percent,file="./esets/summary/recurrence_status_percent.rda")
# save(recurrence_status_numbers,file="./esets/summary/recurrence_status_numbers.rda")

## vital_status summary
vital_status_numbers <- NULL
vital_status_percent <- NULL

for(i in 1:length(esets.mapped)){
	deceased <- length(which(as.character(pData(esets.mapped[[i]])$vital_status)=="deceased"))
	living <- length(which(as.character(pData(esets.mapped[[i]])$vital_status)=="living"))
	missing <- sum(is.na(pData(esets.mapped[[i]])$vital_status))
	total <- length(sampleNames(esets.mapped[[i]]))
	tmp <- c(deceased, living, missing, total)
	tmpp <- tmp/tmp[4]*100
	if(i == 1){
		vital_status_numbers <- tmp
		vital_status_percent <- tmpp

	} else {
		vital_status_numbers <- rbind(vital_status_numbers, tmp)
		vital_status_percent <- rbind(vital_status_percent, tmpp)
	}
}
colnames(vital_status_numbers) <- colnames(vital_status_percent) <- c("deceased", "living", "missing", "total")
rownames(vital_status_numbers) <- rownames(vital_status_percent) <- names(esets.mapped)
vital_status_numbers <- t(vital_status_numbers)
vital_status_percent <- t(vital_status_percent)
total <- c(sum(vital_status_numbers["deceased",])/sum(vital_status_numbers["total",]),
	sum(vital_status_numbers["living",])/sum(vital_status_numbers["total",]),
	 sum(vital_status_numbers["missing",])/sum(vital_status_numbers["total",]),
	  sum(vital_status_numbers["total",])/sum(vital_status_numbers["total",]) )
total <- total*100
vital_status_percent <- cbind(vital_status_percent, total)
#save(vital_status_percent,file="./esets/summary/vital_status_percent.rda")
# save(vital_status_numbers,file="./esets/summary/vital_status_numbers.rda")

## treatment summary

treatment_numbers <- NULL
treatment_percent <- NULL

for(i in 1:length(esets.mapped)){
	untreated <- length(which(as.character(pData(esets.mapped[[i]])$treatment)=="untreated"))
	chemotherapy <- length(which(as.character(pData(esets.mapped[[i]])$treatment)=="chemotherapy"))
	hormonotherapy <- length(which(as.character(pData(esets.mapped[[i]])$treatment)=="hormonotherapy"))
	endocrine <- length(which(as.character(pData(esets.mapped[[i]])$treatment)=="endocrine"))
	chemo.plus.hormono <-length(which(as.character(pData(esets.mapped[[i]])$treatment)=="chemo.plus.hormono"))
	missing <- sum(is.na(pData(esets.mapped[[i]])$treatment))
	total <- length(sampleNames(esets.mapped[[i]]))
	tmp <- c(untreated, chemotherapy, hormonotherapy, endocrine, chemo.plus.hormono, missing, total)
	tmpp <- tmp/tmp[7]*100
	if(i == 1){
		treatment_numbers <- tmp
		treatment_percent <- tmpp

	} else {
		treatment_numbers <- rbind(treatment_numbers, tmp)
		treatment_percent <- rbind(treatment_percent, tmpp)
	}
}
colnames(treatment_numbers) <- colnames(treatment_percent) <- c("untreated", "chemotherapy", "hormonotherapy", "endocrine", "chemo.plus.hormono", "missing", "total")
rownames(treatment_numbers) <- rownames(treatment_percent) <- names(esets.mapped)
treatment_numbers <- t(treatment_numbers)
treatment_percent <- t(treatment_percent)
total <- c(sum(treatment_numbers["untreated",])/sum(treatment_numbers["total",]),
		sum(treatment_numbers["chemotherapy",])/sum(treatment_numbers["total",]),
		sum(treatment_numbers["hormonotherapy",])/sum(treatment_numbers["total",]),
	 	sum(treatment_numbers["endocrine",])/sum(treatment_numbers["total",]),
	 	sum(treatment_numbers["chemo.plus.hormono",])/sum(treatment_numbers["total",]),
	  	sum(treatment_numbers["missing",])/sum(treatment_numbers["total",]),
	  	sum(treatment_numbers["total",])/sum(treatment_numbers["total",]) ) 
total <- total*100
treatment_percent <- cbind(treatment_percent, total)
#save(treatment_percent,file="./esets/summary/treatment_percent.rda")
# save(treatment_numbers,file="./esets/summary/treatment_numbers.rda")


save(list=c("treatment_percent", "vital_status_percent", "recurrence_status_percent", "dmfs_status_percent", "GRADE_percent", "N_percent", "HER2_percent", "PGR_percent", "ER_percent"), file="./esets/summary/breast/pData_categary_percent_breast.rda")
save(list=c("treatment_numbers", "vital_status_numbers", "recurrence_status_numbers", "dmfs_status_numbers", "GRADE_numbers", "N_numbers", "HER2_numbers", "PGR_numbers", "ER_numbers"), file="./esets/summary/breast/pData_categary_numbers_breast.rda")



##-------
## Ranged Data Summary
##------

#rangedID <- c("tumor_size", "age_at_initial_pathologic_diagnosis", "dmfs_days", "days_to_tumor_recurrence", "percent_normal_cells", "percent_stromal_cells", "percent_tumor_cells")
rangedID <- c("tumor_size", "age_at_initial_pathologic_diagnosis", "dmfs_days", "days_to_tumor_recurrence")
ranged.data <- NULL
category <- NULL
eset.specific <- NULL

for(n in 1:length(rangedID)){
	total <- NULL
	for(i in 1:length(esets.mapped)){
		eset.specific <- pData(esets.mapped[[i]])[,rangedID[n]]
		category[[names(esets.mapped)[i]]] <- eset.specific
		total <- c(total, unlist(category[names(esets.mapped)[i]]))
	}
	category[["total"]] <-total
	ranged.data[[rangedID[n]]] <- category
}

save(ranged.data, file="./esets/summary/breast/ranged_pData_values_breast.rda")



