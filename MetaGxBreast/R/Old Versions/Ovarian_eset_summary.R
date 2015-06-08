######################################
### OVARIAN SUMMARY
#####################################

setwd("~/Documents/curatedOvarianData/lwaldron-curatedovariandata/")

if(!file.exists("./esets/summary/ovarian")){
	dir.create("./esets/summary/ovarian", recursive=T)
}

##---------------------------
## load ovarian esets
##---------------------------
# library(logging)
# library(curatedOvarianData)
# source(system.file("extdata","patientselection.config",package="curatedOvarianData"))
# rm(min.number.of.events)
# rm(meta.required)
# rm(add.surv.y)
# rm(rule.1)
# source(system.file("extdata", "createEsetList.R", package ="curatedOvarianData"))

##---------------------------
## Number of samples
##---------------------------


## generate dataframe for summary of number of samples
numSamples <- NULL
for(i in 1:length(esets)){
	numSamples <- c(numSamples, length(sampleNames(esets[[i]])))
}

SampleNumberSummary <- data.frame(NumberOfSamples = numSamples, row.names = names(esets))
SampleNumberSummary <- rbind(SampleNumberSummary, sum(SampleNumberSummary[,"NumberOfSamples"]))

rownames(SampleNumberSummary)[nrow(SampleNumberSummary)] <- "Total"
SampleNumberSummaryPercent <- data.frame(Percent = SampleNumberSummary$NumberOfSamples/SampleNumberSummary$NumberOfSamples[nrow(SampleNumberSummary)])
rownames(SampleNumberSummaryPercent) <- rownames(SampleNumberSummary)
# save(SampleNumberSummary, file="./esets/summary/ovarian/SampleNumberSummary.rda")
# save(SampleNumberSummaryPercent, file="./esets/summary/ovarian/SampleNumberSummaryPercent.rda")
save(list = c(SampleNumberSummary, SampleNumberSummaryPercent), file = "./esets/summary/ovarian/SampleNumberSummaries_breast.rda")

##---------------------------
## Overall pData available
##---------------------------

pDataID <- c("sample_type", "histological_type", "primarysite", "summarygrade", "summarystage", "tumorstage", "grade", "age_at_initial_pathologic_diagnosis", "pltx", "tax", "neo", "days_to_tumor_recurrence", "recurrence_status", "days_to_death", "vital_status")


pDataPercentSummaryTable <- NULL
pDataSummaryNumbersTable <- NULL
for(e in 1:length(esets)){
	print(e)
	eset <- esets[[e]]
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
colnames(pDataSummaryNumbersTable) <- names(esets)
colnames(pDataPercentSummaryTable) <- names(esets)

pDataSummaryNumbersTable <- rbind(pDataSummaryNumbersTable, SampleNumberSummary["Total",])
rownames(pDataSummaryNumbersTable)[nrow(pDataSummaryNumbersTable)] <- "Total"

# save(pDataSummaryNumbersTable, file="./esets/summary/ovarian/pDataSummaryNumbersTable.rda")
# save(pDataPercentSummaryTable, file="./esets/summary/ovarian/pDataPercentSummaryTable.rda")
save(list=c(pDataSummaryNumbersTable, pDataPercentSummaryTable), file = "./esets/summary/ovarian/OverallpDataSummary_breast.rda")


##---------------------------
## Specific pData available
##---------------------------

## sample_type summary
sample_type_numbers <- NULL
sample_type_percent <- NULL
categories <- c("tumor", "metastatic", "cellline", "healthy", "adjacentnormal", "missing", "total")
for(i in 1:length(esets)){
	tumor <- length(which(as.character(pData(esets[[i]])$sample_type)=="tumor"))
	metastatic <- length(which(as.character(pData(esets[[i]])$sample_type)=="metastatic"))
	cellline <- length(which(as.character(pData(esets[[i]])$sample_type)=="cellline"))
	healthy <- length(which(as.character(pData(esets[[i]])$sample_type)=="healthy"))
	adjacentnormal <- length(which(as.character(pData(esets[[i]])$sample_type)=="adjacentnormal"))
	missing <- sum(is.na(pData(esets[[i]])$sample_type))
	total <- length(sampleNames(esets[[i]]))
	tmp <- c(tumor, metastatic, cellline, healthy, adjacentnormal, missing, total)
	tmpp <- tmp/tmp[length(categories)]*100
	if(i == 1){
		sample_type_numbers <- tmp
		sample_type_percent <- tmpp
	} else {
		sample_type_numbers <- rbind(sample_type_numbers, tmp)
		sample_type_percent <- rbind(sample_type_percent, tmpp)
	}
}
colnames(sample_type_numbers) <- colnames(sample_type_percent) <- categories
rownames(sample_type_numbers) <- rownames(sample_type_percent) <- names(esets)

sample_type_numbers <- t(sample_type_numbers)
sample_type_percent <- t(sample_type_percent)
total <- unlist(lapply(rownames(sample_type_numbers), function(x){sum(sample_type_numbers[x,])/sum(sample_type_numbers["total",])}))
total <- total*100
sample_type_percent <- cbind(sample_type_percent, total)
#save(sample_type_percent,file="./esets/summary/ovarian/sample_type_percent.rda")
#save(sample_type_numbers,file="./esets/summary/ovarian/sample_type_numbers.rda")

## histological_type summary
histological_type_numbers <- NULL
histological_type_percent <- NULL
categories <- c("ser", "endo", "clearcell", "mucinous", "other", "mix",  "undifferentiated", "missing", "total")
for(i in 1:length(esets)){
	ser <- length(which(as.character(pData(esets[[i]])$histological_type)=="ser"))
	endo <- length(which(as.character(pData(esets[[i]])$histological_type)=="endo"))
	clearcell <- length(which(as.character(pData(esets[[i]])$histological_type)=="clearcell"))
	mucinous <- length(which(as.character(pData(esets[[i]])$histological_type)=="mucinous"))
	other <- length(which(as.character(pData(esets[[i]])$histological_type)=="other"))
	mix <- length(which(as.character(pData(esets[[i]])$histological_type)=="mix"))
	undifferentiated <- length(which(as.character(pData(esets[[i]])$histological_type)=="undifferentiated"))
	missing <- sum(is.na(pData(esets[[i]])$histological_type))
	total <- length(sampleNames(esets[[i]]))
	tmp <- c(ser, endo, clearcell, mucinous, other, mix, undifferentiated, missing, total)
	tmpp <- tmp/tmp[length(categories)]*100
	if(i == 1){
		histological_type_numbers <- tmp
		histological_type_percent <- tmpp
	} else {
		histological_type_numbers <- rbind(histological_type_numbers, tmp)
		histological_type_percent <- rbind(histological_type_percent, tmpp)
	}
}
colnames(histological_type_numbers) <- colnames(histological_type_percent) <- categories
rownames(histological_type_numbers) <- rownames(histological_type_percent) <- names(esets)

histological_type_numbers <- t(histological_type_numbers)
histological_type_percent <- t(histological_type_percent)
total <- unlist(lapply(rownames(histological_type_numbers), function(x){sum(histological_type_numbers[x,])/sum(histological_type_numbers["total",])}))

total <- total*100
histological_type_percent <- cbind(histological_type_percent, total)
#save(histological_type_percent,file="./esets/summary/ovarian/histological_type_percent.rda")
#save(histological_type_numbers,file="./esets/summary/ovarian/histological_type_numbers.rda")



## primarysite summary
primarysite_numbers <- NULL
primarysite_percent <- NULL
categories <- c("ov", "ft", "other", "missing", "total")
for(i in 1:length(esets)){
	ov <- length(which(as.character(pData(esets[[i]])$primarysite)=="ov"))
	ft <- length(which(as.character(pData(esets[[i]])$primarysite)=="ft"))
	other <- length(which(as.character(pData(esets[[i]])$primarysite)=="other"))
	missing <- sum(is.na(pData(esets[[i]])$primarysite))
	total <- length(sampleNames(esets[[i]]))
	tmp <- c(ov, ft, other, missing, total)
	tmpp <- tmp/tmp[length(categories)]*100
	if(i == 1){
		primarysite_numbers <- tmp
		primarysite_percent <- tmpp

	} else {
		primarysite_numbers <- rbind(primarysite_numbers, tmp)
		primarysite_percent <- rbind(primarysite_percent, tmpp)
	}
}
colnames(primarysite_numbers) <- colnames(primarysite_percent) <- categories
rownames(primarysite_numbers) <- rownames(primarysite_percent) <- names(esets)
primarysite_numbers <- t(primarysite_numbers)
primarysite_percent <- t(primarysite_percent)
total <- unlist(lapply(rownames(primarysite_numbers), function(x){sum(primarysite_numbers[x,])/sum(primarysite_numbers["total",])}))
total <- total*100
primarysite_percent <- cbind(primarysite_percent, total)
#save(primarysite_percent,file="./esets/summary/ovarian/primarysite_percent.rda")
#save(primarysite_numbers,file="./esets/summary/ovarian/primarysite_numbers.rda")


## summarygrade summary

summarygrade_numbers <- NULL
summarygrade_percent <- NULL
categories <- c("low", "high", "missing", "total")
for(i in 1:length(esets)){
	low <- length(which(as.character(pData(esets[[i]])$summarygrade)=="low"))
	high <- length(which(as.character(pData(esets[[i]])$summarygrade)=="high"))
	missing <- sum(is.na(pData(esets[[i]])$summarygrade))
	total <- length(sampleNames(esets[[i]]))
	tmp <- c(low, high, missing, total)
	tmpp <- tmp/tmp[length(categories)]*100
	if(i == 1){
		summarygrade_numbers <- tmp
		summarygrade_percent <- tmpp

	} else {
		summarygrade_numbers <- rbind(summarygrade_numbers, tmp)
		summarygrade_percent <- rbind(summarygrade_percent, tmpp)
	}
}
colnames(summarygrade_numbers) <- colnames(summarygrade_percent) <- categories
rownames(summarygrade_numbers) <- rownames(summarygrade_percent) <- names(esets)
summarygrade_numbers <- t(summarygrade_numbers)
summarygrade_percent <- t(summarygrade_percent)
total <-  unlist(lapply(rownames(summarygrade_numbers), function(x){sum(summarygrade_numbers[x,])/sum(summarygrade_numbers["total",])}))
total <- total*100
summarygrade_percent <- cbind(summarygrade_percent, total)
#save(summarygrade_percent,file="./esets/summary/ovarian/summarygrade_percent.rda")
#save(summarygrade_numbers,file="./esets/summary/ovarian/summarygrade_numbers.rda")



## summarystage summary
summarystage_numbers <- NULL
summarystage_percent <- NULL
categories <- c("early", "late", "missing", "total")
for(i in 1:length(esets)){
	early <- length(which(as.character(pData(esets[[i]])$summarystage)=="early"))
	late <- length(which(as.character(pData(esets[[i]])$summarystage)=="late"))
	missing <- sum(is.na(pData(esets[[i]])$summarystage))
	total <- length(sampleNames(esets[[i]]))
	tmp <- c(early, late, missing, total)
	tmpp <- tmp/tmp[length(categories)]*100
	if(i == 1){
		summarystage_numbers <- tmp
		summarystage_percent <- tmpp

	} else {
		summarystage_numbers <- rbind(summarystage_numbers, tmp)
		summarystage_percent <- rbind(summarystage_percent, tmpp)
	}
}
colnames(summarystage_numbers) <- colnames(summarystage_percent) <- categories
rownames(summarystage_numbers) <- rownames(summarystage_percent) <- names(esets)
summarystage_numbers <- t(summarystage_numbers)
summarystage_percent <- t(summarystage_percent)
total <-  unlist(lapply(rownames(summarystage_numbers), function(x){sum(summarystage_numbers[x,])/sum(summarystage_numbers["total",])}))
total <- total*100
summarystage_percent <- cbind(summarystage_percent, total)
#save(summarystage_percent,file="./esets/summary/ovarian/summarystage_percent.rda")
#save(summarystage_numbers,file="./esets/summary/ovarian/summarystage_numbers.rda")


## tumorstage summary
tumorstage_numbers <- NULL
tumorstage_percent <- NULL
categories <- c("one", "two", "three", "four", "missing", "total")
for(i in 1:length(esets)){
	one <- length(which(as.character(pData(esets[[i]])$tumorstage)=="1"))
	two <- length(which(as.character(pData(esets[[i]])$tumorstage)=="2"))
	three <- length(which(as.character(pData(esets[[i]])$tumorstage)=="3"))
	four <- length(which(as.character(pData(esets[[i]])$tumorstage)=="4"))
	missing <- sum(is.na(pData(esets[[i]])$tumorstage))
	total <- length(sampleNames(esets[[i]]))
	tmp <- c(one, two, three, four, missing, total)
	tmpp <- tmp/tmp[length(categories)]*100
	if(i == 1){
		tumorstage_numbers <- tmp
		tumorstage_percent <- tmpp

	} else {
		tumorstage_numbers <- rbind(tumorstage_numbers, tmp)
		tumorstage_percent <- rbind(tumorstage_percent, tmpp)
	}
}
colnames(tumorstage_numbers) <- colnames(tumorstage_percent) <- categories
rownames(tumorstage_numbers) <- rownames(tumorstage_percent) <- names(esets)
tumorstage_numbers <- t(tumorstage_numbers)
tumorstage_percent <- t(tumorstage_percent)
total <-  unlist(lapply(rownames(tumorstage_numbers), function(x){sum(tumorstage_numbers[x,])/sum(tumorstage_numbers["total",])}))
total <- total*100
tumorstage_percent <- cbind(tumorstage_percent, total)
#save(tumorstage_percent,file="./esets/summary/ovarian/tumorstage_percent.rda")
#save(tumorstage_numbers,file="./esets/summary/ovarian/tumorstage_numbers.rda")


## grade summary
grade_numbers <- NULL
grade_percent <- NULL
categories <- c("one", "two", "three", "missing", "total")
for(i in 1:length(esets)){
	one <- length(which(as.character(pData(esets[[i]])$grade)=="1"))
	two <- length(which(as.character(pData(esets[[i]])$grade)=="2"))
	three <- length(which(as.character(pData(esets[[i]])$grade)=="3"))
	missing <- sum(is.na(pData(esets[[i]])$grade))
	total <- length(sampleNames(esets[[i]]))
	tmp <- c(one, two, three, missing, total)
	tmpp <- tmp/tmp[length(categories)]*100
	if(i == 1){
		grade_numbers <- tmp
		grade_percent <- tmpp

	} else {
		grade_numbers <- rbind(grade_numbers, tmp)
		grade_percent <- rbind(grade_percent, tmpp)
	}
}
colnames(grade_numbers) <- colnames(grade_percent) <- categories
rownames(grade_numbers) <- rownames(grade_percent) <- names(esets)
grade_numbers <- t(grade_numbers)
grade_percent <- t(grade_percent)
total <-  unlist(lapply(rownames(grade_numbers), function(x){sum(grade_numbers[x,])/sum(grade_numbers["total",])}))
total <- total*100
grade_percent <- cbind(grade_percent, total)
#save(grade_percent,file="./esets/summary/ovarian/grade_percent.rda")
#save(grade_numbers,file="./esets/summary/ovarian/grade_numbers.rda")

## pltx summary
pltx_numbers <- NULL
pltx_percent <- NULL
categories <- c("y", "n", "missing", "total")
for(i in 1:length(esets)){
	y <- length(which(as.character(pData(esets[[i]])$pltx)=="y"))
	n <- length(which(as.character(pData(esets[[i]])$pltx)=="n"))
	missing <- sum(is.na(pData(esets[[i]])$pltx))
	total <- length(sampleNames(esets[[i]]))
	tmp <- c(y, n,  missing, total)
	tmpp <- tmp/tmp[length(categories)]*100
	if(i == 1){
		pltx_numbers <- tmp
		pltx_percent <- tmpp

	} else {
		pltx_numbers <- rbind(pltx_numbers, tmp)
		pltx_percent <- rbind(pltx_percent, tmpp)
	}
}
colnames(pltx_numbers) <- colnames(pltx_percent) <- categories
rownames(pltx_numbers) <- rownames(pltx_percent) <- names(esets)
pltx_numbers <- t(pltx_numbers)
pltx_percent <- t(pltx_percent)
total <-  unlist(lapply(rownames(pltx_numbers), function(x){sum(pltx_numbers[x,])/sum(pltx_numbers["total",])}))
total <- total*100
pltx_percent <- cbind(pltx_percent, total)
#save(pltx_percent,file="./esets/summary/ovarian/pltx_percent.rda")
#save(pltx_numbers,file="./esets/summary/ovarian/pltx_numbers.rda")


## tax summary
tax_numbers <- NULL
tax_percent <- NULL
categories <- c("y", "n", "missing", "total")
for(i in 1:length(esets)){
	y <- length(which(as.character(pData(esets[[i]])$tax)=="y"))
	n <- length(which(as.character(pData(esets[[i]])$tax)=="n"))
	missing <- sum(is.na(pData(esets[[i]])$tax))
	total <- length(sampleNames(esets[[i]]))
	tmp <- c(y, n,  missing, total)
	tmpp <- tmp/tmp[length(categories)]*100
	if(i == 1){
		tax_numbers <- tmp
		tax_percent <- tmpp

	} else {
		tax_numbers <- rbind(tax_numbers, tmp)
		tax_percent <- rbind(tax_percent, tmpp)
	}
}
colnames(tax_numbers) <- colnames(tax_percent) <- categories
rownames(tax_numbers) <- rownames(tax_percent) <- names(esets)
tax_numbers <- t(tax_numbers)
tax_percent <- t(tax_percent)
total <-  unlist(lapply(rownames(tax_numbers), function(x){sum(tax_numbers[x,])/sum(tax_numbers["total",])}))
total <- total*100
tax_percent <- cbind(tax_percent, total)
#save(tax_percent,file="./esets/summary/ovarian/tax_percent.rda")
#save(tax_numbers,file="./esets/summary/ovarian/tax_numbers.rda")

## neo summary
neo_numbers <- NULL
neo_percent <- NULL
categories <- c("y", "n", "missing", "total")
for(i in 1:length(esets)){
	y <- length(which(as.character(pData(esets[[i]])$neo)=="y"))
	n <- length(which(as.character(pData(esets[[i]])$neo)=="n"))
	missing <- sum(is.na(pData(esets[[i]])$neo))
	total <- length(sampleNames(esets[[i]]))
	tmp <- c(y, n,  missing, total)
	tmpp <- tmp/tmp[length(categories)]*100
	if(i == 1){
		neo_numbers <- tmp
		neo_percent <- tmpp

	} else {
		neo_numbers <- rbind(neo_numbers, tmp)
		neo_percent <- rbind(neo_percent, tmpp)
	}
}
colnames(neo_numbers) <- colnames(neo_percent) <- categories
rownames(neo_numbers) <- rownames(neo_percent) <- names(esets)
neo_numbers <- t(neo_numbers)
neo_percent <- t(neo_percent)
total <-  unlist(lapply(rownames(neo_numbers), function(x){sum(neo_numbers[x,])/sum(neo_numbers["total",])}))
total <- total*100
neo_percent <- cbind(neo_percent, total)
#save(neo_percent,file="./esets/summary/ovarian/neo_percent.rda")
#save(neo_numbers,file="./esets/summary/ovarian/neo_numbers.rda")


## recurrence_status summary
recurrence_status_numbers <- NULL
recurrence_status_percent <- NULL
categories <- c("recurrence", "no recurrence" ,"missing", "total")
for(i in 1:length(esets)){
	recurrence <- length(which(as.character(pData(esets[[i]])$recurrence_status)=="recurrence"))
	no <- length(which(as.character(pData(esets[[i]])$recurrence_status)=="no recurrence"))
	missing <- sum(is.na(pData(esets[[i]])$recurrence_status))
	total <- length(sampleNames(esets[[i]]))
	tmp <- c(recurrence, no,  missing, total)
	tmpp <- tmp/tmp[length(categories)]*100
	if(i == 1){
		recurrence_status_numbers <- tmp
		recurrence_status_percent <- tmpp

	} else {
		recurrence_status_numbers <- rbind(recurrence_status_numbers, tmp)
		recurrence_status_percent <- rbind(recurrence_status_percent, tmpp)
	}
}
colnames(recurrence_status_numbers) <- colnames(recurrence_status_percent) <- categories
rownames(recurrence_status_numbers) <- rownames(recurrence_status_percent) <- names(esets)
recurrence_status_numbers <- t(recurrence_status_numbers)
recurrence_status_percent <- t(recurrence_status_percent)
total <-  unlist(lapply(rownames(recurrence_status_numbers), function(x){sum(recurrence_status_numbers[x,])/sum(recurrence_status_numbers["total",])}))
total <- total*100
recurrence_status_percent <- cbind(recurrence_status_percent, total)
#save(recurrence_status_percent,file="./esets/summary/ovarian/recurrence_status_percent.rda")
#save(recurrence_status_numbers,file="./esets/summary/ovarian/recurrence_status_numbers.rda")


## vital_status summary
vital_status_numbers <- NULL
vital_status_percent <- NULL
categories <- c("living", "deceased" ,"missing", "total")
for(i in 1:length(esets)){
	living <- length(which(as.character(pData(esets[[i]])$vital_status)=="living"))
	deceased <- length(which(as.character(pData(esets[[i]])$vital_status)=="deceased"))
	missing <- sum(is.na(pData(esets[[i]])$vital_status))
	total <- length(sampleNames(esets[[i]]))
	tmp <- c(living, deceased,  missing, total)
	tmpp <- tmp/tmp[length(categories)]*100
	if(i == 1){
		vital_status_numbers <- tmp
		vital_status_percent <- tmpp

	} else {
		vital_status_numbers <- rbind(vital_status_numbers, tmp)
		vital_status_percent <- rbind(vital_status_percent, tmpp)
	}
}
colnames(vital_status_numbers) <- colnames(vital_status_percent) <- categories
rownames(vital_status_numbers) <- rownames(vital_status_percent) <- names(esets)
vital_status_numbers <- t(vital_status_numbers)
vital_status_percent <- t(vital_status_percent)
total <-  unlist(lapply(rownames(vital_status_numbers), function(x){sum(vital_status_numbers[x,])/sum(vital_status_numbers["total",])}))
total <- total*100
vital_status_percent <- cbind(vital_status_percent, total)
#save(vital_status_percent,file="./esets/summary/ovarian/vital_status_percent.rda")
#save(vital_status_numbers,file="./esets/summary/ovarian/vital_status_numbers.rda")


save(list=paste(c("sample_type", "histological_type", "primarysite", "summarygrade", "summarystage", "tumorstage", "grade", "pltx", "tax", "neo", "recurrence_status", "vital_status"), "_percent", sep=""), file="./esets/summary/ovarian/ovarian_pData_categary_percent.rda")
save(list=paste(c("sample_type", "histological_type", "primarysite", "summarygrade", "summarystage", "tumorstage", "grade", "pltx", "tax", "neo", "recurrence_status", "vital_status"), "_numbers", sep=""), file="./esets/summary/ovarian/ovarian_pData_categary_numbers.rda")



##-------
## Ranged Data Summary
##------

rangedID <- c("age_at_initial_pathologic_diagnosis", "days_to_tumor_recurrence", "days_to_death")
ranged.data <- NULL
category <- NULL
eset.specific <- NULL

for(n in 1:length(rangedID)){
	total <- NULL
	for(i in 1:length(esets)){
		eset.specific <- pData(esets[[i]])[,rangedID[n]]
		category[[names(esets)[i]]] <- eset.specific
		total <- c(total, unlist(category[names(esets)[i]]))
	}
	category[["total"]] <-total
	ranged.data[[rangedID[n]]] <- category
}

save(ranged.data, file="./esets/summary/ovarian/ovarian_ranged_pData_values.rda")



