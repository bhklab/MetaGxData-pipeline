(files <- dir(pattern="txt$"))

##change the index [1] on the next line to determine which file to check:
curated <- read.delim(files[1])

catvars.to.check <- c("er", "pgr", "her2", "N", "grade", "recurrence_status", "vital_status", "tissue", "treatment", "percent_normal_cells", "percent_stromal_cells", "percent_tumor_cells")
for (column in catvars.to.check)
    curated[[column]] <- factor(curated[[column]])

summary(curated[, -match("uncurated_author_metadata", colnames(curated))])

par(mfrow=c(4,3))
for (column in catvars.to.check){
    if(sum(!is.na(curated[[column]])) > 1)
        plot(factor(curated[[column]]), main=column)
}


##Make Kaplan-Meier plots if survival info is available.
library(survival)

(dmfs.events <- sum(curated$dmfs_status=="deceased_or_recurrence", na.rm=TRUE))
(recurrence.events <- sum(curated$recurrence_status=="deceased_or_recurrence", na.rm=TRUE))
(os.events <- sum(curated$vital_status=="deceased", na.rm=TRUE))


if(dmfs.events > 1){
    par(mfrow=c(1, 2))
    surv.obj <- Surv(curated$dmfs_days/30.4, curated$dmfs_status=="deceased_or_recurrence")
    surv.kmfit <- survfit(surv.obj ~ 1)
    plot(surv.kmfit, xlab="Time (months)", ylab="DMFS probability", main=files[1])
    cens.obj <- Surv(curated$dmfs_days/30.4, curated$dmfs_status=="living_norecurrence")
    cens.kmfit <- survfit(cens.obj ~ 1)
    plot(cens.kmfit, xlab="Time (months)", ylab="DMFS follow-up probability (reverse Kaplan-Meier estimate)", main=files[1])
}else{
    print("No DMFS available")
}


if(recurrence.events > 1){
    par(mfrow=c(1, 2))
    surv.obj <- Surv(curated$days_to_tumor_recurrence/30.4, curated$recurrence_status=="deceased_or_recurrence")
    surv.kmfit <- survfit(surv.obj ~ 1)
    plot(surv.kmfit, xlab="Time (months)", ylab="Recurrence probability", main=files[1])
    cens.obj <- Surv(curated$days_to_tumor_recurrence/30.4, curated$recurrence_status=="living_norecurrence")
    cens.kmfit <- survfit(cens.obj ~ 1)
    plot(cens.kmfit, xlab="Time (months)", ylab="Recurrence follow-up probability (reverse Kaplan-Meier estimate)", main=files[1])
}else{
    print("No recurrence available")
}


if(os.events > 1){
    par(mfrow=c(1, 2))
    surv.obj <- Surv(curated$days_to_death/30.4, curated$vital_status=="deceased")
    surv.kmfit <- survfit(surv.obj ~ 1)
    plot(surv.kmfit, xlab="Time (months)", ylab="Overall Survivalprobability", main=files[1])
    cens.obj <- Surv(curated$days_to_death/30.4, curated$vital_status=="living")
    cens.kmfit <- survfit(cens.obj ~ 1)
    plot(cens.kmfit, xlab="Time (months)", ylab="Overall Survival follow-up probability (reverse Kaplan-Meier estimate)", main=files[1])
}else{
    print("No overall survival available")
}

