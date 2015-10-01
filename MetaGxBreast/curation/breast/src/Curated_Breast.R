source("../../functions.R")

if (!file.exists("../curated"))
    dir.create("../curated")

##Read files named xyz1111.csv, xyz2222.csv, etc.
filenames <- list.files(path="../uncurated", pattern="*csv")

##Create list of data frame names without the ".csv" part
ds.names <- sub("\\.csv$", "", filenames)

###Load all files
for(ds.name in ds.names){
  infile <- file.path("../uncurated", paste(ds.name,".csv",sep=""))
  outfile <- file.path("../curated", paste(ds.name,"_curated.txt",sep=""))
  uncurated <- read.csv(infile,as.is=TRUE,row.names=1)
  ##initial creation of curated dataframe
  curated <- initialCuratedDF(rownames(uncurated), template.filename="template_breast.csv")

  if(ds.name == "TCGA"){
	  ##--------------------
	  ##start the mappings
	  ##--------------------
	
	  ##alt_sample_name
	  ##unique_patient_ID
	  ##er
	  ##pgr
	  ##her2
	  ##age
	  ##days_to_death
	  ##vital_status
	  ##sample_type
	  ##batch
	
	  ##alt_sample_name 
	  curated$alt_sample_name <- uncurated$bcr_patient_barcode
	
	  ##unique_patient_ID
	  curated$unique_patient_ID <- uncurated$bcr_patient_uuid
	
	  ##er
	  curated$er[uncurated$er_status_by_ihc == "Negative"] <- "negative"
	  curated$er[uncurated$er_status_by_ihc == "Positive"] <- "positive"
	  curated$er[uncurated$er_status_by_ihc == "[Not Evaluated]"] <- NA
	  curated$er[uncurated$er_status_by_ihc == "Indeterminate"] <- NA
	
	  ##pgr
	  curated$pgr[uncurated$pr_status_by_ihc == "Negative"] <- "negative"
	  curated$pgr[uncurated$pr_status_by_ihc == "Positive"] <- "positive"
	  curated$pgr[uncurated$pr_status_by_ihc == "[Not Evaluated]"] <- NA
	  curated$pgr[uncurated$pr_status_by_ihc == "Indeterminate"] <- NA
	
	
	  ##her2
	  curated$her2[uncurated$her2_status_by_ihc == "Negative"] <- "negative"
	  curated$her2[uncurated$her2_status_by_ihc == "Positive"] <- "positive"
	  curated$her2[uncurated$her2_status_by_ihc == "Equivolcal"] <- NA
	  curated$her2[uncurated$her2_status_by_ihc == "[Not Available]"] <- NA
	  curated$her2[uncurated$her2_status_by_ihc == "[Not Evaluated]"] <- NA
	  curated$her2[uncurated$her2_status_by_ihc == "Indeterminate"] <- NA
	
	  ##age_at_initial_pathologic_diagnosis
	  curated$age_at_initial_pathologic_diagnosis <- uncurated$age_at_diagnosis
	
	  ##vital_status
	  tmp <-uncurated$vital_status
	  tmp[tmp=="Alive"] <- "living"
	  tmp[tmp=="Dead"] <- "deceased"
	  curated$vital_status <- tmp
	  
	  ##days_to_death
	  curated[curated$vital_status=="living",]$days_to_death <- uncurated$last_contact_days_to[curated$vital_status=="living"]
	  curated[curated$vital_status=="deceased",]$days_to_death <- uncurated$death_days_to[curated$vital_status=="deceased"]
	  
	  ##tissue
	  tmp <- rep("tumor", each=length(curated$alt_sample_name))
	  curated$sample_type <- tmp
	  
	  ##batch
	  tmp <- rep("TCGA", each=length(curated$alt_sample_name))
	  curated$batch <- tmp
  } else {
	
	  ##--------------------
	  ##start the mappings
	  ##--------------------
	
	  ##alt_sample_name
	  ##unique_patient_ID
	  ##er
	  ##pgr
	  ##her2
	  ##size
	  ##node
	  ##age
	  ##grade
	  ##t.dmfs
	  ##e.dmfs
	  ##t.rfs
	  ##e.rfs
	  ##t.os
	  ##e.os
	  ##tissue
	  ##treatment
	  ##series
	
	
	  ##id -> alt_sample_name
	  ##Should be unique_patient_id for DFHCC2, MDA4,
	  ##For SUPERTAM_HGU133A.csv and TRANSBIG, ignore the id column
	  if(ds.name %in% c("DFHCC2", "MDA4")){
	      curated$unique_patient_ID <- uncurated$id
	  }else if(!ds.name %in% c("SUPERTAM_HGU133A", "TRANSBIG")) {
	      curated$alt_sample_name <- uncurated$id
	  }  #don't use uncurated$id at all for SUPERTAM_HGU133A and TRANSBIG.
	
	  ##er
	  curated$er[uncurated$er==0] <- "negative"
	  curated$er[uncurated$er==1] <- "positive"
	
	  ##pgr
	  curated$pgr[uncurated$pgr==0] <- "negative"
	  curated$pgr[uncurated$pgr==1] <- "positive"
	
	  ##her2
	  curated$her2[uncurated$her2==0] <- "negative"
	  curated$her2[uncurated$her2==1] <- "positive"
	
	  if(ds.name %in% c("GSE32646")){
	    curated$her2<-uncurated$her2
	    curated$er<-uncurated$er
	    curated$pgr<-uncurated$pgr
	  }
	  
	  if(ds.name %in% c("GSE32646","GSE25066","GSE48091","GSE58644")){
	    curated$tumor_size<-NA}
	  else{
	  ##tumor_size
	  tmp <- uncurated$size
	  tmp[tmp < 0] <- NA
	  curated$tumor_size <- tmp}
	
	  if(ds.name %in% c("GSE32646","GSE25066","GSE48091","GSE58644")){
	    curated$T<-uncurated$T}
	  
	  ##N
	  if(ds.name %in% c("GSE48091")){
	    curated$N <-NA}
	    else{
	  tmp <- uncurated$node
	  tmp[tmp > 1] <- 1
	  curated$N <- tmp}
	
	  ##age_at_initial_pathologic_diagnosis
	  curated$age_at_initial_pathologic_diagnosis <- uncurated$age
	
	  ##grade
	  curated$grade <- uncurated$grade
	
	  ##dmfs_days
	  curated$dmfs_days <- uncurated$t.dmfs
	
	  ##dmfs_status
	  curated$dmfs_status[uncurated$e.dmfs==0] <- "norecurrence"
	  curated$dmfs_status[uncurated$e.dmfs==1] <- "recurrence"
	
	  ##days_to_tumor_recurrence
	  curated$days_to_tumor_recurrence <- uncurated$t.rfs
	
	  ##recurrence_status
	  curated$recurrence_status[uncurated$e.rfs==0] <- "norecurrence"
	  curated$recurrence_status[uncurated$e.rfs==1] <- "recurrence"
	
	  ##days_to_death
	  curated$days_to_death <- uncurated$t.os
	
	  ##vital_status
	  curated$vital_status[uncurated$e.os==0] <- "living"
	  curated$vital_status[uncurated$e.os==1] <- "deceased"
	
	  ##asmple_type
	  curated$sample_type[uncurated$tissue==1] <- "tumor"
	  curated$sample_type[uncurated$tissue==0] <- "healthy"
	
	  ##treatment
	  curated$treatment[uncurated$treatment==0] <- "untreated"
	  curated$treatment[uncurated$treatment==1] <- "chemotherapy"
	  curated$treatment[uncurated$treatment==2] <- "hormonotherapy"
	  curated$treatment[uncurated$treatment==5] <- "endocrine"
	  curated$treatment[uncurated$treatment==6] <- "chemo.plus.hormono"
	  curated$treatment[uncurated$treatment==99] <- NA
	  
	  ##series
	  curated$batch <- uncurated$series
	
	  # For consistency of the definition of dmfs and rfs across datasets: for TRANSBIG, if the days_to_death and recurrence/distant metastasis are the same, then set recurrence/distant metastasis to "no event" (i.e. death is not an event)
	  if(ds.name == "TRANSBIG") {
	  	curated$dmfs_status[curated$dmfs_days == curated$days_to_death & curated$vital_status == "deceased" & curated$dmfs_status == "recurrence"] <- "norecurrence"
	  	curated$recurrence_status[curated$days_to_tumor_recurrence == curated$days_to_death & curated$vital_status == "deceased" & curated$recurrence_status == "recurrence"] <- "norecurrence"
	  }
  }
  curated <- postProcess(curated, uncurated)

  write.table(curated, file=outfile, row.names=FALSE, sep="\t")
  
}

