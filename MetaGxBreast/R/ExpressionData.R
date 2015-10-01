library(GEOquery)

url <- "http://compbio.dfci.harvard.edu/pubs/sbtpaper/data.zip"
data.dir <- "./data"

##################################################
##### FOR FIRST 35 Gene Expression Datasets (From BEN) #####
##################################################

if(file.exists("datafiles.txt"))
  readdata.files <- read.table("datafiles.txt", as.is=TRUE)[, 1]

if(!exists("readdata.files") || !identical(readdata.files, dir(data.dir))){
  ## put data.zip in a temporary directory
  set.seed(1)
  tmp.dir <- tempdir()
  dest.file <- paste(tmp.dir, "data.zip", sep="/")
  if(!exists(dest.file))
    download.file(url=url, destfile=dest.file) #Obtain RData per Expression Dataset
  ## unzip to ./data:
  unzip(dest.file, exdir=".")
  ## make a record of the unzipped files:
  data.files <- dir(data.dir)
  write.table(data.files, file="datafiles.txt", row.names=FALSE, col.names=FALSE)
  readdata.files <- read.table("datafiles.txt")[, 1]
}

##################################################
##### FOR OTHER DATASETS OBTAINED FROM GEO #####
##################################################
# For these datasets, there is now zipped folder containing the RData
# First download the expression as an RData and save in a GEO sub-folder

if(file.exists("datafilesGEO.txt"))
  readdata.files2 <- read.table("datafilesGEO.txt", as.is=TRUE)[, 1]

if(!exists("readdata.files2") || !identical(readdata.files2, dir(data.dir))){
  ## create files and the RData per GSE to ./data:
  source("./R/ObtainGEOdatasets.R")
  
  ## make a record of the files:
  data.files <- dir(data.dir) #Now contains list of all files (JNCI + GEO)
  write.table(data.files, file="datafiles.txt", row.names=FALSE, col.names=FALSE)
  readdata.files <- read.table("datafiles.txt")[, 1]
}

##################################################
##### PROCESS ALL DATA FROM ALL SOURCES #####
#### Create all subdirectories per data file ####
##################################################

# readdata.files <- read.table("datafiles.txt")[,1]
dataset.names <-sub("^([^.]*).*", "\\1", readdata.files)

for (i in 1:length(dataset.names)){
  print(paste("Dataset:", dataset.names[i]))
  load(file.path(data.dir, readdata.files[i]))
  if(is(data, "function")) next
  path <- file.path(getwd(), "DATA", paste(dataset.names)[i], "PROCESSED","DEFAULT")
  dir.create(path, recursive=TRUE)
  out.filename <- file.path(path, paste0(dataset.names[i], "_default_curated_exprs.txt"))
  mydata<- t(data)
  write.table(mydata, file=out.filename, sep="\t")
  rm(data, demo, annot, path, out.filename, mydata)
}

