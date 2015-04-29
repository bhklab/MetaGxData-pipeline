url <- "http://compbio.dfci.harvard.edu/pubs/sbtpaper/data.zip"
data.dir <- "./data"

if(file.exists("datafiles.txt"))
  readdata.files <- read.table("datafiles.txt", as.is=TRUE)[, 1]

if(!exists("readdata.files") || !identical(readdata.files, dir(data.dir))){
  ## put data.zip in a temporary directory
  set.seed(1)
  tmp.dir <- tempdir()
  dest.file <- paste(tmp.dir, "data.zip", sep="/")
  if(!exists(dest.file))
    download.file(url=url, destfile=dest.file)
  ## unzip to ./data:
  unzip(dest.file, exdir=".")
  ## make a record of the unzipped files:
  data.files <- dir(data.dir)
  write.table(data.files, file="datafiles.txt", row.names=FALSE, col.names=FALSE)
  readdata.files <- read.table("datafiles.txt")[, 1]
}

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

