originalDir <- getwd()

##similar to gse_RAW.R from Levi
datasets <- read.csv("./datasets.csv")
GSEdatasets <- grep("GSE", datasets$Dataset, value=T)

### download raw 
if(!file.exists("./DATA")){
  dir.create("./DATA")
}

############################
## Downloading GSE datasets
############################

downloadGSE <- function(gse){
  require(GEOquery)
  x <- try(getGEOSuppFiles(gse, makeDirectory = FALSE,baseDir=celDir))
  if(class(x)=="try-error"){
    warning(paste(gse,"has no supplementary data"))
    file.remove(targetdir)
    stop(paste("failed to download supplemental data for",strInputAccession))
  } else{
    setwd(celDir)
    system("tar xf *.tar")
    system("rm *.tar")
    setwd(originalDir)
    return("Success")
  }
}


for (i in 1:length(GSEdatasets)){
	gse <- GSEdatasets[i]
	celDir <- paste("./DATA/", gse, "/RAW", sep="")
	if(file.exists(celDir)){
		warning(paste(celDir, "alreadyexists. Remove to re-download raw data."))
		next
	} else {
		dir.create(celDir, recursive=T)
	}
	options(download.file.method="wget")
	downloadGSE(gse)

}

############################
## Downloading PMID
############################
options(download.file.method="wget")
## PMID17290060


URL <- "https://discovery.genome.duke.edu/express/resources/1144/PlatinumJCO.zip"
study <- "PMID17290060"

rawDir <- paste("./DATA/", study, "/RAW")

defaultDir <- paste("./DATA/", study, "/PROCESSED/DEFAULT")

if(file.exists(rawDir)){
		warning(paste(rawDir, "alreadyexists. Remove to re-download raw data."))
		next
	} else {
		dir.create(rawDir, recursive=T)
		dir.create(defaultDir, recursive=T)
	}
download.file(url=URL, destfile=rawDir)
unzip(paste(rawDir, "/PlatinumJCO.zip"))



## PMID15897565

URL <- "https://discovery.genome.duke.edu/express/resources/1144/survival_CEL_files.zip"
study <-"PMID15897565"

rawDir <- paste("./DATA/", study, "/RAW")

defaultDir <- paste("./DATA/", study, "/PROCESSED/DEFAULT")

if(file.exists(rawDir)){
		warning(paste(rawDir, "alreadyexists. Remove to re-download raw data."))
		next
	} else {
		dir.create(rawDir, recursive=T)
		dir.create(defaultDir, recursive=T)
	}
download.file(url=URL, destfile=rawDir)
unzip(paste(rawDir, "/survival_CEL_files.zip"))



## PMID19318476

URL <- "https://discovery.genome.duke.edu/express/resources/74/suppmat.zip"
study <-"PMID19318476"

rawDir <- paste("./DATA/", study, "/RAW")

defaultDir <- paste("./DATA/", study, "/PROCESSED/DEFAULT")

if(file.exists(rawDir)){
		warning(paste(rawDir, "alreadyexists. Remove to re-download raw data."))
		next
	} else {
		dir.create(rawDir, recursive=T)
		dir.create(defaultDir, recursive=T)
	}

setwd(rawDir)
shellcommand <- "wget --no-directories --timestamping $URL
unzip -o suppmat.zip
rm suppmat.zip
mv suppmat/* .
mv data/batch1/*.cel .
mv data/batch2/*.cel .
mv data/batch3/*.cel .
mv data/late/*.cel .
rm -rf suppmat data

cp batches.txt ./DATA/PMID19318476/PROCESSED/DEFAULT/PMID19318476_batches.txt
cp NewClin.txt ./Data/PMID19318476/PROCESSED/DEFAULT/PMID19318476_full_pdata.txt
find `pwd` -name "*.cel" -fprint "$defaultDir"/"$strInputAccession"_RAWfilenames.txt"
system(shellcommand)










