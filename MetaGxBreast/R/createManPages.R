
eset.dir <- "./esets/mapped_esets2/"


cleanText <- function(x){
  gsub("%", "\\%", iconv(x, "latin1", "ASCII", sub="?"), fixed=TRUE)
}

eset.files <- dir(eset.dir,pattern="^.*\\.rda$")

for (iFile in 1:length(eset.files)){
  load(paste(eset.dir,eset.files[iFile],sep="/"))  #load the eset
  eset <- get(sub("_eset.rda","",eset.files[iFile],fixed=TRUE))  #call it eset
  accessionID <- sub("_eset.rda","",eset.files[iFile])
  eset.name <- sub(".rda","",eset.files[iFile])
  print(accessionID)
  ##  eset$batch <- factor(eset$batch)    #Make sure this is a factor ##ovarian specific
  do.call(rm, list(sub("_eset.rda","", eset.files[iFile], fixed=TRUE)))  #remove the original
  ##pdata.nonblank will contain pdata columns with any non-missing values:
  pdata.nonblank <- pData(eset)
  pdata.nonblank <- pdata.nonblank[,apply(pdata.nonblank,2,function(x) sum(!is.na(x)) > 0)]
  thisfile <- sub(".rda",".Rd",eset.files[iFile])
  sink(file=paste(eset.dir, thisfile, sep="/"))
  cat(paste("\\name{", eset.name, "}"))
  cat("\n")
  cat(paste("\\alias{", eset.name, "}"))
  cat("\n")
  cat(paste("\\docType{data}"))
  cat("\n")
  cat(paste("\\title{", cleanText(experimentData(eset)@title), "}"))
  cat("\n")
  if (abstract(eset) != ""){
    cat(paste("\\description{", cleanText(experimentData(eset)@abstract), "}"))
    cat("\n")
  }
  cat(paste("\\usage{data(", cleanText(eset.name), ")}"))
  cat("\n")
  cat("\\format{")
  cat("\n")
  cat("\\preformatted{")
  cat("\n")
  cat("experimentData(eset):")
  cat("\n")
  print( experimentData(eset) )
  cat("\n")
#   cat(paste("Preprocessing:", experimentData(eset)@preprocessing[[1]]))
#   cat("\n")
  cat("featureData(eset):")
  cat("\n")
  print( featureData(eset) )
  cat("\n")
  cat("}}")
  cat("\n")
  cat("\\details{")
  cat("\n")
  cat("\\preformatted{")
  cat("\n")
  cat(paste("assayData:", nrow(eset),"features,",ncol(eset),"samples"))
  cat("\n")
  cat(paste("Platform type:", eset@annotation))
  cat("\n")
  if(!all(is.na(eset$vital_status))){
    time <- eset$days_to_death / 365
    cens <- ifelse(eset$vital_status=="deceased",1,0)
    library(survival)
    cat("Overall survival time-to-event summary (in years):")
    cat("\n")
    print(survfit(Surv(time,cens)~-1))
    cat("\n")
  }
#   if(!all(is.na(eset$os_binary))){
#     cat("Binary overall survival summary (definitions of long and short provided by study authors): \n")
#     cat("\n")
#     print(summary(factor(eset$os_binary)))
#     cat("\n")
#   }
  cat( "--------------------------- \n")
  cat( "Available sample meta-data: \n")
  cat( "--------------------------- \n")
  cat( "\n")
  for (iCol in 1:ncol(pdata.nonblank)){
    if(length(unique(pdata.nonblank[,iCol])) < 6)
      pdata.nonblank[,iCol] <- factor(pdata.nonblank[,iCol])
    cat(paste(colnames(pdata.nonblank)[iCol],": \n",sep=""))
    print(summary(pdata.nonblank[,iCol]))
    cat( "\n")
  }
  cat("}}")
  cat("\n")
  if(experimentData(eset)@url != ""){
    cat(paste("\\source{", experimentData(eset)@url, "}"))
    cat("\n")
  }
  cat("\\keyword{datasets}")
  cat("\n")
  sink(NULL)
}