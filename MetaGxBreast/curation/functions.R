mapstage <- function(x){
  output <- x
  output[output=="I"] <- 1
  output[output=="II"] <- 2
  output[output=="III"] <- 3
  output[output=="IV"] <- 4
  output <- as.integer(output)
  return(output)
}

initialCuratedDF <- function(DF.rownames,template.filename){
  template <- read.csv(template.filename,as.is=TRUE)
  output <- matrix(NA,
                   ncol=nrow(template),
                   nrow=length(DF.rownames))
  colnames(output) <- template$col.name
  rownames(output) <- DF.rownames
  output <- data.frame(output)
  for (i in 1:ncol(output)){
    class(output[,i]) <- template[i,"var.class"]
  }
  output$sample_name <- DF.rownames
  return(output)
}

postProcess <- function(curated, uncurated=NULL, celfile.colname=NULL, do.celfile.batch=TRUE){
    ##curated is the dataframe of curated metadata, uncurated=uncurated.
    ##If celfile.colname is specified, this column will be used for
    ##matching to celfiles, instead of rownames(curated)
    if(!is.null(uncurated)){
        if( identical(all.equal(rownames(curated), rownames(uncurated)), TRUE) ){
            output.matrix <- sapply(1:ncol(uncurated), function(i){
                paste(colnames(uncurated)[i], uncurated[, i], sep=": ")
            })
            uncurated.original <- apply(output.matrix, 1, paste, collapse="///")
            ##strip any non-ASCII characters:
            Encoding(uncurated.original) <- "latin1"
            uncurated.ascii <- iconv(uncurated.original, "latin1", "ASCII", sub="")
            curated$uncurated_author_metadata <- uncurated.ascii
        }else{
            stop("rownames(curated) should be identical to rownames(uncurated)")
        }
    }
    if(exists("celfile.dir") & do.celfile.batch){
        if(file.exists(celfile.dir)){
            if(require(affyio)){
                fnames <- dir(celfile.dir,pattern="\\.[cC][eE][lL]\\.gz$|\\.[cC][eE][lL]$")
                if(length(fnames) < 1) return(curated)  ##break if there are no cel files
                print("Looking up celfile run dates...")
                celRunDate <-
                    sapply(fnames,function(fname)
                       {
                           tempDate1 <- read.celfile.header(paste(celfile.dir,fname,sep="/"),info="full")$ScanDate
                           output <- as.character(as.Date(tempDate1, format = "%m/%d/%y %H:%M:%S"))
                           if(length(tempDate1) > 0 && is.na(output))
                               output <- as.character(as.Date(tempDate1))
                           if(length(output)==0) output <- NA
                           return(output)
                       }
                           )
                ##get rid of the CEL file extensions:
                ##If this is a GEO series, then all filenames should start with GSM:
                isGEO <- all(grepl("gsm",names(celRunDate),ignore.case=TRUE))
                ##If it is a GEO series, get rid of extraneous characters after GSM[0-9]+
                if(isGEO){
                    ##Keep only GSM followed by some digits:
                    names(celRunDate) <- sub("(^GSM[0-9]+).*","\\1",names(celRunDate))
                }else{
                    ##get rid of the .CEL extension:
                    names(celRunDate) <- sub("\\.[cC][eE][lL]\\.gz$|\\.[cC][eE][lL]$","",names(celRunDate))
                }
                ##In some cases (specifically TCGA), the CEL file names are
                ##not the sample names, but can be found in another column of
                ##the curated dataframe.  celfile.colname allows you to
                ##specify this column for matching celfile run dates, rather
                ##than the rownames:
                if(identical(celfile.colname %in% colnames(curated),TRUE)){
                    sample.match.names <- curated[,celfile.colname]
                }else{
                    sample.match.names <- rownames(curated)
                }
                if(any(names(celRunDate) %in% sample.match.names)){
                    print("Adding celfile run dates")
                    celRunDate <- celRunDate[names(celRunDate) %in% sample.match.names]
                    curated[match(names(celRunDate),sample.match.names),"batch"] <- celRunDate
                }else{
                    warning("Could not match any celfile names to curated pdata names.")
                }
            }else{
                warning("affyio package not installed, so not attempting to add celfile run dates")
            }
        }
    }
    return(curated)
}
