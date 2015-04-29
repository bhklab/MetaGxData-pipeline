library(doppelgangR)

filenames <- paste("./esets/mapped_esets/", list.files("./esets/mapped_esets"),sep="")
lapply(filenames, load, .GlobalEnv)

datasets <- read.csv("datasets.csv")
dataset.names <- datasets$Dataset
esets <- lapply(as.character(dataset.names), get)
names(esets) <- dataset.names
duplicates <- doppelgangR(esets, outlierFinder.expr.args=list(bonf.prob=0.2, transFun=atanh, tail="upper"), phenoFinder.args = NULL, cache.dir=NULL)
save(duplicates, file ="./esets/duplicates.rda")
duplicates.df <- data.frame(summary(duplicates))

duplicate_list <- list()
duplicate_list[[1]] <- NULL

already.classified <- list()
## number of lists in duplicates_list
i<-1
first <- NULL
second <- NULL

for (n in 1:nrow(duplicates.df)){
  first <- duplicates.df[n, "sample1"]
  second <- duplicates.df[n, "sample2"]
  if (!is.element(first, already.classified)){
    ## first not element, second is -> add first
    
    if(is.element(second, already.classified)){
      
      for(j in 1:i){
        if(is.element(second, duplicate_list[[j]])){
          duplicate_list[[j]] <- c(duplicate_list[[j]], first)
          already.classified <- c(already.classified, first)
        }
      }
      
    } else {
      ## second element is not classified -> add both to new list
      i <- i+1
      duplicate_list[[i]] <- c(first, second)
      already.classified <- c(already.classified, first, second)
    }
    
    
  } else {
    ## first is classified -> classify second if needed
    if(!is.element(second, already.classified)){
      for(j in 1:i){
        if(is.element(first, duplicate_list[[j]])){
          duplicate_list[[j]] <- c(duplicate_list[[j]], second)
          already.classified <- c(already.classified, second)
        }
      }
    }
    
  }
  
}
save(duplicate_list, file="./esets/duplicate_list.rda")
