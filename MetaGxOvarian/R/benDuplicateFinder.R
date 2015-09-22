library(Biobase)

nthread<-1
availcore <- parallel::detectCores()
if (nthread > availcore) { nthread <- availcore }
options("mc.cores"=nthread)

filenames <- paste("./esets/mapped_esets/", list.files(path="./esets/mapped_esets/",pattern="*.rda"),sep="")
lapply(filenames, load, .GlobalEnv)

datasets <- read.csv("datasets.csv")
n <- which(datasets$Dataset == "TCGA-mirna-8x15kv2")
dataset.names <- datasets$Dataset[-n] ## ## not including TCGA.mirna.8x15kv2_eset

esets <- lapply(as.character(dataset.names), get)
names(esets) <- dataset.names

for(i in 1:length(esets)){
  message(i)
  sampleNames(esets[[i]]) <- paste(names(esets)[i], ":", colnames(exprs(esets[[i]])), sep="")
}
 #save(esets, file="./esets/esets.rda")
# load("./esets/esets.rda")
esets.mapped <- list()
for(eset in esets){
  Biobase::exprs(eset) <- exprs(eset)[fData(eset)$best_probe,]
  Biobase::fData(eset) <- fData(eset)[fData(eset)$best_probe,]
  rownames(fData(eset)) <- rownames(exprs(eset)) <- paste("geneid.", fData(eset)$EntrezGene.ID, sep="")
  esets.mapped <- c(esets.mapped, eset)
}
names(esets.mapped) <- names(esets)


source("./R/duplicateFinder.R")
source("./R/datasetMerging.R")

eset.merged <- datasetMerging(esets.mapped, method="union")
duplicates <- duplicateFinder(eset.merged, dupl.cor=0.98, method="spearman", nthread=nthread)
save(duplicates, file= "./esets/BenDuplicate.rda")


## list to be in package
remove <- duplicates
for(i in 1:length(remove)){
  remove[[i]] <- sort(c(names(remove[i]), remove[[i]]))
}
remove <- remove[!duplicated(remove)]
save(remove, file="./esets/removeSamples.rda")




### PLEASE NOTE!!!!!
## Example: A-B, and B-C pairs might be correlated at 0.98 (pass cutoff)
## BUT: This does not mean that A-C is correlated, because the CORRELATION of A-C might be below cutoff!
## MAKE SURE NOT TO CONFUSE THIS ASSOCIATION WITH THE IDEA OF THE MATHEMATICAL CORRELATION!!!!
# 
# ##for createEsetList
# for(i in 1:length(remove)){
#   datasetnames <- gsub(remove[[i]], pattern="\\..*$", replacement="")
#   if(remove.duplicates =="keep.smallest"){
#     remove[[i]] <- remove[[i]][-which.min(lapply(datasetnames, function(x){length(sampleNames(get(x)))}))]
#   } else if(remove.duplicates =="keep.largest"){
#     remove[[i]] <- remove[[i]][-which.max(lapply(datasetnames, function(x){length(sampleNames(get(x)))}))]
#   }
# }
# remove <- unique(gsub("\\.", replacement=":", x=unlist(remove)))

# 
# 
# #




###########################################################
## For use with levi's doppelgangR



# using Levi's
# duplicate_list <- list()
# # duplicate_list[[1]] <- NULL
# 
# already.classified <- list()
# ## number of lists in duplicates_list
# i<-0
# first <- NULL
# second <- NULL
# 
# for (n in 1:length(duplicates)){
#   first <- duplicates.df[n, "sample1"]
#   second <- duplicates.df[n, "sample2"]
#   if (!is.element(first, already.classified)){
#     ## first not element, second is -> add first
#     
#     if(is.element(second, already.classified)){
#       
#       for(j in 1:i){
#         if(is.element(second, duplicate_list[[j]])){
#           duplicate_list[[j]] <- c(duplicate_list[[j]], first)
#           already.classified <- c(already.classified, first)
#         }
#       }
#       
#     } else {
#       ## second element is not classified -> add both to new list
#       i <- i+1
#       duplicate_list[[i]] <- c(first, second)
#       already.classified <- c(already.classified, first, second)
#     }
#     
#     
#   } else {
#     ## first is classified -> classify second if needed
#     if(!is.element(second, already.classified)){
#       for(j in 1:i){
#         if(is.element(first, duplicate_list[[j]])){
#           duplicate_list[[j]] <- c(duplicate_list[[j]], second)
#           already.classified <- c(already.classified, second)
#         }
#       }
#     }
#     
#   }
#   
# }
# save(duplicate_list, file="./esets/duplicate_list.rda")


## remove specified 