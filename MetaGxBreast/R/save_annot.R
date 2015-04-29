datasets <- read.csv("datasets.csv")
dataset.names <- datasets$Dataset
if(!file.exists("./annotations")){
  dir.create("./annotations")
}

for (i in 1:length(dataset.names)){
  dataset.name <- dataset.names[i]
  if(dataset.name == "METABRIC" | dataset.name == "TCGA"){
      next
  }
  load(paste("./data/", dataset.name,".RData", sep=""))
  save(annot, file= paste("./annotations/", dataset.name, "_annot.rda", sep=""))
}

