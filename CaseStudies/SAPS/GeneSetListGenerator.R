library(GSA)
hallmark.genesets<-GSA.read.gmt(filename="h.all.v5.0.entrez.gmt")
cancer.gene.neighbourhoods<-GSA.read.gmt(filename="c4.cgn.v5.0.entrez.gmt")

names(hallmark.genesets$genesets)<-hallmark.genesets$geneset.names
names(cancer.gene.neighbourhoods$genesets)<-cancer.gene.neighbourhoods$geneset.names

hallmark.genesets <- hallmark.genesets$genesets
cancer.gene.neighbourhoods <- cancer.gene.neighbourhoods$genesets
