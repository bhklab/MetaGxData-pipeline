library(GSA)
GO_BP<-GSA.read.gmt(filename="c5.bp.v5.0.entrez.gmt") #825 genesets
CanonicalPathways<-GSA.read.gmt(filename="c2.cp.v5.0.entrez.gmt") #1330 genesets


names(GO_BP)
#     > names(GO_BP)
#     [1] "genesets"             "geneset.names"        "geneset.descriptions"

names(GO_BP$genesets)<-GO_BP$geneset.names
# NOW HAVE A LIST OF GENE ONTOLOGY BIOLOGICAL PROCESSES
# Ex: GO_BP$genesets[1]

#     > GO_BP$genesets[1]
#     $TRNA_PROCESSING
#     [1] "23536" "51095" "10667" "4234"  "6301"  "16"    "54974" "6741"  "10775" "54888"

CompleteLists<-GO_BP$genesets #List of Genesets with EntrezIDs
