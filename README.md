# MetaGxData

To build : 

Create tar.gz file: R CMD BUILD MetaGx______

To install:

R CMD INSTALL MetaGx_______


To get esets in data package:

library(MetaGx_____)

source("extdata", "patientselection.config", package="MetaGx_____")

source("extdata", "createEsetList.R", package="MetaGx______")


Currently manipulates data from A Three-Gene Model to Robustly Identify Breast Cancer Molecular Subtypes (http://compbio.dfci.harvard.edu/pubs/sbtpaper/data.zip)

Includes TCGA and METABRIC

Total number of expression sets: 37

MetaGxOvarian

Currently manipulates data from FULLVcuratedOvarianData (http://bcb.dfci.harvard.edu/ovariancancer/)

Includes TCGA

Total number of expression sets: 22