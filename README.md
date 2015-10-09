# UPDATED OCTOBER 9, 2015
# Gendoo et al

# MetaGxData Package Compendium

#######################################
VERSION CONTROL

V2.2 - Current Draft 

V2.3 - Modification to gene-wise and patient-wise normalization and new datasets added


#######################################
To build : 

Create tar.gz file: R CMD BUILD MetaGx______

To install:

R CMD INSTALL MetaGx_______


To get esets in data package:

library(MetaGx_____)

source(system.file("extdata", "patientselection.config", package="MetaGx_____"))

source(system.file("extdata", "createEsetList.R", package="MetaGx______"))




########################################


Currently manipulates data from A Three-Gene Model to Robustly Identify Breast Cancer Molecular Subtypes (http://compbio.dfci.harvard.edu/pubs/sbtpaper/data.zip)

Includes TCGA and METABRIC

Total number of expression sets: 39

MetaGxOvarian

Currently manipulates data from FULLVcuratedOvarianData (http://bcb.dfci.harvard.edu/ovariancancer/)

Includes TCGA

Total number of expression sets: 25
