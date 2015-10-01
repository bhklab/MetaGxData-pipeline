# UPDATED OCTOBER 1, 2015
# Deena M.A. Gendoo

# MetaGxData Package Compendium

#######################################
VERSION CONTROL
V2.1 - Current Draft 
V2.2 - Modification to gene-wise and patient-wise normalization
V2.3 - Same as 2.2 but 4 new datasets added

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

Total number of expression sets: 37

MetaGxOvarian

Currently manipulates data from FULLVcuratedOvarianData (http://bcb.dfci.harvard.edu/ovariancancer/)

Includes TCGA

Total number of expression sets: 22
