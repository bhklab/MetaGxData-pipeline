\name{ EORTC10994_eset }
\alias{ EORTC10994_eset }
\docType{data}
\title{  }
\usage{data( EORTC10994_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=15897907 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE1561 
  PMIDs: 15897907 
  No abstract available.
  notes:
   summary:     
      The tumors with an apocrine gene expression profile had strong histological apocrine features. These tumors were androgen
 receptor positive and were all ER negative, creating further classifications of tumor cells based on steroid receptor ac
tivity- luminal which are ER and AR positive, basal that are ER and AR negative, and molecular apocrine that are ER negat
ive and AR positive.
   version:     
      2015-04-27 19:00:24
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: 1007_s_at 1053_at ... AFFX-HUMISGF3A/M97935_MB_at (20967 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 20967 features, 49 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
EORTC10994_PF01 EORTC10994_PF02 EORTC10994_PF03 EORTC10994_PF04 EORTC10994_PF05 EORTC10994_PF06 EORTC10994_PF07 
              1               1               1               1               1               1               1 
EORTC10994_PF09 EORTC10994_PF10 EORTC10994_PF11 EORTC10994_PF12 EORTC10994_PF13 EORTC10994_PF14 EORTC10994_PF15 
              1               1               1               1               1               1               1 
EORTC10994_PF16 EORTC10994_PF17 EORTC10994_PF18 EORTC10994_PF19 EORTC10994_PF20 EORTC10994_PF21 EORTC10994_PF22 
              1               1               1               1               1               1               1 
EORTC10994_PF23 EORTC10994_PF24 EORTC10994_PF25 EORTC10994_PF26 EORTC10994_PF27 EORTC10994_PF28 EORTC10994_PF29 
              1               1               1               1               1               1               1 
EORTC10994_PF30 EORTC10994_PF31 EORTC10994_PF32 EORTC10994_PF33 EORTC10994_PF34 EORTC10994_PF35 EORTC10994_PF36 
              1               1               1               1               1               1               1 
EORTC10994_PF37 EORTC10994_PF38 EORTC10994_PF39 EORTC10994_PF40 EORTC10994_PF41 EORTC10994_PF42 EORTC10994_PF43 
              1               1               1               1               1               1               1 
EORTC10994_PF44 EORTC10994_PF45 EORTC10994_PF46 EORTC10994_PF47 EORTC10994_PF48 EORTC10994_PF49 EORTC10994_PF50 
              1               1               1               1               1               1               1 

alt_sample_name: 
PF01 EpPpT2N0G2 PF02 EpPpT4N1G2 PF03 EpPpT3N0Gu PF04 EpPpT2N1G1 PF05 EpPnT2N1G1 PF06 EnPnT2N0G3 PF07 EpPpT2N0G2 
              1               1               1               1               1               1               1 
PF09 EpPpT3N1G2 PF10 EpPpT4N1G2 PF11 EpPpT3N1G2 PF12 EpPpTiN0G2 PF13 EpPpT2N1G1 PF14 EnPnT2N1G2 PF15 EpPnTiN1G3 
              1               1               1               1               1               1               1 
PF16 EnPnT2N0G3 PF17 EpPnT2N1G3 PF18 EnPnT2N1G3 PF19 EnPuT4N0Gu PF20 EnPnT3N1G2 PF21 EpPpT2N1G2 PF22 EnPnT2N1G2 
              1               1               1               1               1               1               1 
PF23 EnPnT2N0G2 PF24 EnPnTiN0G3 PF25 EnPnT3N2G2 PF26 EpPnT3N0G3 PF27 EpPnT4N1G2 PF28 EnPnT2N1G3 PF29 EnPnT3N1G3 
              1               1               1               1               1               1               1 
PF30 EpPpT2N0G3 PF31 EnPnT2N0G3 PF32 EpPnT3N1G2 PF33 EpPnTiN0G2 PF34 EnPnT3N1G3 PF35 EpPpT2N1G3 PF36 EnPnT2N0G2 
              1               1               1               1               1               1               1 
PF37 EpPpT3N1G2 PF38 EnPnT2N1G3 PF39 EnPuT4N0Gu PF40 EpPnT4N0G2 PF41 EnPnT3N1G2 PF42 EnPnT2N2G3 PF43 EpPpT2N1G2 
              1               1               1               1               1               1               1 
PF44 EnPnT3N0G3 PF45 EpPpT4N0G2 PF46 EnPnT4N1G3 PF47 EpPpT3N1G3 PF48 EnPnT2N0G3 PF49 EpPpT2N1G1 PF50 EpPnT3N1G3 
              1               1               1               1               1               1               1 

sample_type: 
tumor 
   49 

er: 
negative positive 
      22       27 

pgr: 
negative positive     NA's 
      29       18        2 

tumor_size: 
 1  2  3  4 
 4 23 14  8 

N: 
 0  1 
19 30 

grade: 
   1    2    3 NA's 
   4   22   20    3 

batch: 
EORTC10994 
        49 

uncurated_author_metadata: 
  samplename: EORTC10994_PF01///dataset: EORTC10994///series: EORTC10994///id: PF01 EpPpT2N0G2///er: 1///pgr: 1///her2: NA///size: 2///node: 0///age: NA///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF02///dataset: EORTC10994///series: EORTC10994///id: PF02 EpPpT4N1G2///er: 1///pgr: 1///her2: NA///size: 4///node: 1///age: NA///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
 samplename: EORTC10994_PF03///dataset: EORTC10994///series: EORTC10994///id: PF03 EpPpT3N0Gu///er: 1///pgr: 1///her2: NA///size: 3///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF04///dataset: EORTC10994///series: EORTC10994///id: PF04 EpPpT2N1G1///er: 1///pgr: 1///her2: NA///size: 2///node: 1///age: NA///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF05///dataset: EORTC10994///series: EORTC10994///id: PF05 EpPnT2N1G1///er: 1///pgr: 0///her2: NA///size: 2///node: 1///age: NA///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF06///dataset: EORTC10994///series: EORTC10994///id: PF06 EnPnT2N0G3///er: 0///pgr: 0///her2: NA///size: 2///node: 0///age: NA///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF07///dataset: EORTC10994///series: EORTC10994///id: PF07 EpPpT2N0G2///er: 1///pgr: 1///her2: NA///size: 2///node: 0///age: NA///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF09///dataset: EORTC10994///series: EORTC10994///id: PF09 EpPpT3N1G2///er: 1///pgr: 1///her2: NA///size: 3///node: 1///age: NA///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF10///dataset: EORTC10994///series: EORTC10994///id: PF10 EpPpT4N1G2///er: 1///pgr: 1///her2: NA///size: 4///node: 1///age: NA///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF11///dataset: EORTC10994///series: EORTC10994///id: PF11 EpPpT3N1G2///er: 1///pgr: 1///her2: NA///size: 3///node: 1///age: NA///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF12///dataset: EORTC10994///series: EORTC10994///id: PF12 EpPpTiN0G2///er: 1///pgr: 1///her2: NA///size: 1///node: 0///age: NA///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF13///dataset: EORTC10994///series: EORTC10994///id: PF13 EpPpT2N1G1///er: 1///pgr: 1///her2: NA///size: 2///node: 1///age: NA///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF14///dataset: EORTC10994///series: EORTC10994///id: PF14 EnPnT2N1G2///er: 0///pgr: 0///her2: NA///size: 2///node: 1///age: NA///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF15///dataset: EORTC10994///series: EORTC10994///id: PF15 EpPnTiN1G3///er: 1///pgr: 0///her2: NA///size: 1///node: 1///age: NA///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF16///dataset: EORTC10994///series: EORTC10994///id: PF16 EnPnT2N0G3///er: 0///pgr: 0///her2: NA///size: 2///node: 0///age: NA///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF17///dataset: EORTC10994///series: EORTC10994///id: PF17 EpPnT2N1G3///er: 1///pgr: 0///her2: NA///size: 2///node: 1///age: NA///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF18///dataset: EORTC10994///series: EORTC10994///id: PF18 EnPnT2N1G3///er: 0///pgr: 0///her2: NA///size: 2///node: 1///age: NA///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
samplename: EORTC10994_PF19///dataset: EORTC10994///series: EORTC10994///id: PF19 EnPuT4N0Gu///er: 0///pgr: NA///her2: NA///size: 4///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF20///dataset: EORTC10994///series: EORTC10994///id: PF20 EnPnT3N1G2///er: 0///pgr: 0///her2: NA///size: 3///node: 1///age: NA///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF21///dataset: EORTC10994///series: EORTC10994///id: PF21 EpPpT2N1G2///er: 1///pgr: 1///her2: NA///size: 2///node: 1///age: NA///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF22///dataset: EORTC10994///series: EORTC10994///id: PF22 EnPnT2N1G2///er: 0///pgr: 0///her2: NA///size: 2///node: 1///age: NA///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF23///dataset: EORTC10994///series: EORTC10994///id: PF23 EnPnT2N0G2///er: 0///pgr: 0///her2: NA///size: 2///node: 0///age: NA///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF24///dataset: EORTC10994///series: EORTC10994///id: PF24 EnPnTiN0G3///er: 0///pgr: 0///her2: NA///size: 1///node: 0///age: NA///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF25///dataset: EORTC10994///series: EORTC10994///id: PF25 EnPnT3N2G2///er: 0///pgr: 0///her2: NA///size: 3///node: 1///age: NA///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF26///dataset: EORTC10994///series: EORTC10994///id: PF26 EpPnT3N0G3///er: 1///pgr: 0///her2: NA///size: 3///node: 0///age: NA///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF27///dataset: EORTC10994///series: EORTC10994///id: PF27 EpPnT4N1G2///er: 1///pgr: 0///her2: NA///size: 4///node: 1///age: NA///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF28///dataset: EORTC10994///series: EORTC10994///id: PF28 EnPnT2N1G3///er: 0///pgr: 0///her2: NA///size: 2///node: 1///age: NA///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF29///dataset: EORTC10994///series: EORTC10994///id: PF29 EnPnT3N1G3///er: 0///pgr: 0///her2: NA///size: 3///node: 1///age: NA///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF30///dataset: EORTC10994///series: EORTC10994///id: PF30 EpPpT2N0G3///er: 1///pgr: 1///her2: NA///size: 2///node: 0///age: NA///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF31///dataset: EORTC10994///series: EORTC10994///id: PF31 EnPnT2N0G3///er: 0///pgr: 0///her2: NA///size: 2///node: 0///age: NA///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF32///dataset: EORTC10994///series: EORTC10994///id: PF32 EpPnT3N1G2///er: 1///pgr: 0///her2: NA///size: 3///node: 1///age: NA///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF33///dataset: EORTC10994///series: EORTC10994///id: PF33 EpPnTiN0G2///er: 1///pgr: 0///her2: NA///size: 1///node: 0///age: NA///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF34///dataset: EORTC10994///series: EORTC10994///id: PF34 EnPnT3N1G3///er: 0///pgr: 0///her2: NA///size: 3///node: 1///age: NA///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF35///dataset: EORTC10994///series: EORTC10994///id: PF35 EpPpT2N1G3///er: 1///pgr: 1///her2: NA///size: 2///node: 1///age: NA///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF36///dataset: EORTC10994///series: EORTC10994///id: PF36 EnPnT2N0G2///er: 0///pgr: 0///her2: NA///size: 2///node: 0///age: NA///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF37///dataset: EORTC10994///series: EORTC10994///id: PF37 EpPpT3N1G2///er: 1///pgr: 1///her2: NA///size: 3///node: 1///age: NA///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF38///dataset: EORTC10994///series: EORTC10994///id: PF38 EnPnT2N1G3///er: 0///pgr: 0///her2: NA///size: 2///node: 1///age: NA///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
samplename: EORTC10994_PF39///dataset: EORTC10994///series: EORTC10994///id: PF39 EnPuT4N0Gu///er: 0///pgr: NA///her2: NA///size: 4///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF40///dataset: EORTC10994///series: EORTC10994///id: PF40 EpPnT4N0G2///er: 1///pgr: 0///her2: NA///size: 4///node: 0///age: NA///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF41///dataset: EORTC10994///series: EORTC10994///id: PF41 EnPnT3N1G2///er: 0///pgr: 0///her2: NA///size: 3///node: 1///age: NA///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF42///dataset: EORTC10994///series: EORTC10994///id: PF42 EnPnT2N2G3///er: 0///pgr: 0///her2: NA///size: 2///node: 1///age: NA///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF43///dataset: EORTC10994///series: EORTC10994///id: PF43 EpPpT2N1G2///er: 1///pgr: 1///her2: NA///size: 2///node: 1///age: NA///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF44///dataset: EORTC10994///series: EORTC10994///id: PF44 EnPnT3N0G3///er: 0///pgr: 0///her2: NA///size: 3///node: 0///age: NA///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF45///dataset: EORTC10994///series: EORTC10994///id: PF45 EpPpT4N0G2///er: 1///pgr: 1///her2: NA///size: 4///node: 0///age: NA///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF46///dataset: EORTC10994///series: EORTC10994///id: PF46 EnPnT4N1G3///er: 0///pgr: 0///her2: NA///size: 4///node: 1///age: NA///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF47///dataset: EORTC10994///series: EORTC10994///id: PF47 EpPpT3N1G3///er: 1///pgr: 1///her2: NA///size: 3///node: 1///age: NA///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF48///dataset: EORTC10994///series: EORTC10994///id: PF48 EnPnT2N0G3///er: 0///pgr: 0///her2: NA///size: 2///node: 0///age: NA///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF49///dataset: EORTC10994///series: EORTC10994///id: PF49 EpPpT2N1G1///er: 1///pgr: 1///her2: NA///size: 2///node: 1///age: NA///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 
  samplename: EORTC10994_PF50///dataset: EORTC10994///series: EORTC10994///id: PF50 EpPnT3N1G3///er: 1///pgr: 0///her2: NA///size: 3///node: 1///age: NA///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                                      1 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE1561 }
\keyword{datasets}
