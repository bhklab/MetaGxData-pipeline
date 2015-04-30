\name{ LUND_eset }
\alias{ LUND_eset }
\docType{data}
\title{  }
\usage{data( LUND_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=18430221 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31863 
  PMIDs: 18430221 
  No abstract available.
  notes:
   summary:     
      A significant difference was found between the ER positive subgroup and ER negative subgroup in the gene expression profi
les.
   version:     
      2015-04-27 19:13:54
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: H200006618 H200006808 ... H300022925 (11154 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 11154 features, 143 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
  LUND_CA10015_1_HP   LUND_CA10017_1_HP   LUND_CA10152_1_HP   LUND_CA10227_1_HP   LUND_CA10372_1_HP   LUND_CA10388_1_HP 
                  1                   1                   1                   1                   1                   1 
  LUND_CA10389_1_HP   LUND_CA10430_1_HP   LUND_CA10533_1_HP   LUND_CA10624_1_HP   LUND_CA10656_1_HP   LUND_CA10692_2_HP 
                  1                   1                   1                   1                   1                   1 
  LUND_CA10742_1_HP   LUND_CA10844_1_HP   LUND_CA10959_1_HP   LUND_CA11014_1_HP   LUND_CA11305_1_HP   LUND_CA11336_1_HP 
                  1                   1                   1                   1                   1                   1 
  LUND_CA11338_1_HP   LUND_CA11370_1_HP   LUND_CA12078_2_HP   LUND_CA12455_1_HP   LUND_CA12545_2_HP   LUND_CA12873_2_HP 
                  1                   1                   1                   1                   1                   1 
  LUND_CA12984_1_HP   LUND_CA13327_2_HP   LUND_CA13335_1_HP   LUND_CA13833_1_HP   LUND_CA13880_2_HP   LUND_CA13943_1_HP 
                  1                   1                   1                   1                   1                   1 
  LUND_CA14113_2_HP    LUND_CA6513_2_HP    LUND_CA6524_1_HP    LUND_CA6858_1_HP    LUND_CA7478_1_HP    LUND_CA8532_1_HP 
                  1                   1                   1                   1                   1                   1 
   LUND_CA8706_1_HP    LUND_CA8731_1_HP    LUND_CA8732_1_HP    LUND_CA8767_1_HP    LUND_CA8798_1_HP    LUND_CA8825_1_HP 
                  1                   1                   1                   1                   1                   1 
   LUND_CA8835_1_HP    LUND_CA8837_1_HP    LUND_CA8854_1_HP    LUND_CA8869_1_HP    LUND_CA8883_1_HP    LUND_CA8888_1_HP 
                  1                   1                   1                   1                   1                   1 
   LUND_CA8950_1_HP    LUND_CA8962_1_HP    LUND_CA8967_1_HP    LUND_CA9019_1_HP    LUND_CA9047_1_HP    LUND_CA9107_1_HP 
                  1                   1                   1                   1                   1                   1 
   LUND_CA9114_1_HP    LUND_CA9152_1_HP    LUND_CA9189_1_HP    LUND_CA9212_1_HP    LUND_CA9222_2_HP    LUND_CA9254_1_HP 
                  1                   1                   1                   1                   1                   1 
   LUND_CA9340_1_HP    LUND_CA9362_1_HP    LUND_CA9367_2_HP    LUND_CA9383_1_HP    LUND_CA9397_1_HP    LUND_CA9426_1_HP 
                  1                   1                   1                   1                   1                   1 
   LUND_CA9447_1_HP    LUND_CA9546_1_HP    LUND_CA9598_1_HP    LUND_CA9612_1_HP    LUND_CA9665_1_HP    LUND_CA9689_1_HP 
                  1                   1                   1                   1                   1                   1 
   LUND_CA9745_1_HP    LUND_CA9865_1_HP    LUND_CA9866_1_HP    LUND_CA9878_1_HP    LUND_CA9879_1_HP    LUND_CA9893_1_HP 
                  1                   1                   1                   1                   1                   1 
   LUND_CA9913_1_HP    LUND_CA9916_1_HP    LUND_CA9952_1_HP  LUND_G11_5029_1_HP LUND_G115_6106_1_HP LUND_G128_6107_1_HP 
                  1                   1                   1                   1                   1                   1 
LUND_G139_6255_1_HP LUND_G141_8535_1_HP LUND_G143_6316_1_HP LUND_G146_6330_2_HP LUND_G169_7944_2_HP LUND_G178_6720_1_HP 
                  1                   1                   1                   1                   1                   1 
LUND_G182_6710_1_HP LUND_G183_9454_1_HP LUND_G185_7434_2_HP LUND_G191_6749_1_HP LUND_G195_6528_1_HP LUND_G202_6877_1_HP 
                  1                   1                   1                   1                   1                   1 
LUND_G215_7046_1_HP LUND_G227_7174_2_HP LUND_G249_7474_2_HP             (Other) 
                  1                   1                   1                  44 

alt_sample_name: 
  CA10015   CA10017   CA10152   CA10227   CA10372   CA10388   CA10389   CA10430   CA10533   CA10624   CA10656   CA10692 
        1         1         1         1         1         1         1         1         1         1         1         1 
  CA10742   CA10844   CA10959   CA11014   CA11305   CA11336   CA11338   CA11370   CA12078   CA12455   CA12545   CA12873 
        1         1         1         1         1         1         1         1         1         1         1         1 
  CA12984   CA13327   CA13335   CA13833   CA13880   CA13943   CA14113    CA6513    CA6524    CA6858    CA7478    CA8532 
        1         1         1         1         1         1         1         1         1         1         1         1 
   CA8706    CA8731    CA8732    CA8767    CA8798    CA8825    CA8835    CA8837    CA8854    CA8869    CA8883    CA8888 
        1         1         1         1         1         1         1         1         1         1         1         1 
   CA8950    CA8962    CA8967    CA9019    CA9047    CA9107    CA9114    CA9152    CA9189    CA9212    CA9222    CA9254 
        1         1         1         1         1         1         1         1         1         1         1         1 
   CA9340    CA9362    CA9367    CA9383    CA9397    CA9426    CA9447    CA9546    CA9598    CA9612    CA9665    CA9689 
        1         1         1         1         1         1         1         1         1         1         1         1 
   CA9745    CA9865    CA9866    CA9878    CA9879    CA9893    CA9913    CA9916    CA9952  G11_5029 G115_6106 G128_6107 
        1         1         1         1         1         1         1         1         1         1         1         1 
G139_6255 G141_8535 G143_6316 G146_6330 G169_7944 G178_6720 G182_6710 G183_9454 G185_7434 G191_6749 G195_6528 G202_6877 
        1         1         1         1         1         1         1         1         1         1         1         1 
G215_7046 G227_7174 G249_7474   (Other) 
        1         1         1        44 

sample_type: 
tumor 
  143 

er: 
negative positive 
      29      114 

pgr: 
negative positive     NA's 
      47       88        8 

tumor_size: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  0.200   1.100   1.500   1.486   1.800   4.000       2 

N: 
  0 
143 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  27.00   47.50   56.00   54.76   63.00   73.00 

batch: 
LUNDS1 LUNDS2 LUNDS3 LUNDS4 
    30     47     22     44 

uncurated_author_metadata: 
    samplename: LUND_CA10015_1_HP///dataset: LUND///series: LUNDS4///id: CA10015///er: 1///pgr: 1///her2: NA///size: 1.8///node: 0///age: 49///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
    samplename: LUND_CA10017_1_HP///dataset: LUND///series: LUNDS2///id: CA10017///er: 1///pgr: 1///her2: NA///size: 1.7///node: 0///age: 57///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
    samplename: LUND_CA10152_1_HP///dataset: LUND///series: LUNDS4///id: CA10152///er: 1///pgr: 0///her2: NA///size: 1.2///node: 0///age: 63///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
    samplename: LUND_CA10227_1_HP///dataset: LUND///series: LUNDS4///id: CA10227///er: 1///pgr: 0///her2: NA///size: 1.3///node: 0///age: 54///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
    samplename: LUND_CA10372_1_HP///dataset: LUND///series: LUNDS2///id: CA10372///er: 1///pgr: 1///her2: NA///size: 1.1///node: 0///age: 72///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
    samplename: LUND_CA10388_1_HP///dataset: LUND///series: LUNDS4///id: CA10388///er: 1///pgr: 1///her2: NA///size: 1.5///node: 0///age: 61///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
    samplename: LUND_CA10389_1_HP///dataset: LUND///series: LUNDS2///id: CA10389///er: 1///pgr: 1///her2: NA///size: 0.8///node: 0///age: 64///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA10430_1_HP///dataset: LUND///series: LUNDS2///id: CA10430///er: 1///pgr: 1///her2: NA///size: 1///node: 0///age: 65///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
    samplename: LUND_CA10533_1_HP///dataset: LUND///series: LUNDS4///id: CA10533///er: 1///pgr: 1///her2: NA///size: 1.1///node: 0///age: 67///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
    samplename: LUND_CA10624_1_HP///dataset: LUND///series: LUNDS4///id: CA10624///er: 1///pgr: 1///her2: NA///size: 1.5///node: 0///age: 50///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
    samplename: LUND_CA10656_1_HP///dataset: LUND///series: LUNDS4///id: CA10656///er: 0///pgr: 0///her2: NA///size: 1.8///node: 0///age: 63///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
    samplename: LUND_CA10692_2_HP///dataset: LUND///series: LUNDS4///id: CA10692///er: 0///pgr: 0///her2: NA///size: 1.5///node: 0///age: 53///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
    samplename: LUND_CA10742_1_HP///dataset: LUND///series: LUNDS3///id: CA10742///er: 0///pgr: 1///her2: NA///size: 1.2///node: 0///age: 60///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA10844_1_HP///dataset: LUND///series: LUNDS2///id: CA10844///er: 1///pgr: 1///her2: NA///size: 2///node: 0///age: 57///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA10959_1_HP///dataset: LUND///series: LUNDS2///id: CA10959///er: 1///pgr: 0///her2: NA///size: 1///node: 0///age: 73///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
    samplename: LUND_CA11014_1_HP///dataset: LUND///series: LUNDS4///id: CA11014///er: 1///pgr: 1///her2: NA///size: 1.8///node: 0///age: 59///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
    samplename: LUND_CA11305_1_HP///dataset: LUND///series: LUNDS1///id: CA11305///er: 1///pgr: 1///her2: NA///size: 1.6///node: 0///age: 60///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA11336_1_HP///dataset: LUND///series: LUNDS4///id: CA11336///er: 1///pgr: 1///her2: NA///size: 1///node: 0///age: 49///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
    samplename: LUND_CA11338_1_HP///dataset: LUND///series: LUNDS3///id: CA11338///er: 0///pgr: 1///her2: NA///size: 1.4///node: 0///age: 46///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
    samplename: LUND_CA11370_1_HP///dataset: LUND///series: LUNDS3///id: CA11370///er: 1///pgr: 1///her2: NA///size: 1.5///node: 0///age: 44///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA12078_2_HP///dataset: LUND///series: LUNDS4///id: CA12078///er: 1///pgr: 1///her2: NA///size: 2///node: 0///age: 47///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
    samplename: LUND_CA12455_1_HP///dataset: LUND///series: LUNDS3///id: CA12455///er: 0///pgr: 0///her2: NA///size: 1.4///node: 0///age: 61///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
    samplename: LUND_CA12545_2_HP///dataset: LUND///series: LUNDS4///id: CA12545///er: 1///pgr: 1///her2: NA///size: 1.2///node: 0///age: 56///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
    samplename: LUND_CA12873_2_HP///dataset: LUND///series: LUNDS1///id: CA12873///er: 1///pgr: 0///her2: NA///size: 1.3///node: 0///age: 55///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
    samplename: LUND_CA12984_1_HP///dataset: LUND///series: LUNDS4///id: CA12984///er: 1///pgr: 1///her2: NA///size: 2.5///node: 0///age: 45///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
    samplename: LUND_CA13327_2_HP///dataset: LUND///series: LUNDS4///id: CA13327///er: 1///pgr: 1///her2: NA///size: 1.6///node: 0///age: 62///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
    samplename: LUND_CA13335_1_HP///dataset: LUND///series: LUNDS4///id: CA13335///er: 1///pgr: 1///her2: NA///size: 1.1///node: 0///age: 59///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
    samplename: LUND_CA13833_1_HP///dataset: LUND///series: LUNDS4///id: CA13833///er: 1///pgr: 1///her2: NA///size: 1.6///node: 0///age: 46///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
    samplename: LUND_CA13880_2_HP///dataset: LUND///series: LUNDS4///id: CA13880///er: 1///pgr: 1///her2: NA///size: 1.2///node: 0///age: 51///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA13943_1_HP///dataset: LUND///series: LUNDS1///id: CA13943///er: 1///pgr: 1///her2: NA///size: 2///node: 0///age: 49///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA14113_2_HP///dataset: LUND///series: LUNDS4///id: CA14113///er: 1///pgr: 1///her2: NA///size: 2///node: 0///age: 40///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA6513_2_HP///dataset: LUND///series: LUNDS1///id: CA6513///er: 0///pgr: 0///her2: NA///size: 1.6///node: 0///age: 50///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA6524_1_HP///dataset: LUND///series: LUNDS1///id: CA6524///er: 0///pgr: 0///her2: NA///size: 1.8///node: 0///age: 57///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
        samplename: LUND_CA6858_1_HP///dataset: LUND///series: LUNDS1///id: CA6858///er: 1///pgr: 1///her2: NA///size: 1///node: 0///age: 37///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
        samplename: LUND_CA7478_1_HP///dataset: LUND///series: LUNDS1///id: CA7478///er: 1///pgr: 1///her2: NA///size: 1///node: 0///age: 42///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA8532_1_HP///dataset: LUND///series: LUNDS1///id: CA8532///er: 0///pgr: 0///her2: NA///size: 1.5///node: 0///age: 50///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA8706_1_HP///dataset: LUND///series: LUNDS2///id: CA8706///er: 0///pgr: 0///her2: NA///size: 2.2///node: 0///age: 49///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA8731_1_HP///dataset: LUND///series: LUNDS2///id: CA8731///er: 1///pgr: 0///her2: NA///size: 0.8///node: 0///age: 68///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA8732_1_HP///dataset: LUND///series: LUNDS4///id: CA8732///er: 1///pgr: 1///her2: NA///size: 1.2///node: 0///age: 59///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA8767_1_HP///dataset: LUND///series: LUNDS2///id: CA8767///er: 1///pgr: 0///her2: NA///size: 0.4///node: 0///age: 65///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA8798_1_HP///dataset: LUND///series: LUNDS4///id: CA8798///er: 1///pgr: 1///her2: NA///size: 1.3///node: 0///age: 61///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA8825_1_HP///dataset: LUND///series: LUNDS4///id: CA8825///er: 1///pgr: 0///her2: NA///size: 1.8///node: 0///age: 66///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA8835_1_HP///dataset: LUND///series: LUNDS2///id: CA8835///er: 0///pgr: 0///her2: NA///size: 1.6///node: 0///age: 60///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
        samplename: LUND_CA8837_1_HP///dataset: LUND///series: LUNDS4///id: CA8837///er: 1///pgr: 1///her2: NA///size: 1///node: 0///age: 55///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA8854_1_HP///dataset: LUND///series: LUNDS4///id: CA8854///er: 0///pgr: 0///her2: NA///size: 1.2///node: 0///age: 61///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
        samplename: LUND_CA8869_1_HP///dataset: LUND///series: LUNDS4///id: CA8869///er: 0///pgr: 0///her2: NA///size: 2///node: 0///age: 52///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA8883_1_HP///dataset: LUND///series: LUNDS2///id: CA8883///er: 1///pgr: 1///her2: NA///size: 1.2///node: 0///age: 42///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA8888_1_HP///dataset: LUND///series: LUNDS3///id: CA8888///er: 1///pgr: 1///her2: NA///size: 0.9///node: 0///age: 73///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA8950_1_HP///dataset: LUND///series: LUNDS2///id: CA8950///er: 1///pgr: 1///her2: NA///size: 1.4///node: 0///age: 64///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA8962_1_HP///dataset: LUND///series: LUNDS4///id: CA8962///er: 1///pgr: 1///her2: NA///size: 0.7///node: 0///age: 70///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA8967_1_HP///dataset: LUND///series: LUNDS4///id: CA8967///er: 1///pgr: 0///her2: NA///size: 1.3///node: 0///age: 55///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
        samplename: LUND_CA9019_1_HP///dataset: LUND///series: LUNDS2///id: CA9019///er: 1///pgr: 1///her2: NA///size: 1///node: 0///age: 67///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
        samplename: LUND_CA9047_1_HP///dataset: LUND///series: LUNDS2///id: CA9047///er: 0///pgr: 0///her2: NA///size: 2///node: 0///age: 43///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
        samplename: LUND_CA9107_1_HP///dataset: LUND///series: LUNDS4///id: CA9107///er: 1///pgr: 1///her2: NA///size: 1///node: 0///age: 67///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA9114_1_HP///dataset: LUND///series: LUNDS4///id: CA9114///er: 1///pgr: 1///her2: NA///size: 1.5///node: 0///age: 66///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
        samplename: LUND_CA9152_1_HP///dataset: LUND///series: LUNDS3///id: CA9152///er: 1///pgr: 1///her2: NA///size: 1///node: 0///age: 45///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA9189_1_HP///dataset: LUND///series: LUNDS2///id: CA9189///er: 1///pgr: 1///her2: NA///size: 1.3///node: 0///age: 63///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
        samplename: LUND_CA9212_1_HP///dataset: LUND///series: LUNDS2///id: CA9212///er: 1///pgr: 1///her2: NA///size: 1///node: 0///age: 43///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
        samplename: LUND_CA9222_2_HP///dataset: LUND///series: LUNDS4///id: CA9222///er: 1///pgr: 0///her2: NA///size: 2///node: 0///age: 66///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA9254_1_HP///dataset: LUND///series: LUNDS4///id: CA9254///er: 1///pgr: 1///her2: NA///size: 1.5///node: 0///age: 55///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA9340_1_HP///dataset: LUND///series: LUNDS3///id: CA9340///er: 1///pgr: 1///her2: NA///size: 2.9///node: 0///age: 52///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
        samplename: LUND_CA9362_1_HP///dataset: LUND///series: LUNDS2///id: CA9362///er: 1///pgr: 1///her2: NA///size: 2///node: 0///age: 64///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA9367_2_HP///dataset: LUND///series: LUNDS3///id: CA9367///er: 1///pgr: 1///her2: NA///size: 1.5///node: 0///age: 54///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA9383_1_HP///dataset: LUND///series: LUNDS2///id: CA9383///er: 1///pgr: 1///her2: NA///size: 1.8///node: 0///age: 55///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA9397_1_HP///dataset: LUND///series: LUNDS2///id: CA9397///er: 1///pgr: 1///her2: NA///size: 1.7///node: 0///age: 68///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA9426_1_HP///dataset: LUND///series: LUNDS2///id: CA9426///er: 0///pgr: 0///her2: NA///size: 0.9///node: 0///age: 69///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
        samplename: LUND_CA9447_1_HP///dataset: LUND///series: LUNDS2///id: CA9447///er: 1///pgr: 1///her2: NA///size: 2///node: 0///age: 67///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA9546_1_HP///dataset: LUND///series: LUNDS3///id: CA9546///er: 1///pgr: 0///her2: NA///size: 1.2///node: 0///age: 62///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA9598_1_HP///dataset: LUND///series: LUNDS2///id: CA9598///er: 1///pgr: 1///her2: NA///size: 1.5///node: 0///age: 66///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA9612_1_HP///dataset: LUND///series: LUNDS1///id: CA9612///er: 0///pgr: 0///her2: NA///size: 1.7///node: 0///age: 49///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA9665_1_HP///dataset: LUND///series: LUNDS2///id: CA9665///er: 1///pgr: 1///her2: NA///size: 1.7///node: 0///age: 56///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
        samplename: LUND_CA9689_1_HP///dataset: LUND///series: LUNDS2///id: CA9689///er: 1///pgr: 1///her2: NA///size: 2///node: 0///age: 56///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA9745_1_HP///dataset: LUND///series: LUNDS2///id: CA9745///er: 1///pgr: 0///her2: NA///size: 1.7///node: 0///age: 68///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA9865_1_HP///dataset: LUND///series: LUNDS2///id: CA9865///er: 1///pgr: 1///her2: NA///size: 1.5///node: 0///age: 67///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA9866_1_HP///dataset: LUND///series: LUNDS2///id: CA9866///er: 1///pgr: 1///her2: NA///size: 1.3///node: 0///age: 50///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA9878_1_HP///dataset: LUND///series: LUNDS4///id: CA9878///er: 1///pgr: 1///her2: NA///size: 2.6///node: 0///age: 46///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA9879_1_HP///dataset: LUND///series: LUNDS2///id: CA9879///er: 1///pgr: 1///her2: NA///size: 1.5///node: 0///age: 54///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
        samplename: LUND_CA9893_1_HP///dataset: LUND///series: LUNDS2///id: CA9893///er: 1///pgr: 1///her2: NA///size: 1///node: 0///age: 65///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA9913_1_HP///dataset: LUND///series: LUNDS1///id: CA9913///er: 0///pgr: 0///her2: NA///size: 0.2///node: 0///age: 37///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA9916_1_HP///dataset: LUND///series: LUNDS2///id: CA9916///er: 1///pgr: 1///her2: NA///size: 1.2///node: 0///age: 57///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
      samplename: LUND_CA9952_1_HP///dataset: LUND///series: LUNDS4///id: CA9952///er: 1///pgr: 1///her2: NA///size: 1.1///node: 0///age: 61///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
  samplename: LUND_G11_5029_1_HP///dataset: LUND///series: LUNDS2///id: G11_5029///er: 1///pgr: 0///her2: NA///size: 0.5///node: 0///age: 61///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
samplename: LUND_G115_6106_1_HP///dataset: LUND///series: LUNDS4///id: G115_6106///er: 0///pgr: 0///her2: NA///size: 1.3///node: 0///age: 55///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
samplename: LUND_G128_6107_1_HP///dataset: LUND///series: LUNDS4///id: G128_6107///er: 1///pgr: 1///her2: NA///size: 1.1///node: 0///age: 52///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
samplename: LUND_G139_6255_1_HP///dataset: LUND///series: LUNDS4///id: G139_6255///er: 1///pgr: 1///her2: NA///size: 1.2///node: 0///age: 47///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
samplename: LUND_G141_8535_1_HP///dataset: LUND///series: LUNDS3///id: G141_8535///er: 1///pgr: 1///her2: NA///size: 1.5///node: 0///age: 48///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
  samplename: LUND_G143_6316_1_HP///dataset: LUND///series: LUNDS3///id: G143_6316///er: 1///pgr: 1///her2: NA///size: 2///node: 0///age: 58///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
samplename: LUND_G146_6330_2_HP///dataset: LUND///series: LUNDS3///id: G146_6330///er: 1///pgr: 1///her2: NA///size: 1.5///node: 0///age: 46///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
  samplename: LUND_G169_7944_2_HP///dataset: LUND///series: LUNDS3///id: G169_7944///er: 0///pgr: 0///her2: NA///size: 2///node: 0///age: 63///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
samplename: LUND_G178_6720_1_HP///dataset: LUND///series: LUNDS3///id: G178_6720///er: 1///pgr: 1///her2: NA///size: 1.7///node: 0///age: 49///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
samplename: LUND_G182_6710_1_HP///dataset: LUND///series: LUNDS2///id: G182_6710///er: 1///pgr: 0///her2: NA///size: 1.4///node: 0///age: 64///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
samplename: LUND_G183_9454_1_HP///dataset: LUND///series: LUNDS3///id: G183_9454///er: 1///pgr: 1///her2: NA///size: 0.7///node: 0///age: 67///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
samplename: LUND_G185_7434_2_HP///dataset: LUND///series: LUNDS3///id: G185_7434///er: 1///pgr: 1///her2: NA///size: 1.8///node: 0///age: 45///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
samplename: LUND_G191_6749_1_HP///dataset: LUND///series: LUNDS2///id: G191_6749///er: 1///pgr: 1///her2: NA///size: 1.3///node: 0///age: 44///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
samplename: LUND_G195_6528_1_HP///dataset: LUND///series: LUNDS2///id: G195_6528///er: 1///pgr: 1///her2: NA///size: 2.2///node: 0///age: 57///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
  samplename: LUND_G202_6877_1_HP///dataset: LUND///series: LUNDS3///id: G202_6877///er: 0///pgr: 0///her2: NA///size: 3///node: 0///age: 59///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
  samplename: LUND_G215_7046_1_HP///dataset: LUND///series: LUNDS4///id: G215_7046///er: 1///pgr: 1///her2: NA///size: 4///node: 0///age: 65///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
samplename: LUND_G227_7174_2_HP///dataset: LUND///series: LUNDS4///id: G227_7174///er: 1///pgr: 0///her2: NA///size: 1.8///node: 0///age: 67///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
samplename: LUND_G249_7474_2_HP///dataset: LUND///series: LUNDS2///id: G249_7474///er: 1///pgr: 1///her2: NA///size: 2.1///node: 0///age: 62///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                           1 
                                                                                                                                                                                                                                                     (Other) 
                                                                                                                                                                                                                                                          44 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31863 }
\keyword{datasets}
