\name{ UNC4_eset }
\alias{ UNC4_eset }
\docType{data}
\title{ UNC4_eset }
\usage{data( UNC4_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=16141321 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse3494 
  PMIDs: 16141321 
  No abstract available.
  notes:
   summary:     
      A 32-gene expression signature of p53 was identified that differentiates p-53 mutant and wild-type tumors. The signature 
is more effective than sequence-based assessments of p53 in predicting prognosis and therapeutic response.  
   version:     
      2015-04-27 19:46:53
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: probe.10 probe.12 ... probe.79701 (5420 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 5420 features, 305 samples
Platform type: 
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

   65 observations deleted due to missingness 
records   n.max n.start  events  median 0.95LCL 0.95UCL 
 240.00  240.00  240.00   51.00    7.73    6.82      NA 

--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   a44k10829.1192    a44k10829.1281    a44k10829.1282    a44k10829.1283    a44k11240.2741    a44k11240.2760    a44k11240.2761 
                1                 1                 1                 1                 1                 1                 1 
      h1a.013na.c        h1a.c.0077        h1a.c.0078        h1a.c.0079        h1a.c.0080        h1a.c.0081       h1a.v2.1313 
                1                 1                 1                 1                 1                 1                 1 
      h1a.v2.1314     h1a.v2.2411na       h1a0120na.c        h1a012na.c        h1a014na.c         h1a2371v2      h1a4110.3579 
                1                 1                 1                 1                 1                 1                 1 
     h1a4110.3580      h1a4110.3582         h1a4594v2      h1a4956v2.na       h1a5365.899       h1a5365.900       h1a5949.019 
                1                 1                 1                 1                 1                 1                 1 
      h1a6149.098       h1a6149.698       h1a6260.731       h1a6260.732    h1a6261.0922v2    h1a6261.0924v2    h1a6261.0925v2 
                1                 1                 1                 1                 1                 1                 1 
   h1a6261.1488v2    h1a6261.1489v2 h1a6261.1489v2.na    h1a6261.1490v2    h1a6261.1499v2    h1a6261.1500v2    h1a6261.2087v2 
                1                 1                 1                 1                 1                 1                 1 
   h1a6261.2140v2    h1a6261.2141v2 h1a6261.2597v2.na h1a6261.2598v2.na    h1a6261.2607v2 h1a6261.2980v2.na    h1a6433.2272v2 
                1                 1                 1                 1                 1                 1                 1 
  h1a6433.2273.v2   h1a6433.2274.v2   h1a6433.2275.v2   h1a6433.2276.v2    h1a6433.2597v2    h1a6433.2598v2    h1a6433.2599v2 
                1                 1                 1                 1                 1                 1                 1 
   h1a6433.2600v2    h1a6433.2601v2    h1a6433.2608v2 h1a6433.2608v2.na    h1a6433.2610v2 h1a6433.2610v2.na    h1a6433.2611v2 
                1                 1                 1                 1                 1                 1                 1 
  h1a6433.2757.v2   h1a6433.2758.v2   h1a6433.2759.v2   h1a6433.2761.v2   h1a6433.2808.v2   h1a6433.2809.v2   h1a6433.2810.v2 
                1                 1                 1                 1                 1                 1                 1 
  h1a6433.2811.v2   h1a6433.2812.v2   h1a6433.2978.v2   h1a6433.2979.v2    h1a6433.2987v2    h1a6433.2988v2    h1a6433.2989v2 
                1                 1                 1                 1                 1                 1                 1 
   h1a6433.2990v2    h1a6433.2992v2   h1a6433.3007.v2   h1a6433.3009.v2   h1a6433.3010.v2   h1a6433.3011.v2       h1a6546.126 
                1                 1                 1                 1                 1                 1                 1 
      h1a6730.644       h1a6730.735       h1a6730.744       h1a6730.745    h1a6912.6851v2    h1a6912.6852v2    h1a6912.6853v2 
                1                 1                 1                 1                 1                 1                 1 
   h1a6912.6955v2    h1a6912.6956v2    h1a6912.6957v2    h1a6912.6958v2    h1a6912.6984v2    h1a6912.6985v2    h1a6912.6986v2 
                1                 1                 1                 1                 1                 1                 1 
   h1a6912.6987v2           (Other) 
                1               206 

sample_type: 
tumor 
  305 

er: 
negative positive     NA's 
      99      154       52 

pgr: 
negative positive     NA's 
     126      109       70 

her2: 
negative positive     NA's 
     203       58       44 

tumor_size: 
   1  1.5    3    6 NA's 
  60    1  129   43   72 

N: 
   0    1 NA's 
 126  135   44 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  24.00   46.00   55.00   56.73   68.00   89.00      59 

grade: 
   1    2    3 NA's 
  25   80  138   62 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
   30.0   450.0   750.0   954.3  1380.0  3540.0      64 

recurrence_status: 
deceased_or_recurrence    living_norecurrence                   NA's 
                    70                    170                     65 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
     30     540     885    1104    1590    5190      65 

vital_status: 
deceased   living     NA's 
      51      189       65 

batch: 
UNC4 
 305 

uncurated_author_metadata: 
        samplename: a44k10829.1192///dataset: UNC4///series: UNC4///id: NA///er: 1///pgr: 0///her2: 0///size: 6///node: 0///age: 69///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 420///e.rfs: 1///t.os: 420///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
       samplename: a44k10829.1281///dataset: UNC4///series: UNC4///id: NA///er: 1///pgr: 1///her2: 0///size: 1///node: 1///age: 54///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: 630///e.rfs: 0///t.os: 630///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
        samplename: a44k10829.1282///dataset: UNC4///series: UNC4///id: NA///er: 1///pgr: 1///her2: 0///size: 3///node: 1///age: 54///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 630///e.rfs: 0///t.os: 630///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
        samplename: a44k10829.1283///dataset: UNC4///series: UNC4///id: NA///er: 1///pgr: 1///her2: 0///size: 3///node: 0///age: 44///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 750///e.rfs: 0///t.os: 750///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
        samplename: a44k11240.2741///dataset: UNC4///series: UNC4///id: NA///er: 1///pgr: 1///her2: 0///size: 3///node: 0///age: 53///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 960///e.rfs: 0///t.os: 960///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
        samplename: a44k11240.2760///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: 0///her2: 0///size: 3///node: 0///age: 82///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 930///e.rfs: 0///t.os: 930///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
       samplename: a44k11240.2761///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: 0///her2: 1///size: 6///node: 1///age: 58///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: 900///e.rfs: 0///t.os: 900///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
          samplename: h1a.013na.c///dataset: UNC4///series: UNC4///id: NA///er: 1///pgr: 1///her2: 0///size: NA///node: 1///age: 43///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 660///e.rfs: 0///t.os: 660///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
              samplename: h1a.c.0077///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: 0///her2: 0///size: 6///node: 1///age: 68///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 90///e.rfs: 1///t.os: 90///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
         samplename: h1a.c.0078///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: 0///her2: 0///size: NA///node: 1///age: 44///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1860///e.rfs: 0///t.os: 1860///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
            samplename: h1a.c.0079///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: 0///her2: 1///size: 3///node: 1///age: 55///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 300///e.rfs: 0///t.os: 300///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
          samplename: h1a.c.0080///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: 0///her2: 0///size: 3///node: 0///age: 43///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2280///e.rfs: 0///t.os: 2280///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
          samplename: h1a.c.0081///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: 0///her2: 0///size: 3///node: 1///age: 43///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1440///e.rfs: 1///t.os: 2160///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
        samplename: h1a.v2.1313///dataset: UNC4///series: UNC4///id: NA///er: NA///pgr: NA///her2: NA///size: 1///node: 0///age: 55///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
       samplename: h1a.v2.1314///dataset: UNC4///series: UNC4///id: NA///er: NA///pgr: NA///her2: NA///size: 3///node: 1///age: 65///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
         samplename: h1a.v2.2411na///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: 0///her2: 0///size: 3///node: 1///age: 36///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
       samplename: h1a0120na.c///dataset: UNC4///series: UNC4///id: NA///er: NA///pgr: 0///her2: 0///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
        samplename: h1a012na.c///dataset: UNC4///series: UNC4///id: NA///er: NA///pgr: 0///her2: 0///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
          samplename: h1a014na.c///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: 0///her2: 0///size: 3///node: 0///age: 44///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1950///e.rfs: 0///t.os: 1950///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
         samplename: h1a2371v2///dataset: UNC4///series: UNC4///id: NA///er: NA///pgr: 1///her2: 1///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
          samplename: h1a4110.3579///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: 0///her2: 0///size: 6///node: 1///age: 68///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 360///e.rfs: 0///t.os: 360///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
          samplename: h1a4110.3580///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: 0///her2: 0///size: 3///node: 0///age: 84///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 810///e.rfs: 0///t.os: 810///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
         samplename: h1a4110.3582///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: NA///her2: 0///size: 3///node: 0///age: 51///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
         samplename: h1a4594v2///dataset: UNC4///series: UNC4///id: NA///er: NA///pgr: 0///her2: 0///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
          samplename: h1a4956v2.na///dataset: UNC4///series: UNC4///id: NA///er: 1///pgr: 1///her2: 1///size: 3///node: 1///age: 46///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 30///e.rfs: 1///t.os: 1320///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
      samplename: h1a5365.899///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: 61///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: 30///e.rfs: 1///t.os: 1530///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
         samplename: h1a5365.900///dataset: UNC4///series: UNC4///id: NA///er: 1///pgr: 0///her2: 0///size: 1///node: 1///age: 32///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 1530///e.rfs: 1///t.os: 1920///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
      samplename: h1a5949.019///dataset: UNC4///series: UNC4///id: NA///er: NA///pgr: NA///her2: 0///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
           samplename: h1a6149.098///dataset: UNC4///series: UNC4///id: NA///er: 1///pgr: 0///her2: 0///size: 3///node: 0///age: 64///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 480///e.rfs: 0///t.os: 480///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
      samplename: h1a6149.698///dataset: UNC4///series: UNC4///id: NA///er: NA///pgr: NA///her2: 1///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
          samplename: h1a6260.731///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: 0///her2: 1///size: NA///node: 1///age: 51///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 240///e.rfs: 1///t.os: 360///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
           samplename: h1a6260.732///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: 0///her2: 0///size: 1///node: 0///age: 24///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 570///e.rfs: 0///t.os: 570///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: h1a6261.0922v2///dataset: UNC4///series: UNC4///id: NA///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
        samplename: h1a6261.0924v2///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: 0///her2: 0///size: 3///node: 0///age: 58///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 900///e.rfs: 0///t.os: 900///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
    samplename: h1a6261.0925v2///dataset: UNC4///series: UNC4///id: NA///er: NA///pgr: 0///her2: 0///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
     samplename: h1a6261.1488v2///dataset: UNC4///series: UNC4///id: NA///er: 1///pgr: 0///her2: NA///size: 3///node: 1///age: 52///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1080///e.rfs: 1///t.os: 1920///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
samplename: h1a6261.1489v2.na///dataset: UNC4///series: UNC4///id: NA///er: NA///pgr: NA///her2: 0///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
        samplename: h1a6261.1489v2///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: 0///her2: 0///size: 6///node: 1///age: NA///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 600///e.rfs: 1///t.os: 630///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
    samplename: h1a6261.1490v2///dataset: UNC4///series: UNC4///id: NA///er: NA///pgr: 0///her2: 0///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
        samplename: h1a6261.1499v2///dataset: UNC4///series: UNC4///id: NA///er: 1///pgr: 1///her2: 0///size: 1///node: 1///age: 63///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 270///e.rfs: 0///t.os: 270///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
        samplename: h1a6261.1500v2///dataset: UNC4///series: UNC4///id: NA///er: 1///pgr: 1///her2: 0///size: 1///node: 0///age: 60///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 510///e.rfs: 0///t.os: 510///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
       samplename: h1a6261.2087v2///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: NA///her2: 0///size: 3///node: 1///age: 25///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 270///e.rfs: 0///t.os: 270///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
        samplename: h1a6261.2140v2///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: 0///her2: 0///size: 3///node: 0///age: 49///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 300///e.rfs: 0///t.os: 300///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
    samplename: h1a6261.2141v2///dataset: UNC4///series: UNC4///id: NA///er: NA///pgr: 0///her2: 0///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
    samplename: h1a6261.2597v2.na///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: NA///her2: NA///size: 3///node: 1///age: 61///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 30///e.rfs: 1///t.os: 390///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
samplename: h1a6261.2598v2.na///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: NA///her2: NA///size: NA///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: 930///e.rfs: 1///t.os: 2220///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
    samplename: h1a6261.2607v2///dataset: UNC4///series: UNC4///id: NA///er: NA///pgr: 0///her2: 1///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
 samplename: h1a6261.2980v2.na///dataset: UNC4///series: UNC4///id: NA///er: NA///pgr: 0///her2: 1///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
        samplename: h1a6433.2272v2///dataset: UNC4///series: UNC4///id: NA///er: 1///pgr: 1///her2: 0///size: 1///node: 1///age: 48///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 570///e.rfs: 0///t.os: 570///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
       samplename: h1a6433.2273.v2///dataset: UNC4///series: UNC4///id: NA///er: 1///pgr: 1///her2: 0///size: 3///node: 1///age: 50///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 720///e.rfs: 0///t.os: 720///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
       samplename: h1a6433.2274.v2///dataset: UNC4///series: UNC4///id: NA///er: 1///pgr: 1///her2: 0///size: 3///node: 1///age: 46///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 630///e.rfs: 0///t.os: 630///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
       samplename: h1a6433.2275.v2///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: 0///her2: 1///size: 3///node: 1///age: 72///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 600///e.rfs: 0///t.os: 600///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
     samplename: h1a6433.2276.v2///dataset: UNC4///series: UNC4///id: NA///er: 1///pgr: 0///her2: 0///size: NA///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
    samplename: h1a6433.2597v2///dataset: UNC4///series: UNC4///id: NA///er: NA///pgr: 1///her2: 0///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
    samplename: h1a6433.2598v2///dataset: UNC4///series: UNC4///id: NA///er: NA///pgr: 0///her2: 0///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
       samplename: h1a6433.2599v2///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: 0///her2: 1///size: 6///node: 1///age: NA///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 330///e.rfs: 1///t.os: 1260///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
    samplename: h1a6433.2600v2///dataset: UNC4///series: UNC4///id: NA///er: NA///pgr: 1///her2: 0///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
    samplename: h1a6433.2601v2///dataset: UNC4///series: UNC4///id: NA///er: NA///pgr: 0///her2: 1///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
     samplename: h1a6433.2608v2.na///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: 0///her2: 0///size: 1///node: 0///age: 50///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 570///e.rfs: 0///t.os: 570///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
        samplename: h1a6433.2608v2///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: 0///her2: 1///size: 3///node: 1///age: 71///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 810///e.rfs: 1///t.os: 840///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
    samplename: h1a6433.2610v2.na///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: 0///her2: 0///size: NA///node: 1///age: 39///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 150///e.rfs: 0///t.os: 150///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
    samplename: h1a6433.2610v2///dataset: UNC4///series: UNC4///id: NA///er: NA///pgr: 0///her2: 0///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
      samplename: h1a6433.2611v2///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: 0///her2: 1///size: 6///node: 0///age: 63///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1290///e.rfs: 0///t.os: 1290///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
   samplename: h1a6433.2757.v2///dataset: UNC4///series: UNC4///id: NA///er: NA///pgr: 0///her2: 0///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
         samplename: h1a6433.2758.v2///dataset: UNC4///series: UNC4///id: NA///er: 1///pgr: 0///her2: 0///size: 3///node: 1///age: 46///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 60///e.rfs: 0///t.os: 60///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
   samplename: h1a6433.2759.v2///dataset: UNC4///series: UNC4///id: NA///er: NA///pgr: 0///her2: 0///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
     samplename: h1a6433.2761.v2///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: NA///her2: NA///size: NA///node: 1///age: 29///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 30///e.rfs: 1///t.os: 480///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
    samplename: h1a6433.2808.v2///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: NA///her2: NA///size: 1///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: 750///e.rfs: 0///t.os: 750///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
       samplename: h1a6433.2809.v2///dataset: UNC4///series: UNC4///id: NA///er: 1///pgr: 1///her2: 0///size: 3///node: 0///age: 80///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 480///e.rfs: 0///t.os: 480///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
      samplename: h1a6433.2810.v2///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: NA///her2: 0///size: 3///node: 0///age: 58///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: h1a6433.2811.v2///dataset: UNC4///series: UNC4///id: NA///er: NA///pgr: NA///her2: 0///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
       samplename: h1a6433.2812.v2///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: 0///her2: 1///size: 6///node: 1///age: 52///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 270///e.rfs: 0///t.os: 270///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
      samplename: h1a6433.2978.v2///dataset: UNC4///series: UNC4///id: NA///er: 1///pgr: NA///her2: 0///size: 3///node: 1///age: 52///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
       samplename: h1a6433.2979.v2///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: NA///her2: NA///size: 6///node: 1///age: 84///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 60///e.rfs: 0///t.os: 60///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
    samplename: h1a6433.2987v2///dataset: UNC4///series: UNC4///id: NA///er: NA///pgr: 0///her2: 1///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
      samplename: h1a6433.2988v2///dataset: UNC4///series: UNC4///id: NA///er: 1///pgr: 1///her2: 0///size: 3///node: 1///age: 88///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 1170///e.rfs: 0///t.os: 1170///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
    samplename: h1a6433.2989v2///dataset: UNC4///series: UNC4///id: NA///er: NA///pgr: 0///her2: 1///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
      samplename: h1a6433.2990v2///dataset: UNC4///series: UNC4///id: NA///er: 1///pgr: 1///her2: 0///size: 3///node: 0///age: NA///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1110///e.rfs: 0///t.os: 1110///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
        samplename: h1a6433.2992v2///dataset: UNC4///series: UNC4///id: NA///er: 1///pgr: 1///her2: 0///size: 6///node: 1///age: 40///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 990///e.rfs: 0///t.os: 990///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
       samplename: h1a6433.3007.v2///dataset: UNC4///series: UNC4///id: NA///er: 1///pgr: 1///her2: 0///size: 1///node: 0///age: 77///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 900///e.rfs: 0///t.os: 900///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
     samplename: h1a6433.3009.v2///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: NA///her2: NA///size: 1///node: 0///age: NA///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 900///e.rfs: 0///t.os: 900///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
       samplename: h1a6433.3010.v2///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: 0///her2: 0///size: 3///node: 0///age: 59///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 750///e.rfs: 0///t.os: 750///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
       samplename: h1a6433.3011.v2///dataset: UNC4///series: UNC4///id: NA///er: 1///pgr: 0///her2: 0///size: 1///node: 0///age: 48///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 750///e.rfs: 0///t.os: 750///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
         samplename: h1a6546.126///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: NA///her2: 0///size: 3///node: 1///age: 43///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 990///e.rfs: 1///t.os: 1080///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
       samplename: h1a6730.644///dataset: UNC4///series: UNC4///id: NA///er: NA///pgr: 0///her2: 0///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
          samplename: h1a6730.735///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: NA///her2: 0///size: 6///node: 0///age: 77///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 150///e.rfs: 1///t.os: 330///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
          samplename: h1a6730.744///dataset: UNC4///series: UNC4///id: NA///er: 1///pgr: 0///her2: 0///size: 3///node: 0///age: 48///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 840///e.rfs: 1///t.os: 1530///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
       samplename: h1a6730.745///dataset: UNC4///series: UNC4///id: NA///er: NA///pgr: 0///her2: 0///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
   samplename: h1a6912.6851v2///dataset: UNC4///series: UNC4///id: NA///er: NA///pgr: NA///her2: 0///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
      samplename: h1a6912.6852v2///dataset: UNC4///series: UNC4///id: NA///er: 1///pgr: NA///her2: 0///size: NA///node: 1///age: 44///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 660///e.rfs: 1///t.os: 660///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
       samplename: h1a6912.6853v2///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: 0///her2: 0///size: NA///node: 1///age: 65///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 780///e.rfs: 1///t.os: 780///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
    samplename: h1a6912.6955v2///dataset: UNC4///series: UNC4///id: NA///er: NA///pgr: 0///her2: 0///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
     samplename: h1a6912.6956v2///dataset: UNC4///series: UNC4///id: NA///er: 1///pgr: NA///her2: NA///size: 1///node: 1///age: 42///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: 660///e.rfs: 0///t.os: 660///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
       samplename: h1a6912.6957v2///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: 1///her2: 0///size: 1///node: 1///age: 48///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 690///e.rfs: 1///t.os: 1410///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
      samplename: h1a6912.6958v2///dataset: UNC4///series: UNC4///id: NA///er: 1///pgr: 1///her2: NA///size: 3///node: 0///age: 82///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 930///e.rfs: 1///t.os: 2820///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
      samplename: h1a6912.6984v2///dataset: UNC4///series: UNC4///id: NA///er: 0///pgr: 0///her2: 0///size: 3///node: 1///age: 74///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 3180///e.rfs: 0///t.os: 3180///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
    samplename: h1a6912.6985v2///dataset: UNC4///series: UNC4///id: NA///er: 1///pgr: NA///her2: NA///size: NA///node: 0///age: 72///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 780///e.rfs: 1///t.os: 3030///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
     samplename: h1a6912.6986v2///dataset: UNC4///series: UNC4///id: NA///er: 1///pgr: 1///her2: NA///size: 3///node: 0///age: 74///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2880///e.rfs: 0///t.os: 2880///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
    samplename: h1a6912.6987v2///dataset: UNC4///series: UNC4///id: NA///er: 1///pgr: NA///her2: 0///size: NA///node: 1///age: 47///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2880///e.rfs: 0///t.os: 2880///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
                                                                                                                                                                                                                                           (Other) 
                                                                                                                                                                                                                                               206 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse3494 }
\keyword{datasets}
