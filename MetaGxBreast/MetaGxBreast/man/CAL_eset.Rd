\name{ CAL_eset }
\alias{ CAL_eset }
\docType{data}
\title{  }
\usage{data( CAL_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/17157792 
  Title:  
  URL: http://www.ebi.ac.uk/arrayexpress/experiments/E-TABM-158/ 
  PMIDs: 17157792 
  No abstract available.
  notes:
   summary:     
      Recurrent copy number abnormalities differ between tumor subtypes as defined by gene expression patterns.  Accuracy of st
ratification by outcome can be improved by combining expression and copy number.
   version:     
      2015-04-27 18:40:40
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: 1007_s_at 1053_at ... AFFX-HUMISGF3A/M97935_MB_at (21169 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 21169 features, 118 samples
Platform type: 
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

   1 observation deleted due to missingness 
records   n.max n.start  events  median 0.95LCL 0.95UCL 
 117.00  117.00  117.00   77.00    8.96    8.33    9.71 

--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
CAL_091604_WN405650_LBLPOPII.A02 CAL_091604_WN405650_LBLPOPII.A03 CAL_091604_WN405650_LBLPOPII.A04 
                               1                                1                                1 
CAL_091604_WN405650_LBLPOPII.A11 CAL_091604_WN405650_LBLPOPII.B02 CAL_091604_WN405650_LBLPOPII.B04 
                               1                                1                                1 
CAL_091604_WN405650_LBLPOPII.B10 CAL_091604_WN405650_LBLPOPII.B11 CAL_091604_WN405650_LBLPOPII.C03 
                               1                                1                                1 
CAL_091604_WN405650_LBLPOPII.C04 CAL_091604_WN405650_LBLPOPII.C11 CAL_091604_WN405650_LBLPOPII.C12 
                               1                                1                                1 
CAL_091604_WN405650_LBLPOPII.D04 CAL_091604_WN405650_LBLPOPII.D10 CAL_091604_WN405650_LBLPOPII.D11 
                               1                                1                                1 
CAL_091604_WN405650_LBLPOPII.E01 CAL_091604_WN405650_LBLPOPII.E10 CAL_091604_WN405650_LBLPOPII.E11 
                               1                                1                                1 
CAL_091604_WN405650_LBLPOPII.F01 CAL_091604_WN405650_LBLPOPII.F02 CAL_091604_WN405650_LBLPOPII.F09 
                               1                                1                                1 
CAL_091604_WN405650_LBLPOPII.F10 CAL_091604_WN405650_LBLPOPII.G01 CAL_091604_WN405650_LBLPOPII.G02 
                               1                                1                                1 
CAL_091604_WN405650_LBLPOPII.G03 CAL_091604_WN405650_LBLPOPII.G09 CAL_091604_WN405650_LBLPOPII.G10 
                               1                                1                                1 
CAL_091604_WN405650_LBLPOPII.G11 CAL_091604_WN405650_LBLPOPII.H02 CAL_091604_WN405650_LBLPOPII.H03 
                               1                                1                                1 
         CAL_LBL_POP_W405609.A02          CAL_LBL_POP_W405609.A03          CAL_LBL_POP_W405609.A04 
                               1                                1                                1 
         CAL_LBL_POP_W405609.A05          CAL_LBL_POP_W405609.A06          CAL_LBL_POP_W405609.A07 
                               1                                1                                1 
         CAL_LBL_POP_W405609.A10          CAL_LBL_POP_W405609.A12          CAL_LBL_POP_W405609.B02 
                               1                                1                                1 
         CAL_LBL_POP_W405609.B03          CAL_LBL_POP_W405609.B04          CAL_LBL_POP_W405609.B06 
                               1                                1                                1 
         CAL_LBL_POP_W405609.B08          CAL_LBL_POP_W405609.B09          CAL_LBL_POP_W405609.B10 
                               1                                1                                1 
         CAL_LBL_POP_W405609.B11          CAL_LBL_POP_W405609.B12          CAL_LBL_POP_W405609.C02 
                               1                                1                                1 
         CAL_LBL_POP_W405609.C04          CAL_LBL_POP_W405609.C08          CAL_LBL_POP_W405609.C09 
                               1                                1                                1 
         CAL_LBL_POP_W405609.C10          CAL_LBL_POP_W405609.C11          CAL_LBL_POP_W405609.C12 
                               1                                1                                1 
         CAL_LBL_POP_W405609.D01          CAL_LBL_POP_W405609.D04          CAL_LBL_POP_W405609.D05 
                               1                                1                                1 
         CAL_LBL_POP_W405609.D10          CAL_LBL_POP_W405609.D11          CAL_LBL_POP_W405609.E01 
                               1                                1                                1 
         CAL_LBL_POP_W405609.E02          CAL_LBL_POP_W405609.E03          CAL_LBL_POP_W405609.E04 
                               1                                1                                1 
         CAL_LBL_POP_W405609.E05          CAL_LBL_POP_W405609.E06          CAL_LBL_POP_W405609.E07 
                               1                                1                                1 
         CAL_LBL_POP_W405609.E08          CAL_LBL_POP_W405609.E09          CAL_LBL_POP_W405609.E10 
                               1                                1                                1 
         CAL_LBL_POP_W405609.E11          CAL_LBL_POP_W405609.F01          CAL_LBL_POP_W405609.F02 
                               1                                1                                1 
         CAL_LBL_POP_W405609.F03          CAL_LBL_POP_W405609.F04          CAL_LBL_POP_W405609.F05 
                               1                                1                                1 
         CAL_LBL_POP_W405609.F06          CAL_LBL_POP_W405609.F07          CAL_LBL_POP_W405609.F08 
                               1                                1                                1 
         CAL_LBL_POP_W405609.F09          CAL_LBL_POP_W405609.F10          CAL_LBL_POP_W405609.F11 
                               1                                1                                1 
         CAL_LBL_POP_W405609.G01          CAL_LBL_POP_W405609.G02          CAL_LBL_POP_W405609.G03 
                               1                                1                                1 
         CAL_LBL_POP_W405609.G04          CAL_LBL_POP_W405609.G05          CAL_LBL_POP_W405609.G06 
                               1                                1                                1 
         CAL_LBL_POP_W405609.G07          CAL_LBL_POP_W405609.G08          CAL_LBL_POP_W405609.G09 
                               1                                1                                1 
         CAL_LBL_POP_W405609.G10          CAL_LBL_POP_W405609.G11          CAL_LBL_POP_W405609.H01 
                               1                                1                                1 
         CAL_LBL_POP_W405609.H02          CAL_LBL_POP_W405609.H03          CAL_LBL_POP_W405609.H04 
                               1                                1                                1 
         CAL_LBL_POP_W405609.H05          CAL_LBL_POP_W405609.H06          CAL_LBL_POP_W405609.H07 
                               1                                1                                1 
                         (Other) 
                              19 

sample_type: 
tumor 
  118 

er: 
negative positive 
      43       75 

pgr: 
negative positive     NA's 
      51       66        1 

tumor_size: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  0.300   1.675   2.300   2.729   3.500   7.500       2 

N: 
 0  1 
51 67 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  31.00   44.00   51.00   55.06   66.00   88.00       1 

grade: 
   1    2    3 NA's 
  10   42   61    5 

dmfs_days: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
      0     767    2059    2094    3336    5183       1 

dmfs_status: 
deceased_or_recurrence    living_norecurrence                   NA's 
                    26                     91                      1 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
      0     767    2059    2094    3336    5183       1 

recurrence_status: 
deceased_or_recurrence    living_norecurrence                   NA's 
                    36                     81                      1 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
     47    1117    2234    2347    3504    5183       1 

vital_status: 
deceased   living 
      77       41 

treatment: 
chemo.plus.hormono       chemotherapy     hormonotherapy          untreated               NA's 
                25                 36                 40                 14                  3 

batch: 
CAL 
118 

uncurated_author_metadata: 
  samplename: CAL_091604_WN405650_LBLPOPII.A02///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 0///her2: NA///size: 1.6///node: 1///age: 39///grade: 3///t.dmfs: 113///e.dmfs: 1///t.rfs: 113///e.rfs: 1///t.os: 1460///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                              1 
   samplename: CAL_091604_WN405650_LBLPOPII.A03///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 1.4///node: 0///age: 51///grade: 3///t.dmfs: 153///e.dmfs: 1///t.rfs: 153///e.rfs: 1///t.os: 861///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                              1 
samplename: CAL_091604_WN405650_LBLPOPII.A04///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 1.5///node: 1///age: 43///grade: 3///t.dmfs: 2055///e.dmfs: 0///t.rfs: 2055///e.rfs: 0///t.os: 2055///e.os: 1///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                              1 
samplename: CAL_091604_WN405650_LBLPOPII.A11///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 0///her2: NA///size: 1.8///node: 1///age: 60///grade: 3///t.dmfs: 2727///e.dmfs: 0///t.rfs: 2727///e.rfs: 0///t.os: 2727///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                              1 
samplename: CAL_091604_WN405650_LBLPOPII.B02///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 1///her2: NA///size: 2.1///node: 0///age: 44///grade: 3///t.dmfs: 1632///e.dmfs: 0///t.rfs: 1632///e.rfs: 0///t.os: 1632///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                              1 
samplename: CAL_091604_WN405650_LBLPOPII.B04///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 0///her2: NA///size: 2.2///node: 0///age: 45///grade: 2///t.dmfs: 1212///e.dmfs: 0///t.rfs: 1212///e.rfs: 1///t.os: 4143///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                              1 
samplename: CAL_091604_WN405650_LBLPOPII.B10///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 0///her2: NA///size: 1.8///node: 1///age: 36///grade: 2///t.dmfs: 5165///e.dmfs: 0///t.rfs: 5165///e.rfs: 0///t.os: 5165///e.os: 1///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                              1 
  samplename: CAL_091604_WN405650_LBLPOPII.B11///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 1///node: 0///age: 42///grade: 3///t.dmfs: 3113///e.dmfs: 0///t.rfs: 3113///e.rfs: 0///t.os: 3113///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                              1 
samplename: CAL_091604_WN405650_LBLPOPII.C03///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 0///her2: NA///size: 3.5///node: 1///age: 50///grade: 3///t.dmfs: 3139///e.dmfs: 0///t.rfs: 3139///e.rfs: 0///t.os: 3139///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                              1 
   samplename: CAL_091604_WN405650_LBLPOPII.C04///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 1.5///node: 1///age: 59///grade: 3///t.dmfs: 412///e.dmfs: 0///t.rfs: 412///e.rfs: 0///t.os: 412///e.os: 1///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                              1 
samplename: CAL_091604_WN405650_LBLPOPII.C11///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 5.5///node: 1///age: 72///grade: 3///t.dmfs: 3913///e.dmfs: 0///t.rfs: 3913///e.rfs: 0///t.os: 3913///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
  samplename: CAL_091604_WN405650_LBLPOPII.C12///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 7///node: 1///age: 57///grade: 3///t.dmfs: 2321///e.dmfs: 0///t.rfs: 2321///e.rfs: 1///t.os: 2387///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
samplename: CAL_091604_WN405650_LBLPOPII.D04///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 0///her2: NA///size: 1.6///node: 0///age: 66///grade: 1///t.dmfs: 3376///e.dmfs: 0///t.rfs: 3376///e.rfs: 0///t.os: 3376///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
samplename: CAL_091604_WN405650_LBLPOPII.D10///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 0///her2: NA///size: 2.5///node: 1///age: 46///grade: 3///t.dmfs: 2048///e.dmfs: 0///t.rfs: 2048///e.rfs: 0///t.os: 2048///e.os: 1///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                              1 
samplename: CAL_091604_WN405650_LBLPOPII.D11///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 1.3///node: 0///age: 80///grade: 3///t.dmfs: 1117///e.dmfs: 0///t.rfs: 1117///e.rfs: 0///t.os: 1117///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
  samplename: CAL_091604_WN405650_LBLPOPII.E01///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 2///node: 0///age: 83///grade: 2///t.dmfs: 2478///e.dmfs: 1///t.rfs: 2478///e.rfs: 1///t.os: 3362///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
samplename: CAL_091604_WN405650_LBLPOPII.E10///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 0///her2: NA///size: 1.9///node: 0///age: 36///grade: 2///t.dmfs: 3723///e.dmfs: 0///t.rfs: 3723///e.rfs: 0///t.os: 3723///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                              1 
samplename: CAL_091604_WN405650_LBLPOPII.E11///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 0///her2: NA///size: 1.4///node: 1///age: 74///grade: 3///t.dmfs: 2986///e.dmfs: 0///t.rfs: 2986///e.rfs: 0///t.os: 2986///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                              1 
   samplename: CAL_091604_WN405650_LBLPOPII.F01///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 0///her2: NA///size: 2.5///node: 0///age: 43///grade: 2///t.dmfs: 777///e.dmfs: 0///t.rfs: 777///e.rfs: 1///t.os: 902///e.os: 1///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                              1 
samplename: CAL_091604_WN405650_LBLPOPII.F02///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 3.6///node: 1///age: 67///grade: 3///t.dmfs: 2949///e.dmfs: 0///t.rfs: 2949///e.rfs: 0///t.os: 2949///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
samplename: CAL_091604_WN405650_LBLPOPII.F09///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 1.4///node: 1///age: 63///grade: 2///t.dmfs: 4698///e.dmfs: 0///t.rfs: 4698///e.rfs: 0///t.os: 4698///e.os: 1///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                              1 
    samplename: CAL_091604_WN405650_LBLPOPII.F10///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 2///node: 0///age: 48///grade: 2///t.dmfs: 989///e.dmfs: 1///t.rfs: 989///e.rfs: 1///t.os: 1380///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
  samplename: CAL_091604_WN405650_LBLPOPII.G01///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 3///node: 0///age: 47///grade: 3///t.dmfs: 4103///e.dmfs: 0///t.rfs: 4103///e.rfs: 0///t.os: 4103///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
samplename: CAL_091604_WN405650_LBLPOPII.G02///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 2.4///node: 0///age: 51///grade: 3///t.dmfs: 3745///e.dmfs: 0///t.rfs: 3745///e.rfs: 0///t.os: 3745///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
     samplename: CAL_091604_WN405650_LBLPOPII.G03///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 0///her2: NA///size: 5///node: 1///age: 82///grade: 2///t.dmfs: 785///e.dmfs: 1///t.rfs: 785///e.rfs: 1///t.os: 920///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
  samplename: CAL_091604_WN405650_LBLPOPII.G09///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 0///her2: NA///size: 6///node: 1///age: 50///grade: 1///t.dmfs: 3628///e.dmfs: 0///t.rfs: 3628///e.rfs: 0///t.os: 3628///e.os: 1///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                              1 
samplename: CAL_091604_WN405650_LBLPOPII.G10///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 2.1///node: 1///age: 71///grade: 2///t.dmfs: 1518///e.dmfs: 0///t.rfs: 1518///e.rfs: 0///t.os: 1518///e.os: 1///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                              1 
   samplename: CAL_091604_WN405650_LBLPOPII.G11///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 3.5///node: 1///age: 48///grade: 3///t.dmfs: 496///e.dmfs: 0///t.rfs: 496///e.rfs: 0///t.os: 496///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                              1 
  samplename: CAL_091604_WN405650_LBLPOPII.H02///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 0///her2: NA///size: 2.3///node: 0///age: 63///grade: NA///t.dmfs: 748///e.dmfs: 0///t.rfs: 748///e.rfs: 0///t.os: 748///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                              1 
samplename: CAL_091604_WN405650_LBLPOPII.H03///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 1.6///node: 0///age: 61///grade: 3///t.dmfs: 1927///e.dmfs: 0///t.rfs: 1927///e.rfs: 0///t.os: 1927///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
         samplename: CAL_LBL_POP_W405609.A02///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 0///her2: NA///size: 1.5///node: 0///age: 43///grade: 1///t.dmfs: 3270///e.dmfs: 0///t.rfs: 3270///e.rfs: 0///t.os: 3270///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                              1 
         samplename: CAL_LBL_POP_W405609.A03///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 1.5///node: 0///age: 68///grade: 2///t.dmfs: 3781///e.dmfs: 0///t.rfs: 3781///e.rfs: 0///t.os: 3781///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
           samplename: CAL_LBL_POP_W405609.A04///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 2///node: 0///age: 40///grade: 2///t.dmfs: 1475///e.dmfs: 1///t.rfs: 1475///e.rfs: 1///t.os: 2102///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                              1 
         samplename: CAL_LBL_POP_W405609.A05///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 0///her2: NA///size: 3.4///node: 0///age: 49///grade: 2///t.dmfs: 3289///e.dmfs: 0///t.rfs: 3289///e.rfs: 0///t.os: 3289///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                              1 
           samplename: CAL_LBL_POP_W405609.A06///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 0///her2: NA///size: 3///node: 0///age: 78///grade: 2///t.dmfs: 2081///e.dmfs: 0///t.rfs: 2081///e.rfs: 0///t.os: 2081///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
         samplename: CAL_LBL_POP_W405609.A07///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 2.5///node: 1///age: 83///grade: 2///t.dmfs: 3577///e.dmfs: 0///t.rfs: 3577///e.rfs: 0///t.os: 3577///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
           samplename: CAL_LBL_POP_W405609.A10///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 0///her2: NA///size: 2///node: 1///age: 72///grade: 2///t.dmfs: 3165///e.dmfs: 0///t.rfs: 3165///e.rfs: 0///t.os: 4219///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                              1 
         samplename: CAL_LBL_POP_W405609.A12///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 1///her2: NA///size: 3.5///node: 1///age: 40///grade: 2///t.dmfs: 3504///e.dmfs: 0///t.rfs: 3504///e.rfs: 0///t.os: 3504///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
         samplename: CAL_LBL_POP_W405609.B02///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 0///her2: NA///size: 2.5///node: 1///age: 46///grade: 3///t.dmfs: 2234///e.dmfs: 0///t.rfs: 2234///e.rfs: 0///t.os: 2234///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                              1 
           samplename: CAL_LBL_POP_W405609.B03///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 0///her2: NA///size: 2///node: 0///age: 57///grade: 1///t.dmfs: 4008///e.dmfs: 0///t.rfs: 4008///e.rfs: 0///t.os: 4008///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                              1 
         samplename: CAL_LBL_POP_W405609.B04///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 1///her2: NA///size: 5.5///node: 0///age: 54///grade: 3///t.dmfs: 1000///e.dmfs: 0///t.rfs: 1000///e.rfs: 0///t.os: 1000///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
              samplename: CAL_LBL_POP_W405609.B06///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 4///node: 1///age: 73///grade: 2///t.dmfs: 767///e.dmfs: 0///t.rfs: 767///e.rfs: 0///t.os: 767///e.os: 1///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                              1 
         samplename: CAL_LBL_POP_W405609.B08///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 0///her2: NA///size: 1.5///node: 0///age: 55///grade: 2///t.dmfs: 2467///e.dmfs: 0///t.rfs: 2467///e.rfs: 0///t.os: 2467///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
              samplename: CAL_LBL_POP_W405609.B09///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 3///node: 1///age: 80///grade: 3///t.dmfs: 562///e.dmfs: 0///t.rfs: 562///e.rfs: 0///t.os: 489///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
             samplename: CAL_LBL_POP_W405609.B10///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: NA///node: 1///age: 39///grade: 3///t.dmfs: 672///e.dmfs: 0///t.rfs: 672///e.rfs: 0///t.os: 672///e.os: 0///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                              1 
           samplename: CAL_LBL_POP_W405609.B11///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 2.5///node: 0///age: 37///grade: 2///t.dmfs: 978///e.dmfs: 1///t.rfs: 978///e.rfs: 1///t.os: 1248///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                              1 
            samplename: CAL_LBL_POP_W405609.B12///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 0///her2: NA///size: 1.9///node: 1///age: 51///grade: 2///t.dmfs: 551///e.dmfs: 0///t.rfs: 551///e.rfs: 0///t.os: 551///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
             samplename: CAL_LBL_POP_W405609.C02///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 0///her2: NA///size: 5///node: 1///age: 51///grade: 2///t.dmfs: 573///e.dmfs: 0///t.rfs: 573///e.rfs: 1///t.os: 1288///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                              1 
        samplename: CAL_LBL_POP_W405609.C04///dataset: CAL///series: CAL///id: NA///er: 0///pgr: NA///her2: NA///size: 1.5///node: 0///age: 50///grade: 2///t.dmfs: 1321///e.dmfs: 1///t.rfs: 1321///e.rfs: 1///t.os: 2070///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                              1 
           samplename: CAL_LBL_POP_W405609.C08///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 3///node: 1///age: 88///grade: 1///t.dmfs: 2924///e.dmfs: 0///t.rfs: 2924///e.rfs: 0///t.os: 2924///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
            samplename: CAL_LBL_POP_W405609.C09///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 0///her2: NA///size: 1.5///node: 0///age: 45///grade: 1///t.dmfs: 314///e.dmfs: 0///t.rfs: 314///e.rfs: 1///t.os: 314///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                              1 
           samplename: CAL_LBL_POP_W405609.C10///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 1///her2: NA///size: 1.8///node: 1///age: 36///grade: 3///t.dmfs: 558///e.dmfs: 1///t.rfs: 558///e.rfs: 1///t.os: 2084///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                              1 
         samplename: CAL_LBL_POP_W405609.C11///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 1///her2: NA///size: 2.1///node: 0///age: 47///grade: 3///t.dmfs: 3059///e.dmfs: 0///t.rfs: 3059///e.rfs: 0///t.os: 3059///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                              1 
         samplename: CAL_LBL_POP_W405609.C12///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 0///her2: NA///size: 1.5///node: 1///age: 39///grade: 1///t.dmfs: 2059///e.dmfs: 0///t.rfs: 2059///e.rfs: 0///t.os: 2059///e.os: 1///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                              1 
           samplename: CAL_LBL_POP_W405609.D01///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 0///her2: NA///size: 1///node: 0///age: 76///grade: 3///t.dmfs: 3336///e.dmfs: 0///t.rfs: 3336///e.rfs: 0///t.os: 3336///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
           samplename: CAL_LBL_POP_W405609.D04///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 5///node: 1///age: 75///grade: 2///t.dmfs: 3760///e.dmfs: 0///t.rfs: 3760///e.rfs: 0///t.os: 3760///e.os: 1///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                              1 
         samplename: CAL_LBL_POP_W405609.D05///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 2.2///node: 1///age: 41///grade: 2///t.dmfs: 2497///e.dmfs: 0///t.rfs: 2497///e.rfs: 0///t.os: 2497///e.os: 1///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                              1 
        samplename: CAL_LBL_POP_W405609.D10///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 0///her2: NA///size: 3.8///node: 0///age: 44///grade: NA///t.dmfs: 3643///e.dmfs: 0///t.rfs: 3643///e.rfs: 0///t.os: 3643///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                              1 
              samplename: CAL_LBL_POP_W405609.D11///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 2///node: 1///age: 71///grade: 3///t.dmfs: 420///e.dmfs: 1///t.rfs: 420///e.rfs: 1///t.os: 818///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
         samplename: CAL_LBL_POP_W405609.E01///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 0///her2: NA///size: 1.5///node: 0///age: 55///grade: 2///t.dmfs: 3322///e.dmfs: 0///t.rfs: 3322///e.rfs: 0///t.os: 3322///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
           samplename: CAL_LBL_POP_W405609.E02///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 0///her2: NA///size: 3///node: 1///age: 49///grade: 3///t.dmfs: 4132///e.dmfs: 0///t.rfs: 4132///e.rfs: 0///t.os: 4132///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                              1 
           samplename: CAL_LBL_POP_W405609.E03///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 3///node: 1///age: 47///grade: 3///t.dmfs: 2759///e.dmfs: 0///t.rfs: 2759///e.rfs: 0///t.os: 2759///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                              1 
         samplename: CAL_LBL_POP_W405609.E04///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 0///her2: NA///size: 2.2///node: 0///age: 56///grade: 2///t.dmfs: 1580///e.dmfs: 1///t.rfs: 1580///e.rfs: 1///t.os: 1580///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
         samplename: CAL_LBL_POP_W405609.E05///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 1.5///node: 0///age: 56///grade: 2///t.dmfs: 4665///e.dmfs: 0///t.rfs: 4665///e.rfs: 0///t.os: 4665///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
           samplename: CAL_LBL_POP_W405609.E06///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 2.4///node: 1///age: 56///grade: 3///t.dmfs: 124///e.dmfs: 1///t.rfs: 124///e.rfs: 1///t.os: 1069///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
               samplename: CAL_LBL_POP_W405609.E07///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 0///her2: NA///size: 7.5///node: 1///age: 33///grade: 3///t.dmfs: 0///e.dmfs: 1///t.rfs: 0///e.rfs: 1///t.os: 1194///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
          samplename: CAL_LBL_POP_W405609.E08///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 0///her2: NA///size: 5///node: 1///age: 37///grade: NA///t.dmfs: 2796///e.dmfs: 0///t.rfs: 2796///e.rfs: 0///t.os: 2796///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                              1 
         samplename: CAL_LBL_POP_W405609.E09///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 1///her2: NA///size: 3.5///node: 1///age: 43///grade: 2///t.dmfs: 1602///e.dmfs: 0///t.rfs: 1602///e.rfs: 0///t.os: 1602///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                              1 
                 samplename: CAL_LBL_POP_W405609.E10///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 0///her2: NA///size: 2///node: 1///age: 37///grade: 2///t.dmfs: 0///e.dmfs: 1///t.rfs: 0///e.rfs: 1///t.os: 1708///e.os: 0///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                              1 
              samplename: CAL_LBL_POP_W405609.E11///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 3///node: 1///age: 86///grade: 3///t.dmfs: 383///e.dmfs: 0///t.rfs: 383///e.rfs: 0///t.os: 383///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
           samplename: CAL_LBL_POP_W405609.F01///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 4.5///node: 1///age: 55///grade: 3///t.dmfs: 770///e.dmfs: 0///t.rfs: 770///e.rfs: 1///t.os: 1091///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                              1 
           samplename: CAL_LBL_POP_W405609.F02///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 0///her2: NA///size: 1///node: 0///age: 69///grade: 2///t.dmfs: 2533///e.dmfs: 0///t.rfs: 2533///e.rfs: 0///t.os: 2533///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
           samplename: CAL_LBL_POP_W405609.F03///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 0///her2: NA///size: 5///node: 0///age: 86///grade: 1///t.dmfs: 3457///e.dmfs: 0///t.rfs: 3457///e.rfs: 0///t.os: 3457///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                              1 
           samplename: CAL_LBL_POP_W405609.F04///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 1///her2: NA///size: 3///node: 1///age: 66///grade: 2///t.dmfs: 3814///e.dmfs: 0///t.rfs: 3814///e.rfs: 0///t.os: 3814///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
         samplename: CAL_LBL_POP_W405609.F05///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 2.4///node: 1///age: 88///grade: 3///t.dmfs: 2183///e.dmfs: 0///t.rfs: 2183///e.rfs: 0///t.os: 2183///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
         samplename: CAL_LBL_POP_W405609.F06///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 1.1///node: 0///age: 47///grade: 2///t.dmfs: 4008///e.dmfs: 0///t.rfs: 4008///e.rfs: 0///t.os: 4008///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                              1 
         samplename: CAL_LBL_POP_W405609.F07///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 2.5///node: 0///age: 44///grade: 3///t.dmfs: 5019///e.dmfs: 0///t.rfs: 5019///e.rfs: 0///t.os: 5019///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                              1 
        samplename: CAL_LBL_POP_W405609.F08///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 1.5///node: 0///age: 42///grade: NA///t.dmfs: 5183///e.dmfs: 0///t.rfs: 5183///e.rfs: 0///t.os: 5183///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                              1 
            samplename: CAL_LBL_POP_W405609.F09///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 2.3///node: 0///age: 57///grade: 2///t.dmfs: 157///e.dmfs: 0///t.rfs: 157///e.rfs: 0///t.os: 157///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                              1 
                 samplename: CAL_LBL_POP_W405609.F10///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 0///her2: NA///size: 3///node: 1///age: 60///grade: 3///t.dmfs: 47///e.dmfs: 0///t.rfs: 47///e.rfs: 0///t.os: 47///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                              1 
         samplename: CAL_LBL_POP_W405609.F11///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 2.3///node: 1///age: 53///grade: 1///t.dmfs: 3774///e.dmfs: 0///t.rfs: 3774///e.rfs: 0///t.os: 3774///e.os: 1///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                              1 
           samplename: CAL_LBL_POP_W405609.G01///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 4.5///node: 0///age: 50///grade: 2///t.dmfs: 150///e.dmfs: 0///t.rfs: 150///e.rfs: 0///t.os: 150///e.os: 1///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                              1 
           samplename: CAL_LBL_POP_W405609.G02///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 0///her2: NA///size: 6.2///node: 1///age: 52///grade: 3///t.dmfs: 829///e.dmfs: 1///t.rfs: 829///e.rfs: 1///t.os: 1405///e.os: 0///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                              1 
           samplename: CAL_LBL_POP_W405609.G03///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 2.7///node: 0///age: 48///grade: 1///t.dmfs: 902///e.dmfs: 1///t.rfs: 902///e.rfs: 1///t.os: 3219///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
            samplename: CAL_LBL_POP_W405609.G04///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 0///her2: NA///size: 2.4///node: 1///age: 54///grade: 3///t.dmfs: 482///e.dmfs: 1///t.rfs: 482///e.rfs: 1///t.os: 818///e.os: 0///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                              1 
           samplename: CAL_LBL_POP_W405609.G05///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 0///her2: NA///size: 4///node: 0///age: 34///grade: 3///t.dmfs: 1767///e.dmfs: 0///t.rfs: 1767///e.rfs: 0///t.os: 2310///e.os: 1///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                              1 
           samplename: CAL_LBL_POP_W405609.G06///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 0///her2: NA///size: 3///node: 1///age: 40///grade: 2///t.dmfs: 3544///e.dmfs: 0///t.rfs: 3544///e.rfs: 0///t.os: 3544///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                              1 
         samplename: CAL_LBL_POP_W405609.G07///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 1///her2: NA///size: 1.6///node: 1///age: 51///grade: 2///t.dmfs: 1719///e.dmfs: 0///t.rfs: 1719///e.rfs: 0///t.os: 1719///e.os: 1///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                              1 
         samplename: CAL_LBL_POP_W405609.G08///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 1.8///node: 1///age: 35///grade: 3///t.dmfs: 4365///e.dmfs: 0///t.rfs: 4365///e.rfs: 0///t.os: 4365///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
               samplename: CAL_LBL_POP_W405609.G09///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 1.2///node: 0///age: NA///grade: 3///t.dmfs: 0///e.dmfs: 1///t.rfs: 0///e.rfs: 1///t.os: 219///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                              1 
            samplename: CAL_LBL_POP_W405609.G10///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 1.8///node: 1///age: 83///grade: 3///t.dmfs: 245///e.dmfs: 1///t.rfs: 245///e.rfs: 1///t.os: 584///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                              1 
           samplename: CAL_LBL_POP_W405609.G11///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 0.3///node: 0///age: 71///grade: 3///t.dmfs: 299///e.dmfs: 0///t.rfs: 299///e.rfs: 1///t.os: 2588///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
            samplename: CAL_LBL_POP_W405609.H01///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 2.3///node: 0///age: 70///grade: 3///t.dmfs: 905///e.dmfs: 0///t.rfs: 905///e.rfs: 1///t.os: 905///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
                  samplename: CAL_LBL_POP_W405609.H02///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 2///node: 0///age: 49///grade: 2///t.dmfs: 0///e.dmfs: 1///t.rfs: 0///e.rfs: 1///t.os: 887///e.os: 0///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                              1 
         samplename: CAL_LBL_POP_W405609.H03///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 2.2///node: 1///age: 53///grade: 3///t.dmfs: 1051///e.dmfs: 1///t.rfs: 1051///e.rfs: 1///t.os: 1223///e.os: 0///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                              1 
        samplename: CAL_LBL_POP_W405609.H04///dataset: CAL///series: CAL///id: NA///er: 0///pgr: 0///her2: NA///size: 1.8///node: 1///age: 47///grade: 3///t.dmfs: 4117///e.dmfs: 0///t.rfs: 4117///e.rfs: 0///t.os: 4117///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                              1 
           samplename: CAL_LBL_POP_W405609.H05///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 3///node: 1///age: 44///grade: 2///t.dmfs: 2197///e.dmfs: 0///t.rfs: 2197///e.rfs: 0///t.os: 2197///e.os: 1///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                              1 
         samplename: CAL_LBL_POP_W405609.H06///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 0.6///node: 0///age: 82///grade: 3///t.dmfs: 3683///e.dmfs: 0///t.rfs: 3683///e.rfs: 0///t.os: 3683///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
         samplename: CAL_LBL_POP_W405609.H07///dataset: CAL///series: CAL///id: NA///er: 1///pgr: 1///her2: NA///size: 4.5///node: 1///age: 86///grade: 3///t.dmfs: 2456///e.dmfs: 0///t.rfs: 2456///e.rfs: 0///t.os: 2456///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                              1 
                                                                                                                                                                                                                                                        (Other) 
                                                                                                                                                                                                                                                             19 

}}
\source{ http://www.ebi.ac.uk/arrayexpress/experiments/E-TABM-158/ }
\keyword{datasets}
