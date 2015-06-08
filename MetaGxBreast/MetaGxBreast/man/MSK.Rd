\name{ MSK_eset }
\alias{ MSK_eset }
\docType{data}
\title{ MSK_eset }
\usage{data( MSK_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=18592372 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse10510 
  PMIDs: 18592372 
  No abstract available.
  notes:
   summary:     
      A method was developed to separate tumor cells and their microenvironment to test the prognostic abilities of the immune 
system. Results showed that lymphatic infiltration is beneficial for ER negative patients, but probably not beneficial fo
r ER positive patients.
   version:     
      2015-04-27 19:20:52
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
assayData: 20967 features, 99 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
  MSK_B1_T  MSK_B10_T MSK_B100_T MSK_B101_T MSK_B102_T MSK_B103_T MSK_B104_T MSK_B105_T MSK_B106_T MSK_B107_T MSK_B108_T 
         1          1          1          1          1          1          1          1          1          1          1 
MSK_B109_T  MSK_B11_T MSK_B110_T MSK_B111_T MSK_B112_T MSK_B113_T MSK_B114_T MSK_B115_T  MSK_B12_T  MSK_B13_T  MSK_B14_T 
         1          1          1          1          1          1          1          1          1          1          1 
 MSK_B15_T  MSK_B16_T  MSK_B17_T  MSK_B18_T   MSK_B2_T  MSK_B23_T  MSK_B24_T  MSK_B26_T  MSK_B27_T  MSK_B28_T  MSK_B29_T 
         1          1          1          1          1          1          1          1          1          1          1 
  MSK_B3_T  MSK_B31_T  MSK_B32_T  MSK_B33_T   MSK_B4_T  MSK_B42_T  MSK_B43_T  MSK_B44_T  MSK_B45_T  MSK_B46_T  MSK_B47_T 
         1          1          1          1          1          1          1          1          1          1          1 
 MSK_B48_T  MSK_B49_T   MSK_B5_T  MSK_B50_T  MSK_B51_T  MSK_B52_T  MSK_B53_T  MSK_B54_T  MSK_B55_T  MSK_B56_T  MSK_B57_T 
         1          1          1          1          1          1          1          1          1          1          1 
 MSK_B58_T  MSK_B59_T   MSK_B6_T  MSK_B60_T  MSK_B61_T  MSK_B62_T  MSK_B63_T  MSK_B66_T  MSK_B67_T  MSK_B68_T  MSK_B69_T 
         1          1          1          1          1          1          1          1          1          1          1 
  MSK_B7_T  MSK_B70_T  MSK_B71_T  MSK_B72_T  MSK_B73_T  MSK_B74_T  MSK_B75_T  MSK_B76_T  MSK_B77_T  MSK_B78_T  MSK_B79_T 
         1          1          1          1          1          1          1          1          1          1          1 
  MSK_B8_T  MSK_B80_T  MSK_B81_T  MSK_B82_T  MSK_B83_T  MSK_B84_T  MSK_B85_T  MSK_B86_T  MSK_B87_T  MSK_B88_T  MSK_B89_T 
         1          1          1          1          1          1          1          1          1          1          1 
  MSK_B9_T  MSK_B90_T  MSK_B91_T  MSK_B92_T  MSK_B93_T  MSK_B94_T  MSK_B95_T  MSK_B96_T  MSK_B97_T  MSK_B98_T  MSK_B99_T 
         1          1          1          1          1          1          1          1          1          1          1 

alt_sample_name: 
  B1-T  B10-T B100-T B101-T B102-T B103-T B104-T B105-T B106-T B107-T B108-T B109-T  B11-T B110-T B111-T B112-T B113-T B114-T 
     1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1 
B115-T  B12-T  B13-T  B14-T  B15-T  B16-T  B17-T  B18-T   B2-T  B23-T  B24-T  B26-T  B27-T  B28-T  B29-T   B3-T  B31-T  B32-T 
     1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1 
 B33-T   B4-T  B42-T  B43-T  B44-T  B45-T  B46-T  B47-T  B48-T  B49-T   B5-T  B50-T  B51-T  B52-T  B53-T  B54-T  B55-T  B56-T 
     1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1 
 B57-T  B58-T  B59-T   B6-T  B60-T  B61-T  B62-T  B63-T  B66-T  B67-T  B68-T  B69-T   B7-T  B70-T  B71-T  B72-T  B73-T  B74-T 
     1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1 
 B75-T  B76-T  B77-T  B78-T  B79-T   B8-T  B80-T  B81-T  B82-T  B83-T  B84-T  B85-T  B86-T  B87-T  B88-T  B89-T   B9-T  B90-T 
     1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1 
 B91-T  B92-T  B93-T  B94-T  B95-T  B96-T  B97-T  B98-T  B99-T 
     1      1      1      1      1      1      1      1      1 

sample_type: 
tumor 
   99 

er: 
negative positive 
      42       57 

pgr: 
negative positive     NA's 
      55       43        1 

her2: 
positive     NA's 
      85       14 

tumor_size: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  1.100   2.450   3.200   3.624   4.300  10.000 

N: 
 0  1 
34 65 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  30.00   46.50   56.00   55.81   63.50   87.00 

dmfs_days: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    245    1279    1971    1888    2575    3924      17 

dmfs_status: 
deceased_or_recurrence    living_norecurrence                   NA's 
                    27                     55                     17 

batch: 
MSK 
 99 

uncurated_author_metadata: 
     samplename: MSK_B1_T///dataset: MSK///series: MSK///id: B1-T///er: 0///pgr: 0///her2: 1///size: 7.5///node: 1///age: 65///grade: NA///t.dmfs: 1307///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B10_T///dataset: MSK///series: MSK///id: B10-T///er: 0///pgr: 0///her2: 1///size: 6.5///node: 1///age: 52///grade: NA///t.dmfs: 2267///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
 samplename: MSK_B100_T///dataset: MSK///series: MSK///id: B100-T///er: 0///pgr: 0///her2: 1///size: 2.5///node: 1///age: 78///grade: NA///t.dmfs: 2124///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
 samplename: MSK_B101_T///dataset: MSK///series: MSK///id: B101-T///er: 0///pgr: 0///her2: 1///size: 3.5///node: 1///age: 47///grade: NA///t.dmfs: 3920///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B102_T///dataset: MSK///series: MSK///id: B102-T///er: 0///pgr: 0///her2: 1///size: 4///node: 1///age: 50///grade: NA///t.dmfs: 2646///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
 samplename: MSK_B103_T///dataset: MSK///series: MSK///id: B103-T///er: 0///pgr: 0///her2: 1///size: 1.3///node: 1///age: 57///grade: NA///t.dmfs: 1175///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
samplename: MSK_B104_T///dataset: MSK///series: MSK///id: B104-T///er: 1///pgr: 0///her2: 1///size: 2.95///node: 1///age: 49///grade: NA///t.dmfs: 2672///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
samplename: MSK_B105_T///dataset: MSK///series: MSK///id: B105-T///er: 1///pgr: 1///her2: 1///size: 2.95///node: 1///age: 49///grade: NA///t.dmfs: 2672///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
  samplename: MSK_B106_T///dataset: MSK///series: MSK///id: B106-T///er: 1///pgr: 1///her2: 1///size: 2.1///node: 1///age: 57///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
  samplename: MSK_B107_T///dataset: MSK///series: MSK///id: B107-T///er: 1///pgr: 1///her2: 1///size: 2.1///node: 1///age: 57///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
 samplename: MSK_B108_T///dataset: MSK///series: MSK///id: B108-T///er: 1///pgr: 1///her2: 1///size: 2.6///node: 1///age: 54///grade: NA///t.dmfs: 2387///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
    samplename: MSK_B109_T///dataset: MSK///series: MSK///id: B109-T///er: 0///pgr: 0///her2: 1///size: 6///node: 1///age: 32///grade: NA///t.dmfs: 504///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B11_T///dataset: MSK///series: MSK///id: B11-T///er: 1///pgr: 1///her2: 1///size: 1.8///node: 1///age: 52///grade: NA///t.dmfs: 1445///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
 samplename: MSK_B110_T///dataset: MSK///series: MSK///id: B110-T///er: 0///pgr: 0///her2: 1///size: 2.2///node: 1///age: 61///grade: NA///t.dmfs: 2646///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
 samplename: MSK_B111_T///dataset: MSK///series: MSK///id: B111-T///er: 1///pgr: 1///her2: 1///size: 1.3///node: 1///age: 46///grade: NA///t.dmfs: 1967///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B112_T///dataset: MSK///series: MSK///id: B112-T///er: 0///pgr: 1///her2: 1///size: 4///node: 1///age: 41///grade: NA///t.dmfs: 2548///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
  samplename: MSK_B113_T///dataset: MSK///series: MSK///id: B113-T///er: 1///pgr: 1///her2: 1///size: 1.5///node: 1///age: 49///grade: NA///t.dmfs: 938///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
  samplename: MSK_B114_T///dataset: MSK///series: MSK///id: B114-T///er: 0///pgr: 0///her2: 1///size: 3.8///node: 1///age: 61///grade: NA///t.dmfs: 672///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B115_T///dataset: MSK///series: MSK///id: B115-T///er: 0///pgr: 0///her2: 1///size: 5///node: 1///age: 40///grade: NA///t.dmfs: 2380///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
     samplename: MSK_B12_T///dataset: MSK///series: MSK///id: B12-T///er: 1///pgr: 1///her2: NA///size: 2///node: 1///age: 57///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B13_T///dataset: MSK///series: MSK///id: B13-T///er: 0///pgr: 0///her2: 1///size: 6.1///node: 1///age: 44///grade: NA///t.dmfs: 1894///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
     samplename: MSK_B14_T///dataset: MSK///series: MSK///id: B14-T///er: 1///pgr: 1///her2: 1///size: 5///node: 1///age: 61///grade: NA///t.dmfs: 2665///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
    samplename: MSK_B15_T///dataset: MSK///series: MSK///id: B15-T///er: 1///pgr: 0///her2: 1///size: 4.5///node: 1///age: 43///grade: NA///t.dmfs: 431///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
    samplename: MSK_B16_T///dataset: MSK///series: MSK///id: B16-T///er: 1///pgr: 1///her2: 1///size: 2.6///node: 1///age: 56///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
    samplename: MSK_B17_T///dataset: MSK///series: MSK///id: B17-T///er: 1///pgr: 1///her2: 1///size: 3.5///node: 1///age: 62///grade: NA///t.dmfs: 475///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B18_T///dataset: MSK///series: MSK///id: B18-T///er: 1///pgr: 1///her2: 1///size: 3.5///node: 1///age: 40///grade: NA///t.dmfs: 1613///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
     samplename: MSK_B2_T///dataset: MSK///series: MSK///id: B2-T///er: 0///pgr: 0///her2: 1///size: 4.5///node: 1///age: 59///grade: NA///t.dmfs: 1424///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
    samplename: MSK_B23_T///dataset: MSK///series: MSK///id: B23-T///er: 1///pgr: 1///her2: 1///size: 1.5///node: 1///age: 38///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
     samplename: MSK_B24_T///dataset: MSK///series: MSK///id: B24-T///er: 1///pgr: 1///her2: 1///size: 3///node: 1///age: 72///grade: NA///t.dmfs: 2292///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
     samplename: MSK_B26_T///dataset: MSK///series: MSK///id: B26-T///er: 0///pgr: 0///her2: 1///size: 4///node: 1///age: 54///grade: NA///t.dmfs: 1526///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B27_T///dataset: MSK///series: MSK///id: B27-T///er: 1///pgr: 1///her2: 1///size: 4.1///node: 1///age: 38///grade: NA///t.dmfs: 1343///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
     samplename: MSK_B28_T///dataset: MSK///series: MSK///id: B28-T///er: 1///pgr: 1///her2: 1///size: 6///node: 1///age: 54///grade: NA///t.dmfs: 1504///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
     samplename: MSK_B29_T///dataset: MSK///series: MSK///id: B29-T///er: 1///pgr: 1///her2: 1///size: 8///node: 1///age: 43///grade: NA///t.dmfs: 1183///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
       samplename: MSK_B3_T///dataset: MSK///series: MSK///id: B3-T///er: 1///pgr: 1///her2: 1///size: 6///node: 1///age: 59///grade: NA///t.dmfs: 1913///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
    samplename: MSK_B31_T///dataset: MSK///series: MSK///id: B31-T///er: 1///pgr: 1///her2: 1///size: 2.9///node: 1///age: 87///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
      samplename: MSK_B32_T///dataset: MSK///series: MSK///id: B32-T///er: 0///pgr: 0///her2: 1///size: 7///node: 1///age: 59///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
    samplename: MSK_B33_T///dataset: MSK///series: MSK///id: B33-T///er: 1///pgr: 1///her2: 1///size: 3.5///node: 1///age: 38///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
        samplename: MSK_B4_T///dataset: MSK///series: MSK///id: B4-T///er: 1///pgr: 1///her2: 1///size: 3///node: 1///age: 70///grade: NA///t.dmfs: 840///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B42_T///dataset: MSK///series: MSK///id: B42-T///er: 1///pgr: 1///her2: 1///size: 2.9///node: 0///age: 73///grade: NA///t.dmfs: 2785///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B43_T///dataset: MSK///series: MSK///id: B43-T///er: 0///pgr: 0///her2: 1///size: 4.5///node: 0///age: 52///grade: NA///t.dmfs: 1267///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
     samplename: MSK_B44_T///dataset: MSK///series: MSK///id: B44-T///er: 0///pgr: 0///her2: 1///size: 6///node: 0///age: 61///grade: NA///t.dmfs: 1445///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B45_T///dataset: MSK///series: MSK///id: B45-T///er: 1///pgr: 1///her2: 1///size: 3.5///node: 0///age: 68///grade: NA///t.dmfs: 1121///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B46_T///dataset: MSK///series: MSK///id: B46-T///er: 0///pgr: 0///her2: 1///size: 3.2///node: 0///age: 49///grade: NA///t.dmfs: 1135///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
    samplename: MSK_B47_T///dataset: MSK///series: MSK///id: B47-T///er: 0///pgr: 0///her2: 1///size: 6.5///node: 0///age: 70///grade: NA///t.dmfs: 343///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
     samplename: MSK_B48_T///dataset: MSK///series: MSK///id: B48-T///er: 1///pgr: 0///her2: 1///size: 3///node: 0///age: 54///grade: NA///t.dmfs: 2121///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
      samplename: MSK_B49_T///dataset: MSK///series: MSK///id: B49-T///er: 0///pgr: 0///her2: 1///size: 2///node: 0///age: 80///grade: NA///t.dmfs: 453///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
       samplename: MSK_B5_T///dataset: MSK///series: MSK///id: B5-T///er: 1///pgr: 0///her2: NA///size: 6///node: 1///age: 39///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
      samplename: MSK_B50_T///dataset: MSK///series: MSK///id: B50-T///er: 0///pgr: 0///her2: 1///size: 3///node: 0///age: 78///grade: NA///t.dmfs: 281///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B51_T///dataset: MSK///series: MSK///id: B51-T///er: 1///pgr: 0///her2: 1///size: 1.7///node: 0///age: 57///grade: NA///t.dmfs: 1270///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
     samplename: MSK_B52_T///dataset: MSK///series: MSK///id: B52-T///er: 1///pgr: 1///her2: 1///size: 3///node: 0///age: 40///grade: NA///t.dmfs: 3340///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B53_T///dataset: MSK///series: MSK///id: B53-T///er: 0///pgr: 0///her2: NA///size: 3.5///node: 0///age: 33///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
  samplename: MSK_B54_T///dataset: MSK///series: MSK///id: B54-T///er: 1///pgr: 1///her2: NA///size: 2.8///node: 0///age: 40///grade: NA///t.dmfs: 2504///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B55_T///dataset: MSK///series: MSK///id: B55-T///er: 1///pgr: 0///her2: 1///size: 2.2///node: 0///age: 75///grade: NA///t.dmfs: 2362///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B56_T///dataset: MSK///series: MSK///id: B56-T///er: 1///pgr: 1///her2: 1///size: 2.1///node: 0///age: 55///grade: NA///t.dmfs: 2584///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B57_T///dataset: MSK///series: MSK///id: B57-T///er: 0///pgr: 0///her2: 1///size: 2.8///node: 0///age: 57///grade: NA///t.dmfs: 2763///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
     samplename: MSK_B58_T///dataset: MSK///series: MSK///id: B58-T///er: 1///pgr: 0///her2: NA///size: 5///node: 0///age: 70///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
  samplename: MSK_B59_T///dataset: MSK///series: MSK///id: B59-T///er: 0///pgr: NA///her2: NA///size: 3.5///node: 0///age: 73///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
       samplename: MSK_B6_T///dataset: MSK///series: MSK///id: B6-T///er: 0///pgr: 0///her2: 1///size: 3///node: 1///age: 54///grade: NA///t.dmfs: 1226///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
      samplename: MSK_B60_T///dataset: MSK///series: MSK///id: B60-T///er: 0///pgr: 0///her2: 1///size: 4///node: 0///age: 46///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B61_T///dataset: MSK///series: MSK///id: B61-T///er: 0///pgr: 0///her2: 1///size: 3.5///node: 0///age: 73///grade: NA///t.dmfs: 1610///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
     samplename: MSK_B62_T///dataset: MSK///series: MSK///id: B62-T///er: 1///pgr: 0///her2: 1///size: 5///node: 0///age: 63///grade: NA///t.dmfs: 1909///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B63_T///dataset: MSK///series: MSK///id: B63-T///er: 1///pgr: 0///her2: 1///size: 3.2///node: 0///age: 64///grade: NA///t.dmfs: 2051///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
  samplename: MSK_B66_T///dataset: MSK///series: MSK///id: B66-T///er: 0///pgr: 0///her2: NA///size: 2.8///node: 0///age: 58///grade: NA///t.dmfs: 3924///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B67_T///dataset: MSK///series: MSK///id: B67-T///er: 0///pgr: 0///her2: 1///size: 1.8///node: 0///age: 32///grade: NA///t.dmfs: 2935///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B68_T///dataset: MSK///series: MSK///id: B68-T///er: 1///pgr: 1///her2: 1///size: 3.2///node: 0///age: 86///grade: NA///t.dmfs: 1975///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B69_T///dataset: MSK///series: MSK///id: B69-T///er: 0///pgr: 0///her2: 1///size: 1.5///node: 0///age: 49///grade: NA///t.dmfs: 2216///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
       samplename: MSK_B7_T///dataset: MSK///series: MSK///id: B7-T///er: 0///pgr: 0///her2: 1///size: 10///node: 1///age: 60///grade: NA///t.dmfs: 577///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
      samplename: MSK_B70_T///dataset: MSK///series: MSK///id: B70-T///er: 0///pgr: 0///her2: 1///size: 4///node: 0///age: 51///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B71_T///dataset: MSK///series: MSK///id: B71-T///er: 1///pgr: 1///her2: 1///size: 2.1///node: 0///age: 46///grade: NA///t.dmfs: 2008///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B72_T///dataset: MSK///series: MSK///id: B72-T///er: 0///pgr: 0///her2: NA///size: 1.1///node: 0///age: 63///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B73_T///dataset: MSK///series: MSK///id: B73-T///er: 0///pgr: 0///her2: 1///size: 3.2///node: 0///age: 73///grade: NA///t.dmfs: 1398///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
     samplename: MSK_B74_T///dataset: MSK///series: MSK///id: B74-T///er: 1///pgr: 1///her2: 1///size: 6///node: 0///age: 49///grade: NA///t.dmfs: 1376///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
    samplename: MSK_B75_T///dataset: MSK///series: MSK///id: B75-T///er: 0///pgr: 0///her2: NA///size: 5///node: 0///age: 47///grade: NA///t.dmfs: 1664///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
     samplename: MSK_B76_T///dataset: MSK///series: MSK///id: B76-T///er: 1///pgr: 1///her2: 1///size: 2///node: 0///age: 82///grade: NA///t.dmfs: 1449///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
    samplename: MSK_B77_T///dataset: MSK///series: MSK///id: B77-T///er: 0///pgr: 0///her2: NA///size: 3///node: 0///age: 73///grade: NA///t.dmfs: 1478///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
     samplename: MSK_B78_T///dataset: MSK///series: MSK///id: B78-T///er: 1///pgr: 1///her2: 1///size: 2///node: 1///age: 41///grade: NA///t.dmfs: 2646///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
  samplename: MSK_B79_T///dataset: MSK///series: MSK///id: B79-T///er: 1///pgr: 1///her2: NA///size: 2.4///node: 1///age: 81///grade: NA///t.dmfs: 1851///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
     samplename: MSK_B8_T///dataset: MSK///series: MSK///id: B8-T///er: 1///pgr: 1///her2: 1///size: 3.9///node: 1///age: 42///grade: NA///t.dmfs: 1675///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B80_T///dataset: MSK///series: MSK///id: B80-T///er: 1///pgr: 1///her2: 1///size: 2.2///node: 1///age: 49///grade: NA///t.dmfs: 2595///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B81_T///dataset: MSK///series: MSK///id: B81-T///er: 1///pgr: 1///her2: NA///size: 2.5///node: 1///age: 31///grade: NA///t.dmfs: 475///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B82_T///dataset: MSK///series: MSK///id: B82-T///er: 1///pgr: 1///her2: 1///size: 2.5///node: 1///age: 61///grade: NA///t.dmfs: 2657///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
     samplename: MSK_B83_T///dataset: MSK///series: MSK///id: B83-T///er: 1///pgr: 1///her2: 1///size: 4///node: 1///age: 39///grade: NA///t.dmfs: 2336///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B84_T///dataset: MSK///series: MSK///id: B84-T///er: 1///pgr: 1///her2: 1///size: 2.3///node: 1///age: 37///grade: NA///t.dmfs: 2373///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B85_T///dataset: MSK///series: MSK///id: B85-T///er: 1///pgr: 0///her2: 1///size: 4.3///node: 1///age: 57///grade: NA///t.dmfs: 3033///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B86_T///dataset: MSK///series: MSK///id: B86-T///er: 1///pgr: 0///her2: 1///size: 4.3///node: 1///age: 57///grade: NA///t.dmfs: 3033///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
    samplename: MSK_B87_T///dataset: MSK///series: MSK///id: B87-T///er: 1///pgr: 0///her2: 1///size: 2.5///node: 1///age: 51///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B88_T///dataset: MSK///series: MSK///id: B88-T///er: 1///pgr: 1///her2: 1///size: 4.8///node: 1///age: 61///grade: NA///t.dmfs: 2121///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
     samplename: MSK_B89_T///dataset: MSK///series: MSK///id: B89-T///er: 0///pgr: 0///her2: 1///size: 3///node: 1///age: 30///grade: NA///t.dmfs: 3190///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
      samplename: MSK_B9_T///dataset: MSK///series: MSK///id: B9-T///er: 1///pgr: 0///her2: 1///size: 3.9///node: 1///age: 74///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
  samplename: MSK_B90_T///dataset: MSK///series: MSK///id: B90-T///er: 1///pgr: 1///her2: NA///size: 4.3///node: 1///age: 34///grade: NA///t.dmfs: 2478///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
    samplename: MSK_B91_T///dataset: MSK///series: MSK///id: B91-T///er: 1///pgr: 0///her2: 1///size: 10///node: 1///age: 84///grade: NA///t.dmfs: 2624///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B92_T///dataset: MSK///series: MSK///id: B92-T///er: 1///pgr: 1///her2: 1///size: 2.5///node: 1///age: 72///grade: NA///t.dmfs: 2321///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B93_T///dataset: MSK///series: MSK///id: B93-T///er: 0///pgr: 0///her2: 1///size: 2.3///node: 1///age: 50///grade: NA///t.dmfs: 2354///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
      samplename: MSK_B94_T///dataset: MSK///series: MSK///id: B94-T///er: 0///pgr: 0///her2: 1///size: 4///node: 1///age: 48///grade: NA///t.dmfs: 245///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
    samplename: MSK_B95_T///dataset: MSK///series: MSK///id: B95-T///er: 0///pgr: 0///her2: 1///size: 3.8///node: 1///age: 61///grade: NA///t.dmfs: 672///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
   samplename: MSK_B96_T///dataset: MSK///series: MSK///id: B96-T///er: 1///pgr: 1///her2: 1///size: 2.9///node: 1///age: 65///grade: NA///t.dmfs: 2471///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
     samplename: MSK_B97_T///dataset: MSK///series: MSK///id: B97-T///er: 0///pgr: 0///her2: 1///size: 2///node: 1///age: 64///grade: NA///t.dmfs: 2869///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
     samplename: MSK_B98_T///dataset: MSK///series: MSK///id: B98-T///er: 0///pgr: 0///her2: 1///size: 2///node: 1///age: 58///grade: NA///t.dmfs: 1226///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 
  samplename: MSK_B99_T///dataset: MSK///series: MSK///id: B99-T///er: 1///pgr: 0///her2: NA///size: 4.5///node: 1///age: 54///grade: NA///t.dmfs: 2624///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                            1 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse10510 }
\keyword{datasets}
