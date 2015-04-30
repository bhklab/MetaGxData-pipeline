\name{ STNO2_eset }
\alias{ STNO2_eset }
\docType{data}
\title{  }
\usage{data( STNO2_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=18498629; http://www.ncbi.nlm.nih.gov/pubmed/?term=20697068; http://www.ncbi.nlm.nih.gov/pubmed/?term=18821012 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6532; http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE12093; http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17700; http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17705 
  PMIDs: 18498629, 20697068, 18821012 
  No abstract available.
  notes:
   summary:     
      
   version:     
      2015-04-27 19:31:37
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: IMAGE:1020315 IMAGE:1030271 ... IMAGE:971399 (3663 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 3663 features, 118 samples
Platform type: 
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

records   n.max n.start  events  median 0.95LCL 0.95UCL 
 118.00  118.00  118.00   46.00    4.67    3.34      NA 

--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
         STNO2_BC-A       STNO2_BC-HBC2       STNO2_BC-HBC3    STNO2_BC-HBC4-T1       STNO2_BC-HBC5       STNO2_BC-HBC6 
                  1                   1                   1                   1                   1                   1 
    STNO2_BC102B-BE     STNO2_BC104A-BE     STNO2_BC105A-BE     STNO2_BC106B-BE     STNO2_BC107B-BE     STNO2_BC108A-BE 
                  1                   1                   1                   1                   1                   1 
      STNO2_BC11-FA     STNO2_BC110B-BE     STNO2_BC111A-BE     STNO2_BC111B-BE     STNO2_BC112B-BE     STNO2_BC114A-BE 
                  1                   1                   1                   1                   1                   1 
    STNO2_BC115B-BE     STNO2_BC116A-BE     STNO2_BC117A-BE     STNO2_BC118B-BE     STNO2_BC119A-BE     STNO2_BC120A-BE 
                  1                   1                   1                   1                   1                   1 
    STNO2_BC121B-BE     STNO2_BC123B-BE     STNO2_BC124A-BE        STNO2_BC1257     STNO2_BC125A-BE        STNO2_BC1369 
                  1                   1                   1                   1                   1                   1 
         STNO2_BC14          STNO2_BC16          STNO2_BC17          STNO2_BC18           STNO2_BC2       STNO2_BC20-FA 
                  1                   1                   1                   1                   1                   1 
    STNO2_BC201B-BE     STNO2_BC205A-BE     STNO2_BC206A-BE     STNO2_BC208A-BE     STNO2_BC210B-AF     STNO2_BC213B-BE 
                  1                   1                   1                   1                   1                   1 
    STNO2_BC214B-BE          STNO2_BC23          STNO2_BC24     STNO2_BC303B-BE     STNO2_BC305A-BE     STNO2_BC307B-BE 
                  1                   1                   1                   1                   1                   1 
    STNO2_BC308B-BE     STNO2_BC309A-BE        STNO2_BC31-0        STNO2_BC35-0       STNO2_BC37-FA          STNO2_BC38 
                  1                   1                   1                   1                   1                   1 
      STNO2_BC4-LN4          STNO2_BC40     STNO2_BC402B-BE     STNO2_BC404B-BE     STNO2_BC405A-BE STNO2_BC406A-2ndTUM 
                  1                   1                   1                   1                   1                   1 
         STNO2_BC44          STNO2_BC45     STNO2_BC46-LN46        STNO2_BC48-0     STNO2_BC503B-BE           STNO2_BC6 
                  1                   1                   1                   1                   1                   1 
    STNO2_BC601A-BE     STNO2_BC605B-BE     STNO2_BC606B-AF     STNO2_BC608B-BE     STNO2_BC610A-BE     STNO2_BC702B-BE 
                  1                   1                   1                   1                   1                   1 
    STNO2_BC703B-BE     STNO2_BC704B-AF     STNO2_BC706A-BE     STNO2_BC708B-BE     STNO2_BC709B-BE     STNO2_BC710A-BE 
                  1                   1                   1                   1                   1                   1 
    STNO2_BC711B-BE     STNO2_BC713A-BE         STNO2_BC790     STNO2_BC805A-BE     STNO2_BC807A-BE     STNO2_BC808A-BE 
                  1                   1                   1                   1                   1                   1 
     STNO2_FU_01-BE      STNO2_FU_02-BE      STNO2_FU_04-BE      STNO2_FU_05-BE      STNO2_FU_06-BE      STNO2_FU_07-BE 
                  1                   1                   1                   1                   1                   1 
     STNO2_FU_08-BE      STNO2_FU_09-BE      STNO2_FU_10-BE      STNO2_FU_11-BE      STNO2_FU_12-BE      STNO2_FU_14-BE 
                  1                   1                   1                   1                   1                   1 
     STNO2_FU_15-BE      STNO2_FU_16-BE      STNO2_FU_17-AF             (Other) 
                  1                   1                   1                  19 

alt_sample_name: 
         BC-A       BC-HBC2       BC-HBC3    BC-HBC4-T1       BC-HBC5       BC-HBC6     BC102B-BE     BC104A-BE     BC105A-BE 
            1             1             1             1             1             1             1             1             1 
    BC106B-BE     BC107B-BE     BC108A-BE       BC11-FA     BC110B-BE     BC111A-BE     BC111B-BE     BC112B-BE     BC114A-BE 
            1             1             1             1             1             1             1             1             1 
    BC115B-BE     BC116A-BE     BC117A-BE     BC118B-BE     BC119A-BE     BC120A-BE     BC121B-BE     BC123B-BE     BC124A-BE 
            1             1             1             1             1             1             1             1             1 
       BC1257     BC125A-BE        BC1369          BC14          BC16          BC17          BC18           BC2       BC20-FA 
            1             1             1             1             1             1             1             1             1 
    BC201B-BE     BC205A-BE     BC206A-BE     BC208A-BE     BC210B-AF     BC213B-BE     BC214B-BE          BC23          BC24 
            1             1             1             1             1             1             1             1             1 
    BC303B-BE     BC305A-BE     BC307B-BE     BC308B-BE     BC309A-BE        BC31-0        BC35-0       BC37-FA          BC38 
            1             1             1             1             1             1             1             1             1 
      BC4-LN4          BC40     BC402B-BE     BC404B-BE     BC405A-BE BC406A-2ndTUM          BC44          BC45     BC46-LN46 
            1             1             1             1             1             1             1             1             1 
       BC48-0     BC503B-BE           BC6     BC601A-BE     BC605B-BE     BC606B-AF     BC608B-BE     BC610A-BE     BC702B-BE 
            1             1             1             1             1             1             1             1             1 
    BC703B-BE     BC704B-AF     BC706A-BE     BC708B-BE     BC709B-BE     BC710A-BE     BC711B-BE     BC713A-BE         BC790 
            1             1             1             1             1             1             1             1             1 
    BC805A-BE     BC807A-BE     BC808A-BE      FU_01-BE      FU_02-BE      FU_04-BE      FU_05-BE      FU_06-BE      FU_07-BE 
            1             1             1             1             1             1             1             1             1 
     FU_08-BE      FU_09-BE      FU_10-BE      FU_11-BE      FU_12-BE      FU_14-BE      FU_15-BE      FU_16-BE      FU_17-AF 
            1             1             1             1             1             1             1             1             1 
      (Other) 
           19 

sample_type: 
tumor 
  118 

er: 
negative positive     NA's 
      31       82        5 

tumor_size: 
   1    2    3    4 NA's 
   6   13   62   32    5 

N: 
   0    1 NA's 
  34   79    5 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  21.00   46.25   58.00   58.47   71.75   85.00 

grade: 
   1    2    3 NA's 
  11   49   53    5 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  122.0   396.0   761.0   927.9  1233.0  2800.0      23 

recurrence_status: 
deceased_or_recurrence    living_norecurrence 
                    60                     58 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
     91     426     898    1019    1392    5722 

vital_status: 
deceased   living 
      46       72 

treatment: 
  chemotherapy hormonotherapy      untreated 
            23             73             22 

batch: 
STNO2 
  118 

uncurated_author_metadata: 
                samplename: STNO2_BC-A///dataset: STNO2///series: STNO2///id: BC-A///er: 1///pgr: NA///her2: NA///size: 2///node: 1///age: 38///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 487///e.rfs: 1///t.os: 2892///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                           1 
          samplename: STNO2_BC-HBC2///dataset: STNO2///series: STNO2///id: BC-HBC2///er: 0///pgr: NA///her2: NA///size: 1///node: 0///age: 39///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 457///e.rfs: 1///t.os: 1644///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                           1 
          samplename: STNO2_BC-HBC3///dataset: STNO2///series: STNO2///id: BC-HBC3///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 31///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 883///e.rfs: 1///t.os: 1431///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
   samplename: STNO2_BC-HBC4-T1///dataset: STNO2///series: STNO2///id: BC-HBC4-T1///er: 1///pgr: NA///her2: NA///size: 2///node: 0///age: 38///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 1826///e.rfs: 0///t.os: 1826///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                           1 
           samplename: STNO2_BC-HBC5///dataset: STNO2///series: STNO2///id: BC-HBC5///er: 0///pgr: NA///her2: NA///size: 2///node: 0///age: 28///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 274///e.rfs: 1///t.os: 517///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                           1 
      samplename: STNO2_BC-HBC6///dataset: STNO2///series: STNO2///id: BC-HBC6///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: 41///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: 1735///e.rfs: 0///t.os: 1735///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                           1 
     samplename: STNO2_BC102B-BE///dataset: STNO2///series: STNO2///id: BC102B-BE///er: 1///pgr: NA///her2: NA///size: 4///node: 0///age: 75///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 1248///e.rfs: 0///t.os: 1248///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
     samplename: STNO2_BC104A-BE///dataset: STNO2///series: STNO2///id: BC104A-BE///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 68///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2800///e.rfs: 0///t.os: 2800///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_BC105A-BE///dataset: STNO2///series: STNO2///id: BC105A-BE///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 68///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 274///e.rfs: 1///t.os: 548///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
     samplename: STNO2_BC106B-BE///dataset: STNO2///series: STNO2///id: BC106B-BE///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 78///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 1704///e.rfs: 0///t.os: 1704///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
     samplename: STNO2_BC107B-BE///dataset: STNO2///series: STNO2///id: BC107B-BE///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 66///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 1187///e.rfs: 1///t.os: 2709///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_BC108A-BE///dataset: STNO2///series: STNO2///id: BC108A-BE///er: 0///pgr: NA///her2: NA///size: 3///node: 1///age: 47///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 670///e.rfs: 1///t.os: 670///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                           1 
     samplename: STNO2_BC11-FA///dataset: STNO2///series: STNO2///id: BC11-FA///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: 21///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: 1400///e.rfs: 0///t.os: 1400///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_BC110B-BE///dataset: STNO2///series: STNO2///id: BC110B-BE///er: 1///pgr: NA///her2: NA///size: 4///node: 0///age: 58///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 700///e.rfs: 1///t.os: 913///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
        samplename: STNO2_BC111A-BE///dataset: STNO2///series: STNO2///id: BC111A-BE///er: 1///pgr: NA///her2: NA///size: 4///node: 1///age: 67///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: 1///t.os: 457///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
        samplename: STNO2_BC111B-BE///dataset: STNO2///series: STNO2///id: BC111B-BE///er: 1///pgr: NA///her2: NA///size: 4///node: 1///age: 74///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: 1///t.os: 609///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_BC112B-BE///dataset: STNO2///series: STNO2///id: BC112B-BE///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 75///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 213///e.rfs: 1///t.os: 244///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_BC114A-BE///dataset: STNO2///series: STNO2///id: BC114A-BE///er: 1///pgr: NA///her2: NA///size: 4///node: 1///age: 82///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 335///e.rfs: 1///t.os: 365///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
     samplename: STNO2_BC115B-BE///dataset: STNO2///series: STNO2///id: BC115B-BE///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 47///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2405///e.rfs: 0///t.os: 2405///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
     samplename: STNO2_BC116A-BE///dataset: STNO2///series: STNO2///id: BC116A-BE///er: 1///pgr: NA///her2: NA///size: 4///node: 1///age: 74///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2374///e.rfs: 0///t.os: 2374///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
     samplename: STNO2_BC117A-BE///dataset: STNO2///series: STNO2///id: BC117A-BE///er: 1///pgr: NA///her2: NA///size: 4///node: 0///age: 76///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2435///e.rfs: 0///t.os: 2435///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
     samplename: STNO2_BC118B-BE///dataset: STNO2///series: STNO2///id: BC118B-BE///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 66///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 1004///e.rfs: 0///t.os: 1004///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_BC119A-BE///dataset: STNO2///series: STNO2///id: BC119A-BE///er: 0///pgr: NA///her2: NA///size: 4///node: 1///age: 79///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 365///e.rfs: 1///t.os: 426///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                           1 
      samplename: STNO2_BC120A-BE///dataset: STNO2///series: STNO2///id: BC120A-BE///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 43///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 304///e.rfs: 1///t.os: 1035///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_BC121B-BE///dataset: STNO2///series: STNO2///id: BC121B-BE///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 62///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 700///e.rfs: 1///t.os: 730///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
      samplename: STNO2_BC123B-BE///dataset: STNO2///series: STNO2///id: BC123B-BE///er: 1///pgr: NA///her2: NA///size: 3///node: 0///age: 72///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 274///e.rfs: 1///t.os: 1309///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
     samplename: STNO2_BC124A-BE///dataset: STNO2///series: STNO2///id: BC124A-BE///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 77///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 1796///e.rfs: 1///t.os: 2252///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
           samplename: STNO2_BC1257///dataset: STNO2///series: STNO2///id: BC1257///er: NA///pgr: NA///her2: NA///size: NA///node: 1///age: 46///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: 1///t.os: 5722///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                           1 
     samplename: STNO2_BC125A-BE///dataset: STNO2///series: STNO2///id: BC125A-BE///er: 1///pgr: NA///her2: NA///size: 3///node: 0///age: 46///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2252///e.rfs: 0///t.os: 2252///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
              samplename: STNO2_BC1369///dataset: STNO2///series: STNO2///id: BC1369///er: 0///pgr: NA///her2: NA///size: 3///node: 1///age: 54///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: 1///t.os: 365///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                           1 
               samplename: STNO2_BC14///dataset: STNO2///series: STNO2///id: BC14///er: 0///pgr: NA///her2: NA///size: 2///node: 0///age: 57///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 1218///e.rfs: 0///t.os: 1218///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                           1 
                 samplename: STNO2_BC16///dataset: STNO2///series: STNO2///id: BC16///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 77///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 426///e.rfs: 1///t.os: 426///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
               samplename: STNO2_BC17///dataset: STNO2///series: STNO2///id: BC17///er: 1///pgr: NA///her2: NA///size: 2///node: 1///age: 53///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 1431///e.rfs: 0///t.os: 1431///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                           1 
               samplename: STNO2_BC18///dataset: STNO2///series: STNO2///id: BC18///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 35///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1187///e.rfs: 0///t.os: 1187///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
                   samplename: STNO2_BC2///dataset: STNO2///series: STNO2///id: BC2///er: 0///pgr: NA///her2: NA///size: 3///node: 1///age: 44///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 122///e.rfs: 1///t.os: 152///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_BC20-FA///dataset: STNO2///series: STNO2///id: BC20-FA///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: 40///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: 548///e.rfs: 0///t.os: 548///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                           1 
     samplename: STNO2_BC201B-BE///dataset: STNO2///series: STNO2///id: BC201B-BE///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 55///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2222///e.rfs: 0///t.os: 2222///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_BC205A-BE///dataset: STNO2///series: STNO2///id: BC205A-BE///er: 0///pgr: NA///her2: NA///size: 3///node: 1///age: 71///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 365///e.rfs: 1///t.os: 670///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                           1 
     samplename: STNO2_BC206A-BE///dataset: STNO2///series: STNO2///id: BC206A-BE///er: 1///pgr: NA///her2: NA///size: 3///node: 0///age: 85///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2039///e.rfs: 0///t.os: 2039///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
        samplename: STNO2_BC208A-BE///dataset: STNO2///series: STNO2///id: BC208A-BE///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 49///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: 1///t.os: 274///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
     samplename: STNO2_BC210B-AF///dataset: STNO2///series: STNO2///id: BC210B-AF///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 45///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2100///e.rfs: 0///t.os: 2100///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
    samplename: STNO2_BC213B-BE///dataset: STNO2///series: STNO2///id: BC213B-BE///er: 1///pgr: NA///her2: NA///size: 3///node: 0///age: 51///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: 2009///e.rfs: 0///t.os: 2009///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_BC214B-BE///dataset: STNO2///series: STNO2///id: BC214B-BE///er: 0///pgr: NA///her2: NA///size: 3///node: 0///age: 54///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 122///e.rfs: 1///t.os: 213///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                           1 
                 samplename: STNO2_BC23///dataset: STNO2///series: STNO2///id: BC23///er: 0///pgr: NA///her2: NA///size: 1///node: 0///age: 34///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 761///e.rfs: 0///t.os: 761///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                           1 
                 samplename: STNO2_BC24///dataset: STNO2///series: STNO2///id: BC24///er: 1///pgr: NA///her2: NA///size: 2///node: 1///age: 69///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 700///e.rfs: 0///t.os: 700///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                           1 
     samplename: STNO2_BC303B-BE///dataset: STNO2///series: STNO2///id: BC303B-BE///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 53///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1826///e.rfs: 0///t.os: 1826///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
      samplename: STNO2_BC305A-BE///dataset: STNO2///series: STNO2///id: BC305A-BE///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 32///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 700///e.rfs: 1///t.os: 1035///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_BC307B-BE///dataset: STNO2///series: STNO2///id: BC307B-BE///er: 1///pgr: NA///her2: NA///size: 4///node: 1///age: 71///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 304///e.rfs: 1///t.os: 487///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
     samplename: STNO2_BC308B-BE///dataset: STNO2///series: STNO2///id: BC308B-BE///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 64///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 1857///e.rfs: 0///t.os: 1857///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
        samplename: STNO2_BC309A-BE///dataset: STNO2///series: STNO2///id: BC309A-BE///er: 0///pgr: NA///her2: NA///size: 4///node: 1///age: 46///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: 1///t.os: 304///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                           1 
             samplename: STNO2_BC31-0///dataset: STNO2///series: STNO2///id: BC31-0///er: 1///pgr: NA///her2: NA///size: 1///node: 1///age: 48///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 487///e.rfs: 0///t.os: 487///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                           1 
             samplename: STNO2_BC35-0///dataset: STNO2///series: STNO2///id: BC35-0///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 52///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 487///e.rfs: 1///t.os: 578///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_BC37-FA///dataset: STNO2///series: STNO2///id: BC37-FA///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: 34///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: 365///e.rfs: 0///t.os: 365///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                           1 
                 samplename: STNO2_BC38///dataset: STNO2///series: STNO2///id: BC38///er: 1///pgr: NA///her2: NA///size: 1///node: 1///age: 45///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 548///e.rfs: 0///t.os: 548///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                           1 
           samplename: STNO2_BC4-LN4///dataset: STNO2///series: STNO2///id: BC4-LN4///er: 1///pgr: NA///her2: NA///size: 2///node: 1///age: 79///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 396///e.rfs: 0///t.os: 396///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                           1 
                 samplename: STNO2_BC40///dataset: STNO2///series: STNO2///id: BC40///er: 1///pgr: NA///her2: NA///size: 2///node: 1///age: 54///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 396///e.rfs: 0///t.os: 396///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                           1 
     samplename: STNO2_BC402B-BE///dataset: STNO2///series: STNO2///id: BC402B-BE///er: 1///pgr: NA///her2: NA///size: 3///node: 0///age: 58///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1735///e.rfs: 0///t.os: 1735///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
      samplename: STNO2_BC404B-BE///dataset: STNO2///series: STNO2///id: BC404B-BE///er: 1///pgr: NA///her2: NA///size: 3///node: 0///age: 44///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 761///e.rfs: 1///t.os: 1248///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
     samplename: STNO2_BC405A-BE///dataset: STNO2///series: STNO2///id: BC405A-BE///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 61///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 1704///e.rfs: 0///t.os: 1704///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
samplename: STNO2_BC406A-2ndTUM///dataset: STNO2///series: STNO2///id: BC406A-2ndTUM///er: 0///pgr: NA///her2: NA///size: 4///node: 1///age: 66///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: 1///t.os: 700///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                           1 
                 samplename: STNO2_BC44///dataset: STNO2///series: STNO2///id: BC44///er: 0///pgr: NA///her2: NA///size: 1///node: 0///age: 62///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 244///e.rfs: 0///t.os: 244///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                           1 
                 samplename: STNO2_BC45///dataset: STNO2///series: STNO2///id: BC45///er: 0///pgr: NA///her2: NA///size: 2///node: 1///age: 33///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 335///e.rfs: 0///t.os: 335///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_BC46-LN46///dataset: STNO2///series: STNO2///id: BC46-LN46///er: 0///pgr: NA///her2: NA///size: 1///node: 1///age: 74///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 274///e.rfs: 0///t.os: 274///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                           1 
            samplename: STNO2_BC48-0///dataset: STNO2///series: STNO2///id: BC48-0///er: NA///pgr: NA///her2: NA///size: 3///node: 0///age: 54///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 122///e.rfs: 0///t.os: 122///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                           1 
     samplename: STNO2_BC503B-BE///dataset: STNO2///series: STNO2///id: BC503B-BE///er: 1///pgr: NA///her2: NA///size: 3///node: 0///age: 79///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 1218///e.rfs: 0///t.os: 1218///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
                 samplename: STNO2_BC6///dataset: STNO2///series: STNO2///id: BC6///er: 1///pgr: NA///her2: NA///size: 2///node: 0///age: 64///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1187///e.rfs: 0///t.os: 1187///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_BC601A-BE///dataset: STNO2///series: STNO2///id: BC601A-BE///er: 0///pgr: NA///her2: NA///size: 3///node: 0///age: 56///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: 1///t.os: 1218///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_BC605B-BE///dataset: STNO2///series: STNO2///id: BC605B-BE///er: 0///pgr: NA///her2: NA///size: 3///node: 1///age: 58///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 274///e.rfs: 1///t.os: 335///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_BC606B-AF///dataset: STNO2///series: STNO2///id: BC606B-AF///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 52///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 244///e.rfs: 1///t.os: 396///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_BC608B-BE///dataset: STNO2///series: STNO2///id: BC608B-BE///er: 1///pgr: NA///her2: NA///size: 4///node: 0///age: 79///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: 1///t.os: 1370///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
     samplename: STNO2_BC610A-BE///dataset: STNO2///series: STNO2///id: BC610A-BE///er: 1///pgr: NA///her2: NA///size: 4///node: 1///age: 79///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 1187///e.rfs: 0///t.os: 1187///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
     samplename: STNO2_BC702B-BE///dataset: STNO2///series: STNO2///id: BC702B-BE///er: 1///pgr: NA///her2: NA///size: 4///node: 1///age: 64///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 1096///e.rfs: 0///t.os: 1096///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
        samplename: STNO2_BC703B-BE///dataset: STNO2///series: STNO2///id: BC703B-BE///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 64///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: 1///t.os: 883///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_BC704B-AF///dataset: STNO2///series: STNO2///id: BC704B-AF///er: 1///pgr: NA///her2: NA///size: 3///node: 0///age: 40///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 913///e.rfs: 0///t.os: 913///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_BC706A-BE///dataset: STNO2///series: STNO2///id: BC706A-BE///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 57///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 974///e.rfs: 0///t.os: 974///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_BC708B-BE///dataset: STNO2///series: STNO2///id: BC708B-BE///er: 1///pgr: NA///her2: NA///size: 3///node: 0///age: 47///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 304///e.rfs: 1///t.os: 517///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
      samplename: STNO2_BC709B-BE///dataset: STNO2///series: STNO2///id: BC709B-BE///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 73///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 609///e.rfs: 1///t.os: 1004///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_BC710A-BE///dataset: STNO2///series: STNO2///id: BC710A-BE///er: 1///pgr: NA///her2: NA///size: 4///node: 1///age: 65///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 883///e.rfs: 0///t.os: 883///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_BC711B-BE///dataset: STNO2///series: STNO2///id: BC711B-BE///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 44///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 944///e.rfs: 0///t.os: 944///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_BC713A-BE///dataset: STNO2///series: STNO2///id: BC713A-BE///er: 1///pgr: NA///her2: NA///size: 4///node: 1///age: 76///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: 1///t.os: 1004///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
                samplename: STNO2_BC790///dataset: STNO2///series: STNO2///id: BC790///er: 0///pgr: NA///her2: NA///size: 2///node: 0///age: 80///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: 1///t.os: 183///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_BC805A-BE///dataset: STNO2///series: STNO2///id: BC805A-BE///er: 0///pgr: NA///her2: NA///size: 3///node: 1///age: 48///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 791///e.rfs: 0///t.os: 791///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_BC807A-BE///dataset: STNO2///series: STNO2///id: BC807A-BE///er: 1///pgr: NA///her2: NA///size: 3///node: 0///age: 55///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 913///e.rfs: 0///t.os: 913///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_BC808A-BE///dataset: STNO2///series: STNO2///id: BC808A-BE///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 37///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 670///e.rfs: 0///t.os: 670///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_FU_01-BE///dataset: STNO2///series: STNO2///id: FU_01-BE///er: 0///pgr: NA///her2: NA///size: 4///node: 1///age: 58///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1644///e.rfs: 0///t.os: 1644///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_FU_02-BE///dataset: STNO2///series: STNO2///id: FU_02-BE///er: 1///pgr: NA///her2: NA///size: 4///node: 1///age: 44///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1613///e.rfs: 0///t.os: 1613///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
         samplename: STNO2_FU_04-BE///dataset: STNO2///series: STNO2///id: FU_04-BE///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 37///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 396///e.rfs: 1///t.os: 548///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_FU_05-BE///dataset: STNO2///series: STNO2///id: FU_05-BE///er: 0///pgr: NA///her2: NA///size: 4///node: 1///age: 51///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1218///e.rfs: 1///t.os: 1583///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                           1 
         samplename: STNO2_FU_06-BE///dataset: STNO2///series: STNO2///id: FU_06-BE///er: 0///pgr: NA///her2: NA///size: 4///node: 1///age: 42///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 396///e.rfs: 1///t.os: 913///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_FU_07-BE///dataset: STNO2///series: STNO2///id: FU_07-BE///er: 0///pgr: NA///her2: NA///size: 3///node: 1///age: 51///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 1400///e.rfs: 0///t.os: 1400///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                           1 
          samplename: STNO2_FU_08-BE///dataset: STNO2///series: STNO2///id: FU_08-BE///er: 1///pgr: NA///her2: NA///size: 4///node: 0///age: 59///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: 1///t.os: 304///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_FU_09-BE///dataset: STNO2///series: STNO2///id: FU_09-BE///er: 1///pgr: NA///her2: NA///size: 4///node: 0///age: 63///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 1157///e.rfs: 0///t.os: 1157///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_FU_10-BE///dataset: STNO2///series: STNO2///id: FU_10-BE///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 74///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 1126///e.rfs: 1///t.os: 1400///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_FU_11-BE///dataset: STNO2///series: STNO2///id: FU_11-BE///er: 1///pgr: NA///her2: NA///size: 3///node: 0///age: 59///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 1157///e.rfs: 1///t.os: 1400///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
          samplename: STNO2_FU_12-BE///dataset: STNO2///series: STNO2///id: FU_12-BE///er: 0///pgr: NA///her2: NA///size: 4///node: 1///age: 73///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: 1///t.os: 152///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                           1 
       samplename: STNO2_FU_14-BE///dataset: STNO2///series: STNO2///id: FU_14-BE///er: 0///pgr: NA///her2: NA///size: 3///node: 0///age: 70///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 1218///e.rfs: 0///t.os: 1218///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                           1 
        samplename: STNO2_FU_15-BE///dataset: STNO2///series: STNO2///id: FU_15-BE///er: 0///pgr: NA///her2: NA///size: 4///node: 0///age: 70///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 974///e.rfs: 1///t.os: 1004///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                           1 
         samplename: STNO2_FU_16-BE///dataset: STNO2///series: STNO2///id: FU_16-BE///er: 1///pgr: NA///her2: NA///size: 4///node: 1///age: 74///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: 1///t.os: 1126///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
         samplename: STNO2_FU_17-AF///dataset: STNO2///series: STNO2///id: FU_17-AF///er: 1///pgr: NA///her2: NA///size: 4///node: 1///age: 70///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 122///e.rfs: 0///t.os: 122///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                           1 
                                                                                                                                                                                                                                                     (Other) 
                                                                                                                                                                                                                                                          19 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6532; http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE12093; http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17700; http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17705 }
\keyword{datasets}
