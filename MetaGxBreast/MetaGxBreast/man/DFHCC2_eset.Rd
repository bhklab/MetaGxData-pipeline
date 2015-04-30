\name{ DFHCC2_eset }
\alias{ DFHCC2_eset }
\docType{data}
\title{  }
\description{ Test the efficacy of treating TNBC with neoadjuvant cisplatin; explore biomarkers to identify predictors of response }
\usage{data( DFHCC2_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2834466/ 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18864 
  PMIDs: 20100965 

  Abstract: A 16 word abstract is available. Use 'abstract' method.
  notes:
   summary:     
      A subset of the patients experienced a response induced by cisplatin and biomarkers were identified that could predict re
sponse to cisplatin.  
   version:     
      2015-04-27 18:47:07
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: 1007_s_at 1053_at ... AFFX-HUMISGF3A/M97935_MB_at (42447 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 42447 features, 84 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
       DFHCC2_P1       DFHCC2_P10       DFHCC2_P11       DFHCC2_P12       DFHCC2_P14       DFHCC2_P15       DFHCC2_P16 
               1                1                1                1                1                1                1 
      DFHCC2_P17        DFHCC2_P2       DFHCC2_P21       DFHCC2_P22       DFHCC2_P23       DFHCC2_P24       DFHCC2_P25 
               1                1                1                1                1                1                1 
      DFHCC2_P26       DFHCC2_P27       DFHCC2_P28       DFHCC2_P29        DFHCC2_P3        DFHCC2_P4        DFHCC2_P5 
               1                1                1                1                1                1                1 
       DFHCC2_P6        DFHCC2_P7        DFHCC2_P8      DFHCC2_REF1     DFHCC2_REF10     DFHCC2_REF11    DFHCC2_REF113 
               1                1                1                1                1                1                1 
   DFHCC2_REF115    DFHCC2_REF116 DFHCC2_REF116rep    DFHCC2_REF117    DFHCC2_REF118     DFHCC2_REF12    DFHCC2_REF120 
               1                1                1                1                1                1                1 
   DFHCC2_REF122    DFHCC2_REF123    DFHCC2_REF124    DFHCC2_REF125  DFHCC2_REF12rep    DFHCC2_REF134    DFHCC2_REF136 
               1                1                1                1                1                1                1 
   DFHCC2_REF140    DFHCC2_REF141    DFHCC2_REF145    DFHCC2_REF147     DFHCC2_REF17    DFHCC2_REF187   DFHCC2_REF1rep 
               1                1                1                1                1                1                1 
     DFHCC2_REF2     DFHCC2_REF20     DFHCC2_REF21  DFHCC2_REF21rep     DFHCC2_REF22     DFHCC2_REF25  DFHCC2_REF25rep 
               1                1                1                1                1                1                1 
    DFHCC2_REF28     DFHCC2_REF29     DFHCC2_REF30     DFHCC2_REF35     DFHCC2_REF38  DFHCC2_REF38rep     DFHCC2_REF40 
               1                1                1                1                1                1                1 
    DFHCC2_REF43     DFHCC2_REF44     DFHCC2_REF47      DFHCC2_REF5     DFHCC2_REF54     DFHCC2_REF55  DFHCC2_REF55rep 
               1                1                1                1                1                1                1 
    DFHCC2_REF63  DFHCC2_REF63rep     DFHCC2_REF66      DFHCC2_REF7     DFHCC2_REF70     DFHCC2_REF71     DFHCC2_REF72 
               1                1                1                1                1                1                1 
    DFHCC2_REF80  DFHCC2_REF80rep     DFHCC2_REF85     DFHCC2_REF90     DFHCC2_REF91     DFHCC2_REF96     DFHCC2_REF98 
               1                1                1                1                1                1                1 

unique_patient_ID: 
  P1  P10  P11  P12  P14  P15  P16  P17   P2  P21  P22  P23  P24  P25  P26  P27  P28  P29   P3   P4   P5   P6   P7   P8   T1 
   1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    2 
 T10  T11 T113 T115 T116 T117 T118  T12 T120 T122 T123 T124 T125 T134 T136 T140 T141 T145 T147  T17 T187   T2  T20  T21  T22 
   1    1    1    1    2    1    1    2    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    2    1 
 T25  T28  T29  T30  T35  T38  T40  T43  T44  T47   T5  T54  T55  T63  T66   T7  T70  T71  T72  T80  T85  T90  T91  T96  T98 
   2    1    1    1    1    2    1    1    1    1    1    1    2    2    1    1    1    1    1    2    1    1    1    1    1 

sample_type: 
tumor 
   84 

er: 
negative positive 
      53       31 

pgr: 
negative positive 
      53       31 

her2: 
negative positive 
      66       18 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  29.00   45.00   53.00   52.89   59.00   85.00 

grade: 
 1  2  3 
10 16 58 

treatment: 
chemotherapy 
          84 

batch: 
DFHCC2_CISPLATIN DFHCC2_REFERENCE 
              24               60 

uncurated_author_metadata: 
         samplename: DFHCC2_P1///dataset: DFHCC2///series: DFHCC2_CISPLATIN///id: P1///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 59///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
       samplename: DFHCC2_P10///dataset: DFHCC2///series: DFHCC2_CISPLATIN///id: P10///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 43///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
       samplename: DFHCC2_P11///dataset: DFHCC2///series: DFHCC2_CISPLATIN///id: P11///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 41///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
       samplename: DFHCC2_P12///dataset: DFHCC2///series: DFHCC2_CISPLATIN///id: P12///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 53///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
       samplename: DFHCC2_P14///dataset: DFHCC2///series: DFHCC2_CISPLATIN///id: P14///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 43///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
       samplename: DFHCC2_P15///dataset: DFHCC2///series: DFHCC2_CISPLATIN///id: P15///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 57///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
       samplename: DFHCC2_P16///dataset: DFHCC2///series: DFHCC2_CISPLATIN///id: P16///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 45///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
       samplename: DFHCC2_P17///dataset: DFHCC2///series: DFHCC2_CISPLATIN///id: P17///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 52///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
         samplename: DFHCC2_P2///dataset: DFHCC2///series: DFHCC2_CISPLATIN///id: P2///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 49///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
       samplename: DFHCC2_P21///dataset: DFHCC2///series: DFHCC2_CISPLATIN///id: P21///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 59///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
       samplename: DFHCC2_P22///dataset: DFHCC2///series: DFHCC2_CISPLATIN///id: P22///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 67///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
       samplename: DFHCC2_P23///dataset: DFHCC2///series: DFHCC2_CISPLATIN///id: P23///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 29///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
       samplename: DFHCC2_P24///dataset: DFHCC2///series: DFHCC2_CISPLATIN///id: P24///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 50///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
       samplename: DFHCC2_P25///dataset: DFHCC2///series: DFHCC2_CISPLATIN///id: P25///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 40///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
       samplename: DFHCC2_P26///dataset: DFHCC2///series: DFHCC2_CISPLATIN///id: P26///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 39///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
       samplename: DFHCC2_P27///dataset: DFHCC2///series: DFHCC2_CISPLATIN///id: P27///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 63///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
       samplename: DFHCC2_P28///dataset: DFHCC2///series: DFHCC2_CISPLATIN///id: P28///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 60///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
       samplename: DFHCC2_P29///dataset: DFHCC2///series: DFHCC2_CISPLATIN///id: P29///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 44///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
         samplename: DFHCC2_P3///dataset: DFHCC2///series: DFHCC2_CISPLATIN///id: P3///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 39///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
         samplename: DFHCC2_P4///dataset: DFHCC2///series: DFHCC2_CISPLATIN///id: P4///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 68///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
         samplename: DFHCC2_P5///dataset: DFHCC2///series: DFHCC2_CISPLATIN///id: P5///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 44///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
         samplename: DFHCC2_P6///dataset: DFHCC2///series: DFHCC2_CISPLATIN///id: P6///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 62///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
         samplename: DFHCC2_P7///dataset: DFHCC2///series: DFHCC2_CISPLATIN///id: P7///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 39///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
         samplename: DFHCC2_P8///dataset: DFHCC2///series: DFHCC2_CISPLATIN///id: P8///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 51///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
       samplename: DFHCC2_REF1///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T1///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 62///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
     samplename: DFHCC2_REF10///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T10///er: 1///pgr: 1///her2: 1///size: NA///node: NA///age: 54///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
     samplename: DFHCC2_REF11///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T11///er: 0///pgr: 0///her2: 1///size: NA///node: NA///age: 49///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
   samplename: DFHCC2_REF113///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T113///er: 0///pgr: 0///her2: 1///size: NA///node: NA///age: 53///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
   samplename: DFHCC2_REF115///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T115///er: 0///pgr: 0///her2: 1///size: NA///node: NA///age: 37///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
   samplename: DFHCC2_REF116///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T116///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 65///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
samplename: DFHCC2_REF116rep///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T116///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 65///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
   samplename: DFHCC2_REF117///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T117///er: 0///pgr: 0///her2: 1///size: NA///node: NA///age: 48///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
   samplename: DFHCC2_REF118///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T118///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 58///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
     samplename: DFHCC2_REF12///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T12///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 84///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
   samplename: DFHCC2_REF120///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T120///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 49///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
   samplename: DFHCC2_REF122///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T122///er: 0///pgr: 0///her2: 1///size: NA///node: NA///age: 85///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
   samplename: DFHCC2_REF123///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T123///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 42///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
   samplename: DFHCC2_REF124///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T124///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 74///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
   samplename: DFHCC2_REF125///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T125///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 66///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
  samplename: DFHCC2_REF12rep///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T12///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 84///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
   samplename: DFHCC2_REF134///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T134///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 40///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
   samplename: DFHCC2_REF136///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T136///er: 0///pgr: 0///her2: 1///size: NA///node: NA///age: 63///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
   samplename: DFHCC2_REF140///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T140///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 48///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
   samplename: DFHCC2_REF141///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T141///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 78///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
   samplename: DFHCC2_REF145///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T145///er: 0///pgr: 0///her2: 1///size: NA///node: NA///age: 49///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
   samplename: DFHCC2_REF147///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T147///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 58///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
     samplename: DFHCC2_REF17///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T17///er: 1///pgr: 1///her2: 1///size: NA///node: NA///age: 35///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
   samplename: DFHCC2_REF187///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T187///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 54///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
    samplename: DFHCC2_REF1rep///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T1///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 62///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
       samplename: DFHCC2_REF2///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T2///er: 1///pgr: 1///her2: 1///size: NA///node: NA///age: 43///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
     samplename: DFHCC2_REF20///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T20///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 45///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
     samplename: DFHCC2_REF21///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T21///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 53///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
  samplename: DFHCC2_REF21rep///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T21///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 53///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
     samplename: DFHCC2_REF22///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T22///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 60///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
     samplename: DFHCC2_REF25///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T25///er: 0///pgr: 0///her2: 1///size: NA///node: NA///age: 45///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
  samplename: DFHCC2_REF25rep///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T25///er: 0///pgr: 0///her2: 1///size: NA///node: NA///age: 45///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
     samplename: DFHCC2_REF28///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T28///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 55///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
     samplename: DFHCC2_REF29///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T29///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 41///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
     samplename: DFHCC2_REF30///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T30///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 46///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
     samplename: DFHCC2_REF35///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T35///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 45///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
     samplename: DFHCC2_REF38///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T38///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 48///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
  samplename: DFHCC2_REF38rep///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T38///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 48///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
     samplename: DFHCC2_REF40///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T40///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 61///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
     samplename: DFHCC2_REF43///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T43///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 57///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
     samplename: DFHCC2_REF44///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T44///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 46///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
     samplename: DFHCC2_REF47///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T47///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 58///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
       samplename: DFHCC2_REF5///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T5///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 46///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
     samplename: DFHCC2_REF54///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T54///er: 1///pgr: 1///her2: 1///size: NA///node: NA///age: 37///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
     samplename: DFHCC2_REF55///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T55///er: 0///pgr: 0///her2: 1///size: NA///node: NA///age: 55///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
  samplename: DFHCC2_REF55rep///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T55///er: 0///pgr: 0///her2: 1///size: NA///node: NA///age: 55///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
     samplename: DFHCC2_REF63///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T63///er: 0///pgr: 0///her2: 1///size: NA///node: NA///age: 53///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
  samplename: DFHCC2_REF63rep///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T63///er: 0///pgr: 0///her2: 1///size: NA///node: NA///age: 53///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
     samplename: DFHCC2_REF66///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T66///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 56///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
       samplename: DFHCC2_REF7///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T7///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 46///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
     samplename: DFHCC2_REF70///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T70///er: 1///pgr: 1///her2: 1///size: NA///node: NA///age: 44///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
     samplename: DFHCC2_REF71///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T71///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 66///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
     samplename: DFHCC2_REF72///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T72///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 38///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
     samplename: DFHCC2_REF80///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T80///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 55///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
  samplename: DFHCC2_REF80rep///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T80///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 55///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
     samplename: DFHCC2_REF85///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T85///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 46///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
     samplename: DFHCC2_REF90///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T90///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 50///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
     samplename: DFHCC2_REF91///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T91///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 74///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
     samplename: DFHCC2_REF96///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T96///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 54///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 
     samplename: DFHCC2_REF98///dataset: DFHCC2///series: DFHCC2_REFERENCE///id: T98///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 53///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                            1 

duplicates: 
   Length     Class      Mode 
       84 character character 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18864 }
\keyword{datasets}
