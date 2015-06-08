\name{ PNC_eset }
\alias{ PNC_eset }
\docType{data}
\title{ PNC_eset }
\usage{data( PNC_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=16280042 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse1456 
  PMIDs: 16280042 
  No abstract available.
  notes:
   summary:     
      Expression profiling was able to better predict prognosis compared to histological staging. 
   version:     
      2015-04-27 19:28:16
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
assayData: 42447 features, 92 samples
Platform type: 
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

   4 observations deleted due to missingness 
records   n.max n.start  events  median 0.95LCL 0.95UCL 
   88.0    88.0    88.0    25.0      NA    11.3      NA 

--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
  PNC_P_1  PNC_P_10 PNC_P_101 PNC_P_102 PNC_P_103 PNC_P_105 PNC_P_106 PNC_P_107 PNC_P_108 PNC_P_109  PNC_P_11 PNC_P_112 
        1         1         1         1         1         1         1         1         1         1         1         1 
PNC_P_113 PNC_P_114 PNC_P_115 PNC_P_116 PNC_P_119  PNC_P_12 PNC_P_120 PNC_P_121 PNC_P_122 PNC_P_123 PNC_P_124 PNC_P_125 
        1         1         1         1         1         1         1         1         1         1         1         1 
PNC_P_126  PNC_P_13  PNC_P_15  PNC_P_16  PNC_P_18  PNC_P_19   PNC_P_2  PNC_P_20  PNC_P_23  PNC_P_24  PNC_P_25  PNC_P_27 
        1         1         1         1         1         1         1         1         1         1         1         1 
 PNC_P_28   PNC_P_3  PNC_P_30  PNC_P_31  PNC_P_32  PNC_P_34  PNC_P_35  PNC_P_37  PNC_P_38  PNC_P_39   PNC_P_4  PNC_P_40 
        1         1         1         1         1         1         1         1         1         1         1         1 
 PNC_P_42  PNC_P_43  PNC_P_45  PNC_P_48  PNC_P_49  PNC_P_50  PNC_P_51  PNC_P_53  PNC_P_54  PNC_P_55  PNC_P_56  PNC_P_57 
        1         1         1         1         1         1         1         1         1         1         1         1 
 PNC_P_58   PNC_P_6  PNC_P_60  PNC_P_61  PNC_P_62  PNC_P_65  PNC_P_66  PNC_P_67  PNC_P_68  PNC_P_69   PNC_P_7  PNC_P_70 
        1         1         1         1         1         1         1         1         1         1         1         1 
 PNC_P_71  PNC_P_72  PNC_P_73  PNC_P_74  PNC_P_75  PNC_P_78  PNC_P_79  PNC_P_80  PNC_P_84  PNC_P_85  PNC_P_86  PNC_P_88 
        1         1         1         1         1         1         1         1         1         1         1         1 
 PNC_P_90  PNC_P_92  PNC_P_94  PNC_P_95  PNC_P_96  PNC_P_97  PNC_P_98  PNC_P_99 
        1         1         1         1         1         1         1         1 

alt_sample_name: 
  P_1  P_10 P_101 P_102 P_103 P_105 P_106 P_107 P_108 P_109  P_11 P_112 P_113 P_114 P_115 P_116 P_119  P_12 P_120 P_121 P_122 
    1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1 
P_123 P_124 P_125 P_126  P_13  P_15  P_16  P_18  P_19   P_2  P_20  P_23  P_24  P_25  P_27  P_28   P_3  P_30  P_31  P_32  P_34 
    1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1 
 P_35  P_37  P_38  P_39   P_4  P_40  P_42  P_43  P_45  P_48  P_49  P_50  P_51  P_53  P_54  P_55  P_56  P_57  P_58   P_6  P_60 
    1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1 
 P_61  P_62  P_65  P_66  P_67  P_68  P_69   P_7  P_70  P_71  P_72  P_73  P_74  P_75  P_78  P_79  P_80  P_84  P_85  P_86  P_88 
    1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1 
 P_90  P_92  P_94  P_95  P_96  P_97  P_98  P_99 
    1     1     1     1     1     1     1     1 

sample_type: 
tumor 
   92 

er: 
negative positive     NA's 
      43       45        4 

pgr: 
negative positive     NA's 
      43       40        9 

her2: 
negative positive     NA's 
      64       26        2 

tumor_size: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  0.900   1.700   2.500   2.758   3.000  10.000       6 

N: 
   0    1 NA's 
  43   40    9 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  32.16   48.57   53.90   55.97   64.84   82.13       4 

grade: 
   1    2    3 NA's 
  13    5   70    4 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
     29     967    2216    2122    2931    5139       7 

recurrence_status: 
deceased_or_recurrence    living_norecurrence                   NA's 
                    36                     49                      7 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    318    1940    2372    2525    3043    5139       4 

vital_status: 
deceased   living     NA's 
      25       63        4 

batch: 
PNC 
 92 

uncurated_author_metadata: 
      samplename: PNC_P_1///dataset: PNC///series: PNC///id: P_1///er: 0///pgr: 0///her2: 1///size: 2.2///node: 1///age: 74.94///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 339///e.rfs: 1///t.os: 3022///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
     samplename: PNC_P_10///dataset: PNC///series: PNC///id: P_10///er: 0///pgr: 0///her2: 1///size: 5.1///node: 1///age: 50.12///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 774///e.rfs: 0///t.os: 920///e.os: 1///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
samplename: PNC_P_101///dataset: PNC///series: PNC///id: P_101///er: 1///pgr: 1///her2: 0///size: 1.5///node: NA///age: 67.57///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2884///e.rfs: 0///t.os: 2884///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_102///dataset: PNC///series: PNC///id: P_102///er: 0///pgr: 1///her2: 1///size: 5///node: 1///age: 52.94///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 4898///e.rfs: 0///t.os: 4898///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_103///dataset: PNC///series: PNC///id: P_103///er: 0///pgr: 1///her2: 0///size: 3///node: 1///age: 46.83///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 4844///e.rfs: 0///t.os: 4844///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
 samplename: PNC_P_105///dataset: PNC///series: PNC///id: P_105///er: 0///pgr: 0///her2: 0///size: 1.3///node: 1///age: 41.54///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 3190///e.rfs: 0///t.os: 3190///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
samplename: PNC_P_106///dataset: PNC///series: PNC///id: P_106///er: 1///pgr: NA///her2: 0///size: 2.4///node: 0///age: 51.81///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 1956///e.rfs: 0///t.os: 1956///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_107///dataset: PNC///series: PNC///id: P_107///er: 1///pgr: 1///her2: 0///size: NA///node: 0///age: 53.8///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 3354///e.rfs: 1///t.os: 4530///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_108///dataset: PNC///series: PNC///id: P_108///er: 0///pgr: 0///her2: 1///size: 2.6///node: 1///age: 55.22///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 270///e.rfs: 1///t.os: 318///e.os: 1///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_109///dataset: PNC///series: PNC///id: P_109///er: 1///pgr: 1///her2: 0///size: 3///node: 1///age: 67.77///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2365///e.rfs: 1///t.os: 2891///e.os: 1///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_11///dataset: PNC///series: PNC///id: P_11///er: 0///pgr: 1///her2: 0///size: 2.2///node: 0///age: 52.64///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 3037///e.rfs: 0///t.os: 3037///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_112///dataset: PNC///series: PNC///id: P_112///er: 0///pgr: 0///her2: 1///size: 6///node: 0///age: 47.82///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 3730///e.rfs: 1///t.os: 3949///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
    samplename: PNC_P_113///dataset: PNC///series: PNC///id: P_113///er: 0///pgr: 1///her2: 1///size: 3///node: 1///age: 64.83///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1438///e.rfs: 1///t.os: 442///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
  samplename: PNC_P_114///dataset: PNC///series: PNC///id: P_114///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
  samplename: PNC_P_115///dataset: PNC///series: PNC///id: P_115///er: 0///pgr: 0///her2: 0///size: 2.5///node: 0///age: 59.08///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 934///e.rfs: 1///t.os: 1785///e.os: 1///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
samplename: PNC_P_116///dataset: PNC///series: PNC///id: P_116///er: 0///pgr: 0///her2: 0///size: 3.7///node: NA///age: 57.32///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 3570///e.rfs: 0///t.os: 3570///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
 samplename: PNC_P_119///dataset: PNC///series: PNC///id: P_119///er: 0///pgr: 1///her2: 0///size: 3.5///node: 0///age: 62.62///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1551///e.rfs: 1///t.os: 1712///e.os: 1///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_12///dataset: PNC///series: PNC///id: P_12///er: 1///pgr: 1///her2: 0///size: 3.8///node: 1///age: 46.58///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 3219///e.rfs: 0///t.os: 3219///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
    samplename: PNC_P_120///dataset: PNC///series: PNC///id: P_120///er: 1///pgr: 0///her2: 0///size: 3///node: 1///age: 57.7///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 3646///e.rfs: 0///t.os: 3646///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
  samplename: PNC_P_121///dataset: PNC///series: PNC///id: P_121///er: 1///pgr: 1///her2: 1///size: 4.6///node: 1///age: 52.35///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 967///e.rfs: 1///t.os: 2037///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_122///dataset: PNC///series: PNC///id: P_122///er: 0///pgr: 1///her2: 1///size: 3///node: NA///age: 69.25///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 234///e.rfs: 1///t.os: 1296///e.os: 1///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
 samplename: PNC_P_123///dataset: PNC///series: PNC///id: P_123///er: 1///pgr: 0///her2: 0///size: 1.3///node: 0///age: 66.14///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 1369///e.rfs: 1///t.os: 1891///e.os: 1///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
  samplename: PNC_P_124///dataset: PNC///series: PNC///id: P_124///er: 0///pgr: 1///her2: 0///size: 1.7///node: 0///age: 53.7///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2259///e.rfs: 0///t.os: 2259///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
 samplename: PNC_P_125///dataset: PNC///series: PNC///id: P_125///er: 1///pgr: 0///her2: 0///size: 2.3///node: 0///age: 67.41///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2230///e.rfs: 0///t.os: 2230///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
  samplename: PNC_P_126///dataset: PNC///series: PNC///id: P_126///er: 1///pgr: 0///her2: 0///size: 3.2///node: 1///age: 71.68///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 2227///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_13///dataset: PNC///series: PNC///id: P_13///er: 1///pgr: 0///her2: 0///size: 1.3///node: 0///age: 63.54///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 3296///e.rfs: 0///t.os: 3296///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_15///dataset: PNC///series: PNC///id: P_15///er: 1///pgr: 1///her2: 0///size: 1.7///node: 0///age: 45.07///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 3278///e.rfs: 0///t.os: 3278///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_16///dataset: PNC///series: PNC///id: P_16///er: 1///pgr: 1///her2: 0///size: 1.7///node: 1///age: 60.57///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2778///e.rfs: 0///t.os: 2778///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_18///dataset: PNC///series: PNC///id: P_18///er: 0///pgr: 1///her2: 0///size: 2.8///node: 0///age: 59.56///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 3449///e.rfs: 0///t.os: 3449///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
    samplename: PNC_P_19///dataset: PNC///series: PNC///id: P_19///er: 1///pgr: 1///her2: 0///size: 2.6///node: 1///age: 59.1///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1562///e.rfs: 1///t.os: 3402///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
      samplename: PNC_P_2///dataset: PNC///series: PNC///id: P_2///er: 0///pgr: 0///her2: 1///size: 4.8///node: 1///age: 59.79///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 445///e.rfs: 1///t.os: 1237///e.os: 1///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_20///dataset: PNC///series: PNC///id: P_20///er: 1///pgr: 1///her2: 0///size: 2.7///node: 1///age: 68.03///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2216///e.rfs: 0///t.os: 2216///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
    samplename: PNC_P_23///dataset: PNC///series: PNC///id: P_23///er: 0///pgr: 1///her2: 0///size: 1.2///node: 0///age: 74.24///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 934///e.rfs: 1///t.os: 1486///e.os: 1///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_24///dataset: PNC///series: PNC///id: P_24///er: 1///pgr: 1///her2: 0///size: 3.4///node: 0///age: 53.78///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 3227///e.rfs: 0///t.os: 3227///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_25///dataset: PNC///series: PNC///id: P_25///er: 1///pgr: 0///her2: 0///size: 1.8///node: 0///age: 49.04///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2931///e.rfs: 0///t.os: 2931///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
     samplename: PNC_P_27///dataset: PNC///series: PNC///id: P_27///er: 1///pgr: 1///her2: 1///size: 2.1///node: 1///age: 53.66///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 675///e.rfs: 0///t.os: 675///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
    samplename: PNC_P_28///dataset: PNC///series: PNC///id: P_28///er: 0///pgr: 1///her2: 0///size: 2.9///node: 0///age: 49.61///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 186///e.rfs: 1///t.os: 1821///e.os: 1///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
     samplename: PNC_P_3///dataset: PNC///series: PNC///id: P_3///er: 0///pgr: 0///her2: 1///size: 1.8///node: 1///age: 52.58///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 3154///e.rfs: 0///t.os: 3154///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_30///dataset: PNC///series: PNC///id: P_30///er: 1///pgr: 0///her2: 0///size: 2.1///node: 1///age: 51.36///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2088///e.rfs: 0///t.os: 2088///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_31///dataset: PNC///series: PNC///id: P_31///er: 1///pgr: 0///her2: 0///size: 1.2///node: 0///age: 60.34///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2716///e.rfs: 0///t.os: 2716///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_32///dataset: PNC///series: PNC///id: P_32///er: 1///pgr: 0///her2: 0///size: 2.7///node: 0///age: 51.61///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1635///e.rfs: 0///t.os: 1635///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_34///dataset: PNC///series: PNC///id: P_34///er: 0///pgr: 0///her2: 0///size: 1.5///node: 0///age: 37.66///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1219///e.rfs: 1///t.os: 2690///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_35///dataset: PNC///series: PNC///id: P_35///er: 0///pgr: 0///her2: 0///size: 2.5///node: NA///age: 54.16///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 2851///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_37///dataset: PNC///series: PNC///id: P_37///er: 0///pgr: 0///her2: 0///size: 2.8///node: 0///age: 51.95///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2862///e.rfs: 0///t.os: 2862///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
    samplename: PNC_P_38///dataset: PNC///series: PNC///id: P_38///er: 1///pgr: 0///her2: 0///size: 2.5///node: 1///age: 77.17///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 748///e.rfs: 1///t.os: 1376///e.os: 1///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
    samplename: PNC_P_39///dataset: PNC///series: PNC///id: P_39///er: 1///pgr: 0///her2: 1///size: 10///node: 1///age: 65.32///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1234///e.rfs: 1///t.os: 1927///e.os: 1///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
     samplename: PNC_P_4///dataset: PNC///series: PNC///id: P_4///er: 1///pgr: 0///her2: 0///size: 4.5///node: 1///age: 44.06///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 3066///e.rfs: 0///t.os: 3066///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
  samplename: PNC_P_40///dataset: PNC///series: PNC///id: P_40///er: 1///pgr: 1///her2: 0///size: 2.65///node: 0///age: 51.96///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 3062///e.rfs: 0///t.os: 3062///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
    samplename: PNC_P_42///dataset: PNC///series: PNC///id: P_42///er: 0///pgr: 0///her2: 1///size: 3.5///node: 1///age: 47.75///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 2785///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
    samplename: PNC_P_43///dataset: PNC///series: PNC///id: P_43///er: 1///pgr: NA///her2: 0///size: 5.5///node: 1///age: 82.13///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 511///e.rfs: 0///t.os: 555///e.os: 1///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_45///dataset: PNC///series: PNC///id: P_45///er: 0///pgr: 0///her2: 0///size: 0.9///node: 0///age: 54.86///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2613///e.rfs: 0///t.os: 2613///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_48///dataset: PNC///series: PNC///id: P_48///er: 1///pgr: 0///her2: 0///size: 9.5///node: 1///age: 45.42///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1424///e.rfs: 1///t.os: 2639///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_49///dataset: PNC///series: PNC///id: P_49///er: 1///pgr: 0///her2: 0///size: 1.6///node: 0///age: 60.46///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2449///e.rfs: 0///t.os: 2449///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_50///dataset: PNC///series: PNC///id: P_50///er: 1///pgr: 0///her2: 1///size: 2.5///node: 0///age: 70.34///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2321///e.rfs: 0///t.os: 2321///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
     samplename: PNC_P_51///dataset: PNC///series: PNC///id: P_51///er: 0///pgr: 0///her2: 0///size: 2///node: 0///age: 73.69///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1829///e.rfs: 1///t.os: 2391///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
      samplename: PNC_P_53///dataset: PNC///series: PNC///id: P_53///er: 1///pgr: NA///her2: 0///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
      samplename: PNC_P_54///dataset: PNC///series: PNC///id: P_54///er: 1///pgr: NA///her2: 0///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
     samplename: PNC_P_55///dataset: PNC///series: PNC///id: P_55///er: 0///pgr: 1///her2: 1///size: 3///node: 0///age: 34.14///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1230///e.rfs: 1///t.os: 1818///e.os: 1///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_56///dataset: PNC///series: PNC///id: P_56///er: 0///pgr: 1///her2: 0///size: 3.9///node: 1///age: 36.17///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2263///e.rfs: 0///t.os: 2263///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
    samplename: PNC_P_57///dataset: PNC///series: PNC///id: P_57///er: 0///pgr: 1///her2: 0///size: 1.8///node: 0///age: 48.7///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2267///e.rfs: 0///t.os: 2267///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
    samplename: PNC_P_58///dataset: PNC///series: PNC///id: P_58///er: 0///pgr: 1///her2: 0///size: 3.5///node: 0///age: 68.48///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 537///e.rfs: 1///t.os: 2029///e.os: 1///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
     samplename: PNC_P_6///dataset: PNC///series: PNC///id: P_6///er: 0///pgr: 1///her2: 1///size: 1.5///node: 0///age: 35.86///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2632///e.rfs: 0///t.os: 2632///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
  samplename: PNC_P_60///dataset: PNC///series: PNC///id: P_60///er: 1///pgr: NA///her2: 0///size: 1.7///node: 1///age: 48.64///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2493///e.rfs: 0///t.os: 2493///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
    samplename: PNC_P_61///dataset: PNC///series: PNC///id: P_61///er: 1///pgr: 0///her2: 0///size: 2.5///node: 0///age: 61.7///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1672///e.rfs: 1///t.os: 2584///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
    samplename: PNC_P_62///dataset: PNC///series: PNC///id: P_62///er: 1///pgr: 1///her2: 1///size: 2.7///node: 0///age: 33.44///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 960///e.rfs: 1///t.os: 2354///e.os: 1///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
 samplename: PNC_P_65///dataset: PNC///series: PNC///id: P_65///er: 0///pgr: 0///her2: 0///size: 1.25///node: NA///age: 54.01///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2270///e.rfs: 0///t.os: 2270///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
     samplename: PNC_P_66///dataset: PNC///series: PNC///id: P_66///er: 1///pgr: 0///her2: 1///size: 2.1///node: 0///age: 49.18///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 204///e.rfs: 0///t.os: 719///e.os: 1///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_67///dataset: PNC///series: PNC///id: P_67///er: 0///pgr: 0///her2: 0///size: 1.4///node: 1///age: 32.16///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 5139///e.rfs: 0///t.os: 5139///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_68///dataset: PNC///series: PNC///id: P_68///er: 1///pgr: 1///her2: 1///size: 1.6///node: 0///age: 66.72///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2252///e.rfs: 0///t.os: 2252///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
     samplename: PNC_P_69///dataset: PNC///series: PNC///id: P_69///er: 1///pgr: 1///her2: 0///size: 5.5///node: 1///age: 78.29///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 29///e.rfs: 1///t.os: 1084///e.os: 1///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
     samplename: PNC_P_7///dataset: PNC///series: PNC///id: P_7///er: 0///pgr: 1///her2: 0///size: 2.2///node: 0///age: 65.42///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2559///e.rfs: 1///t.os: 2989///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
    samplename: PNC_P_70///dataset: PNC///series: PNC///id: P_70///er: 1///pgr: 0///her2: 0///size: NA///node: 0///age: 39.78///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2723///e.rfs: 1///t.os: 4672///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_71///dataset: PNC///series: PNC///id: P_71///er: 1///pgr: 1///her2: 0///size: 1.6///node: 0///age: 42.07///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2172///e.rfs: 0///t.os: 2172///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
    samplename: PNC_P_72///dataset: PNC///series: PNC///id: P_72///er: 1///pgr: 1///her2: 0///size: 2.1///node: 1///age: 53.82///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 617///e.rfs: 1///t.os: 1672///e.os: 1///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_73///dataset: PNC///series: PNC///id: P_73///er: 0///pgr: 0///her2: 0///size: 1.7///node: 1///age: 75.13///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2095///e.rfs: 0///t.os: 2095///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
  samplename: PNC_P_74///dataset: PNC///series: PNC///id: P_74///er: 0///pgr: NA///her2: 0///size: 2.3///node: 1///age: 79.86///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1617///e.rfs: 1///t.os: 1650///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_75///dataset: PNC///series: PNC///id: P_75///er: 1///pgr: 0///her2: 1///size: 5.5///node: 1///age: 48.34///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2234///e.rfs: 0///t.os: 2234///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
     samplename: PNC_P_78///dataset: PNC///series: PNC///id: P_78///er: 0///pgr: 0///her2: 1///size: 1.6///node: 1///age: 50.29///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 394///e.rfs: 1///t.os: 577///e.os: 1///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_79///dataset: PNC///series: PNC///id: P_79///er: 0///pgr: 0///her2: 1///size: 2.2///node: 1///age: 45.98///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1982///e.rfs: 0///t.os: 1982///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
    samplename: PNC_P_80///dataset: PNC///series: PNC///id: P_80///er: 0///pgr: 0///her2: 0///size: 1.1///node: 0///age: 53.98///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 748///e.rfs: 1///t.os: 1945///e.os: 1///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
     samplename: PNC_P_84///dataset: PNC///series: PNC///id: P_84///er: 0///pgr: 1///her2: 0///size: 1.3///node: 0///age: 36.3///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 883///e.rfs: 1///t.os: 2584///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
  samplename: PNC_P_85///dataset: PNC///series: PNC///id: P_85///er: NA///pgr: 0///her2: 1///size: 2.5///node: 1///age: 63.61///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 5095///e.rfs: 0///t.os: 5095///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
     samplename: PNC_P_86///dataset: PNC///series: PNC///id: P_86///er: 0///pgr: 1///her2: 1///size: 3///node: 0///age: 46.95///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 4887///e.rfs: 0///t.os: 4887///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
    samplename: PNC_P_88///dataset: PNC///series: PNC///id: P_88///er: 1///pgr: NA///her2: 0///size: 3///node: 1///age: 64.86///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2161///e.rfs: 1///t.os: 2894///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_90///dataset: PNC///series: PNC///id: P_90///er: 0///pgr: 0///her2: 0///size: 2.4///node: 0///age: 54.19///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1394///e.rfs: 0///t.os: 2708///e.os: 1///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
     samplename: PNC_P_92///dataset: PNC///series: PNC///id: P_92///er: NA///pgr: 1///her2: 1///size: 3///node: 1///age: 76.43///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 763///e.rfs: 1///t.os: 2150///e.os: 1///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_94///dataset: PNC///series: PNC///id: P_94///er: 0///pgr: 1///her2: 1///size: 3.5///node: 1///age: 43.96///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 3165///e.rfs: 1///t.os: 4132///e.os: 1///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
     samplename: PNC_P_95///dataset: PNC///series: PNC///id: P_95///er: 0///pgr: 1///her2: 0///size: 2.5///node: 1///age: 56.55///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 540///e.rfs: 1///t.os: 701///e.os: 1///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_96///dataset: PNC///series: PNC///id: P_96///er: 1///pgr: 1///her2: 0///size: 1.5///node: 0///age: 47.88///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2329///e.rfs: 0///t.os: 2329///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_97///dataset: PNC///series: PNC///id: P_97///er: 1///pgr: 1///her2: 0///size: 1.6///node: 0///age: 52.56///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2175///e.rfs: 0///t.os: 2175///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
    samplename: PNC_P_98///dataset: PNC///series: PNC///id: P_98///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: PNC_P_99///dataset: PNC///series: PNC///id: P_99///er: 1///pgr: 0///her2: 0///size: 1.5///node: 0///age: 59.36///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 4909///e.rfs: 0///t.os: 4909///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse1456 }
\keyword{datasets}
