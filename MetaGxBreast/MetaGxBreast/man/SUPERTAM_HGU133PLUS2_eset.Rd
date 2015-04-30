\name{ SUPERTAM_HGU133PLUS2_eset }
\alias{ SUPERTAM_HGU133PLUS2_eset }
\docType{data}
\title{  }
\usage{data( SUPERTAM_HGU133PLUS2_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=17545524 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gsE7390 
  PMIDs: 17545524 
  No abstract available.
  notes:
   summary:     
      The 76-gene signature was validated. The results supports the hypothesis that utilizing the gene signature could reduce t
he number of patients who receive unnecessary adjuvant therapy. 
   version:     
      2015-04-27 19:44:50
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: 1007_s_at 1053_at ... AFFX-HUMISGF3A/M97935_MB_at (42369 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 42369 features, 517 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
GSM102439 GSM102440 GSM102441 GSM102468 GSM102482 GSM102487 GSM102499 GSM102502 GSM102506 GSM102523 GSM102525 GSM102526 
        1         1         1         1         1         1         1         1         1         1         1         1 
GSM102530 GSM102534 GSM102536 GSM102544 GSM102545 GSM102563 GSM102564 GSM102566 GSM102567 GSM102569 GSM102571 GSM102573 
        1         1         1         1         1         1         1         1         1         1         1         1 
GSM102576 GSM102580 GSM102583 GSM117578 GSM117580 GSM117581 GSM117584 GSM117585 GSM117588 GSM117639 GSM117659 GSM117660 
        1         1         1         1         1         1         1         1         1         1         1         1 
GSM117661 GSM117675 GSM117685 GSM117686 GSM117687 GSM117692 GSM117695 GSM117719 GSM117730 GSM117745 GSM117758 GSM117759 
        1         1         1         1         1         1         1         1         1         1         1         1 
GSM117761 GSM117766 GSM117777 GSM137895 GSM137896 GSM137898 GSM137900 GSM137903 GSM137924 GSM137925 GSM137936 GSM137937 
        1         1         1         1         1         1         1         1         1         1         1         1 
GSM137939 GSM137943 GSM137944 GSM137950 GSM137953 GSM137963 GSM137978 GSM137984 GSM137987 GSM137988 GSM137989 GSM137991 
        1         1         1         1         1         1         1         1         1         1         1         1 
GSM137994 GSM138011 GSM138014 GSM138027 GSM138028 GSM138031 GSM138035 GSM138041 GSM138049 GSM138051 GSM138054 GSM152564 
        1         1         1         1         1         1         1         1         1         1         1         1 
GSM152569 GSM152570 GSM152571 GSM152576 GSM152585 GSM152590 GSM152605 GSM152615 GSM152616 GSM152618 GSM152620 GSM152621 
        1         1         1         1         1         1         1         1         1         1         1         1 
GSM152627 GSM152645 GSM152661   (Other) 
        1         1         1       418 

alt_sample_name: 
  50102   50105   50106   50107   50108   50109   50110   50111   50112   50115   50116   50117   50118   50119   50120 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
  50122   50124   50125   50126   50128   50130   50131   50132   50134   50135   50136   50137   50138   50143   50146 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
  50148   50149   50150   50151   50153   50154   50155   50156   50158   50159   50160   50162   50163   50166   50168 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
  50169   50172   50173   50174   50176   50177   50178   50179   50181   50182   50183   50184   50188   50190   50192 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
  50193   50200   50204   50210   50211   50213   50218   50219   50221   50223   50225   50228   50230   50233   50235 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
  50236   50237   50238   50239   50243   50246   50247   50248   50251   50252   50254   50257   50261   50262   50263 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
  50264   50265   50266   50267   50269   50270   50271   50273   50282 (Other) 
      1       1       1       1       1       1       1       1       1     418 

sample_type: 
tumor 
  517 

er: 
negative positive     NA's 
      85      325      107 

pgr: 
negative positive     NA's 
     153      252      112 

her2: 
negative positive     NA's 
     166       61      290 

tumor_size: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  1.070   1.775   2.100   2.425   3.000   7.500     353 

N: 
   0    1 NA's 
  70   94  353 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  25.00   53.00   61.00   60.72   70.00   95.00       1 

grade: 
   1    2    3 NA's 
  63  171  191   92 

dmfs_days: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    188    2552    3278    3294    4336    6151     353 

dmfs_status: 
deceased_or_recurrence    living_norecurrence                   NA's 
                    38                    126                    353 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    188    2484    3264    3266    4336    6151     353 

recurrence_status: 
deceased_or_recurrence    living_norecurrence                   NA's 
                    41                    123                    353 

treatment: 
hormonotherapy           NA's 
           164            353 

batch: 
 EXPO  GUYT GUYT2 
  353    87    77 

uncurated_author_metadata: 
 samplename: GSM102439///dataset: EXPO///series: EXPO///id: GSM102439///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: 45///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM102440///dataset: EXPO///series: EXPO///id: GSM102440///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 45///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
   samplename: GSM102441///dataset: EXPO///series: EXPO///id: GSM102441///er: 0///pgr: 0///her2: NA///size: NA///node: NA///age: 55///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM102468///dataset: EXPO///series: EXPO///id: GSM102468///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 35///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
 samplename: GSM102482///dataset: EXPO///series: EXPO///id: GSM102482///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: 65///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
 samplename: GSM102487///dataset: EXPO///series: EXPO///id: GSM102487///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: 75///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM102499///dataset: EXPO///series: EXPO///id: GSM102499///er: 1///pgr: 1///her2: 1///size: NA///node: NA///age: 55///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM102502///dataset: EXPO///series: EXPO///id: GSM102502///er: 1///pgr: 0///her2: 1///size: NA///node: NA///age: 75///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM102506///dataset: EXPO///series: EXPO///id: GSM102506///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 35///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM102523///dataset: EXPO///series: EXPO///id: GSM102523///er: 1///pgr: 0///her2: 0///size: NA///node: NA///age: 65///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM102525///dataset: EXPO///series: EXPO///id: GSM102525///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 55///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
 samplename: GSM102526///dataset: EXPO///series: EXPO///id: GSM102526///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: 45///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM102530///dataset: EXPO///series: EXPO///id: GSM102530///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 85///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM102534///dataset: EXPO///series: EXPO///id: GSM102534///er: 0///pgr: 0///her2: 1///size: NA///node: NA///age: 55///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM102536///dataset: EXPO///series: EXPO///id: GSM102536///er: 1///pgr: 1///her2: 1///size: NA///node: NA///age: 55///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM102544///dataset: EXPO///series: EXPO///id: GSM102544///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 45///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
 samplename: GSM102545///dataset: EXPO///series: EXPO///id: GSM102545///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: 65///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
 samplename: GSM102563///dataset: EXPO///series: EXPO///id: GSM102563///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: 55///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM102564///dataset: EXPO///series: EXPO///id: GSM102564///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 45///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM102566///dataset: EXPO///series: EXPO///id: GSM102566///er: 1///pgr: 0///her2: 0///size: NA///node: NA///age: 75///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
 samplename: GSM102567///dataset: EXPO///series: EXPO///id: GSM102567///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: 35///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM102569///dataset: EXPO///series: EXPO///id: GSM102569///er: 1///pgr: 1///her2: 1///size: NA///node: NA///age: 45///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM102571///dataset: EXPO///series: EXPO///id: GSM102571///er: 1///pgr: 1///her2: 1///size: NA///node: NA///age: 45///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM102573///dataset: EXPO///series: EXPO///id: GSM102573///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 35///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
samplename: GSM102576///dataset: EXPO///series: EXPO///id: GSM102576///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: 65///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM102580///dataset: EXPO///series: EXPO///id: GSM102580///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 45///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM102583///dataset: EXPO///series: EXPO///id: GSM102583///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 65///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM117578///dataset: EXPO///series: EXPO///id: GSM117578///er: 0///pgr: 0///her2: 1///size: NA///node: NA///age: 65///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM117580///dataset: EXPO///series: EXPO///id: GSM117580///er: 1///pgr: 0///her2: 0///size: NA///node: NA///age: 65///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM117581///dataset: EXPO///series: EXPO///id: GSM117581///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 45///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM117584///dataset: EXPO///series: EXPO///id: GSM117584///er: 1///pgr: 0///her2: 0///size: NA///node: NA///age: 65///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM117585///dataset: EXPO///series: EXPO///id: GSM117585///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 35///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM117588///dataset: EXPO///series: EXPO///id: GSM117588///er: 1///pgr: 0///her2: 1///size: NA///node: NA///age: 45///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
 samplename: GSM117639///dataset: EXPO///series: EXPO///id: GSM117639///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: 75///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM117659///dataset: EXPO///series: EXPO///id: GSM117659///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 65///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM117660///dataset: EXPO///series: EXPO///id: GSM117660///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 65///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM117661///dataset: EXPO///series: EXPO///id: GSM117661///er: 0///pgr: 0///her2: 1///size: NA///node: NA///age: 55///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM117675///dataset: EXPO///series: EXPO///id: GSM117675///er: 1///pgr: 1///her2: 1///size: NA///node: NA///age: 85///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM117685///dataset: EXPO///series: EXPO///id: GSM117685///er: 1///pgr: 1///her2: 1///size: NA///node: NA///age: 75///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM117686///dataset: EXPO///series: EXPO///id: GSM117686///er: 1///pgr: 1///her2: 1///size: NA///node: NA///age: 45///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM117687///dataset: EXPO///series: EXPO///id: GSM117687///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 45///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM117692///dataset: EXPO///series: EXPO///id: GSM117692///er: 0///pgr: 0///her2: 1///size: NA///node: NA///age: 65///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM117695///dataset: EXPO///series: EXPO///id: GSM117695///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 75///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
   samplename: GSM117719///dataset: EXPO///series: EXPO///id: GSM117719///er: 1///pgr: 1///her2: NA///size: NA///node: NA///age: 55///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM117730///dataset: EXPO///series: EXPO///id: GSM117730///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 45///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM117745///dataset: EXPO///series: EXPO///id: GSM117745///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 55///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM117758///dataset: EXPO///series: EXPO///id: GSM117758///er: 0///pgr: 0///her2: 1///size: NA///node: NA///age: 45///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM117759///dataset: EXPO///series: EXPO///id: GSM117759///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 55///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM117761///dataset: EXPO///series: EXPO///id: GSM117761///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 45///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM117766///dataset: EXPO///series: EXPO///id: GSM117766///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 45///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM117777///dataset: EXPO///series: EXPO///id: GSM117777///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 55///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
 samplename: GSM137895///dataset: EXPO///series: EXPO///id: GSM137895///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: 35///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM137896///dataset: EXPO///series: EXPO///id: GSM137896///er: 1///pgr: 0///her2: 0///size: NA///node: NA///age: 65///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM137898///dataset: EXPO///series: EXPO///id: GSM137898///er: 1///pgr: 1///her2: 1///size: NA///node: NA///age: 25///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
  samplename: GSM137900///dataset: EXPO///series: EXPO///id: GSM137900///er: 0///pgr: 0///her2: NA///size: NA///node: NA///age: 75///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM137903///dataset: EXPO///series: EXPO///id: GSM137903///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 65///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM137924///dataset: EXPO///series: EXPO///id: GSM137924///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 35///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM137925///dataset: EXPO///series: EXPO///id: GSM137925///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 85///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM137936///dataset: EXPO///series: EXPO///id: GSM137936///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 55///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM137937///dataset: EXPO///series: EXPO///id: GSM137937///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 35///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
   samplename: GSM137939///dataset: EXPO///series: EXPO///id: GSM137939///er: 1///pgr: 0///her2: NA///size: NA///node: NA///age: 85///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM137943///dataset: EXPO///series: EXPO///id: GSM137943///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 85///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM137944///dataset: EXPO///series: EXPO///id: GSM137944///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 55///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM137950///dataset: EXPO///series: EXPO///id: GSM137950///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 45///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM137953///dataset: EXPO///series: EXPO///id: GSM137953///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 65///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM137963///dataset: EXPO///series: EXPO///id: GSM137963///er: 1///pgr: 0///her2: 0///size: NA///node: NA///age: 65///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM137978///dataset: EXPO///series: EXPO///id: GSM137978///er: 1///pgr: 0///her2: 0///size: NA///node: NA///age: 45///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
   samplename: GSM137984///dataset: EXPO///series: EXPO///id: GSM137984///er: 0///pgr: NA///her2: 0///size: NA///node: NA///age: 55///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM137987///dataset: EXPO///series: EXPO///id: GSM137987///er: 1///pgr: 1///her2: 1///size: NA///node: NA///age: 55///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
   samplename: GSM137988///dataset: EXPO///series: EXPO///id: GSM137988///er: 0///pgr: 0///her2: NA///size: NA///node: NA///age: 75///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM137989///dataset: EXPO///series: EXPO///id: GSM137989///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 55///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM137991///dataset: EXPO///series: EXPO///id: GSM137991///er: 0///pgr: 0///her2: 1///size: NA///node: NA///age: 75///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM137994///dataset: EXPO///series: EXPO///id: GSM137994///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 95///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
 samplename: GSM138011///dataset: EXPO///series: EXPO///id: GSM138011///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: 85///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM138014///dataset: EXPO///series: EXPO///id: GSM138014///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 65///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
 samplename: GSM138027///dataset: EXPO///series: EXPO///id: GSM138027///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: 95///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM138028///dataset: EXPO///series: EXPO///id: GSM138028///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 45///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM138031///dataset: EXPO///series: EXPO///id: GSM138031///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 45///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM138035///dataset: EXPO///series: EXPO///id: GSM138035///er: 0///pgr: 0///her2: 1///size: NA///node: NA///age: 75///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM138041///dataset: EXPO///series: EXPO///id: GSM138041///er: 0///pgr: 0///her2: 1///size: NA///node: NA///age: 75///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
   samplename: GSM138049///dataset: EXPO///series: EXPO///id: GSM138049///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 45///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM138051///dataset: EXPO///series: EXPO///id: GSM138051///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 45///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM138054///dataset: EXPO///series: EXPO///id: GSM138054///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 55///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
  samplename: GSM152564///dataset: EXPO///series: EXPO///id: GSM152564///er: 1///pgr: 0///her2: NA///size: NA///node: NA///age: 75///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
samplename: GSM152569///dataset: EXPO///series: EXPO///id: GSM152569///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: 75///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
 samplename: GSM152570///dataset: EXPO///series: EXPO///id: GSM152570///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: 55///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM152571///dataset: EXPO///series: EXPO///id: GSM152571///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 55///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
 samplename: GSM152576///dataset: EXPO///series: EXPO///id: GSM152576///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: 45///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM152585///dataset: EXPO///series: EXPO///id: GSM152585///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 75///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
   samplename: GSM152590///dataset: EXPO///series: EXPO///id: GSM152590///er: 1///pgr: 0///her2: NA///size: NA///node: NA///age: 85///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
 samplename: GSM152605///dataset: EXPO///series: EXPO///id: GSM152605///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: 65///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM152615///dataset: EXPO///series: EXPO///id: GSM152615///er: 1///pgr: 1///her2: 1///size: NA///node: NA///age: 45///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
 samplename: GSM152616///dataset: EXPO///series: EXPO///id: GSM152616///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: 55///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM152618///dataset: EXPO///series: EXPO///id: GSM152618///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 65///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
 samplename: GSM152620///dataset: EXPO///series: EXPO///id: GSM152620///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: 65///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM152621///dataset: EXPO///series: EXPO///id: GSM152621///er: 1///pgr: 0///her2: 0///size: NA///node: NA///age: 55///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM152627///dataset: EXPO///series: EXPO///id: GSM152627///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 55///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
    samplename: GSM152645///dataset: EXPO///series: EXPO///id: GSM152645///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 55///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
 samplename: GSM152661///dataset: EXPO///series: EXPO///id: GSM152661///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: 65///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                 1 
                                                                                                                                                                                                                                           (Other) 
                                                                                                                                                                                                                                               418 

duplicates: 
   Length     Class      Mode 
      517 character character 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gsE7390 }
\keyword{datasets}
