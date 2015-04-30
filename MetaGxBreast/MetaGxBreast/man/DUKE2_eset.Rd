\name{ DUKE2_eset }
\alias{ DUKE2_eset }
\docType{data}
\title{  }
\description{ Predicting response with gene signature }
\usage{data( DUKE2_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/18024211 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse6861 
  PMIDs: 18024211 

  Abstract: A 5 word abstract is available. Use 'abstract' method.
  notes:
   summary:     
      Retraction in Lancet Feb 2011 (21277543); Regimen specific signatures were able to predict pathological complete response
. Selecting patients with these gene signataures could increase the proportion of patients with pCR than by basing clinic
al decisions on clinical factors. 
   version:     
      2015-04-27 18:55:17
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: 1053_3p_at 117_3p_at ... X79510cds_3p_s_at (45490 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 45490 features, 160 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
 DUKE2_HB01  DUKE2_HB02  DUKE2_HB03  DUKE2_HB04  DUKE2_HB05  DUKE2_HB06  DUKE2_HB07  DUKE2_HB08  DUKE2_HB09  DUKE2_HB10 
          1           1           1           1           1           1           1           1           1           1 
DUKE2_HB100 DUKE2_HB101 DUKE2_HB102 DUKE2_HB103 DUKE2_HB104 DUKE2_HB105 DUKE2_HB106 DUKE2_HB107 DUKE2_HB108 DUKE2_HB109 
          1           1           1           1           1           1           1           1           1           1 
 DUKE2_HB11 DUKE2_HB110 DUKE2_HB111 DUKE2_HB112 DUKE2_HB113 DUKE2_HB114 DUKE2_HB115 DUKE2_HB116 DUKE2_HB117 DUKE2_HB118 
          1           1           1           1           1           1           1           1           1           1 
DUKE2_HB119  DUKE2_HB12 DUKE2_HB120 DUKE2_HB121 DUKE2_HB122 DUKE2_HB123 DUKE2_HB124 DUKE2_HB125 DUKE2_HB126 DUKE2_HB127 
          1           1           1           1           1           1           1           1           1           1 
DUKE2_HB128 DUKE2_HB129  DUKE2_HB13 DUKE2_HB132 DUKE2_HB133 DUKE2_HB134 DUKE2_HB135 DUKE2_HB136 DUKE2_HB137 DUKE2_HB138 
          1           1           1           1           1           1           1           1           1           1 
DUKE2_HB139  DUKE2_HB14 DUKE2_HB140 DUKE2_HB141 DUKE2_HB142 DUKE2_HB143 DUKE2_HB144 DUKE2_HB145 DUKE2_HB146 DUKE2_HB147 
          1           1           1           1           1           1           1           1           1           1 
DUKE2_HB148 DUKE2_HB149  DUKE2_HB15 DUKE2_HB150 DUKE2_HB151 DUKE2_HB152 DUKE2_HB153 DUKE2_HB154 DUKE2_HB155 DUKE2_HB156 
          1           1           1           1           1           1           1           1           1           1 
DUKE2_HB157 DUKE2_HB158 DUKE2_HB159  DUKE2_HB16 DUKE2_HB160 DUKE2_HB161 DUKE2_HB162  DUKE2_HB17  DUKE2_HB18  DUKE2_HB19 
          1           1           1           1           1           1           1           1           1           1 
 DUKE2_HB20  DUKE2_HB21  DUKE2_HB22  DUKE2_HB23  DUKE2_HB24  DUKE2_HB25  DUKE2_HB26  DUKE2_HB27  DUKE2_HB28  DUKE2_HB29 
          1           1           1           1           1           1           1           1           1           1 
 DUKE2_HB30  DUKE2_HB31  DUKE2_HB32  DUKE2_HB33  DUKE2_HB34  DUKE2_HB35  DUKE2_HB36  DUKE2_HB37  DUKE2_HB38     (Other) 
          1           1           1           1           1           1           1           1           1          61 

alt_sample_name: 
     HB01      HB02      HB03      HB04      HB05      HB06      HB07      HB08 HB09septa      HB10     HB100     HB101 
        1         1         1         1         1         1         1         1         1         1         1         1 
   HB102A     HB103     HB104     HB105     HB106     HB107     HB108     HB109      HB11     HB110     HB111     HB112 
        1         1         1         1         1         1         1         1         1         1         1         1 
    HB113     HB114     HB115     HB116     HB117     HB118     HB119      HB12     HB120     HB121     HB122     HB123 
        1         1         1         1         1         1         1         1         1         1         1         1 
    HB124     HB125     HB126     HB127     HB128     HB129      HB13     HB132     HB133     HB134     HB135     HB136 
        1         1         1         1         1         1         1         1         1         1         1         1 
    HB137     HB138     HB139      HB14     HB140     HB141     HB142     HB143     HB144     HB145     HB146     HB147 
        1         1         1         1         1         1         1         1         1         1         1         1 
    HB148     HB149      HB15     HB150     HB151     HB152     HB153     HB154     HB155     HB156     HB157     HB158 
        1         1         1         1         1         1         1         1         1         1         1         1 
    HB159      HB16     HB160     HB161     HB162     HB17R      HB18      HB19      HB20      HB21      HB22      HB23 
        1         1         1         1         1         1         1         1         1         1         1         1 
     HB24      HB25      HB26      HB27      HB28      HB29      HB30      HB31      HB32      HB33      HB34      HB35 
        1         1         1         1         1         1         1         1         1         1         1         1 
     HB36      HB37      HB38   (Other) 
        1         1         1        61 

sample_type: 
tumor 
  160 

er: 
negative positive 
     123       37 

pgr: 
negative positive     NA's 
     133       25        2 

N: 
   0    1 NA's 
  58   95    7 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  26.00   43.00   49.00   49.41   56.00   70.00      35 

grade: 
   1    2    3 NA's 
   2   37   70   51 

treatment: 
chemotherapy 
         160 

batch: 
DUKE2 
  160 

uncurated_author_metadata: 
    samplename: DUKE2_HB01///dataset: DUKE2///series: DUKE2///id: HB01///er: 1///pgr: NA///her2: NA///size: NA///node: 1///age: 54///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
     samplename: DUKE2_HB02///dataset: DUKE2///series: DUKE2///id: HB02///er: 0///pgr: 0///her2: NA///size: NA///node: 0///age: 47///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
     samplename: DUKE2_HB03///dataset: DUKE2///series: DUKE2///id: HB03///er: 1///pgr: 1///her2: NA///size: NA///node: 0///age: 54///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
     samplename: DUKE2_HB04///dataset: DUKE2///series: DUKE2///id: HB04///er: 1///pgr: 1///her2: NA///size: NA///node: 1///age: 55///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
     samplename: DUKE2_HB05///dataset: DUKE2///series: DUKE2///id: HB05///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: 55///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
     samplename: DUKE2_HB06///dataset: DUKE2///series: DUKE2///id: HB06///er: 1///pgr: 1///her2: NA///size: NA///node: 1///age: 61///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
    samplename: DUKE2_HB07///dataset: DUKE2///series: DUKE2///id: HB07///er: 0///pgr: 0///her2: NA///size: NA///node: 0///age: 48///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
     samplename: DUKE2_HB08///dataset: DUKE2///series: DUKE2///id: HB08///er: 1///pgr: 0///her2: NA///size: NA///node: 1///age: 55///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
samplename: DUKE2_HB09///dataset: DUKE2///series: DUKE2///id: HB09septa///er: 1///pgr: 0///her2: NA///size: NA///node: 0///age: 47///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
     samplename: DUKE2_HB10///dataset: DUKE2///series: DUKE2///id: HB10///er: 0///pgr: 1///her2: NA///size: NA///node: 0///age: 54///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB100///dataset: DUKE2///series: DUKE2///id: HB100///er: 0///pgr: 0///her2: NA///size: NA///node: NA///age: 41///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
   samplename: DUKE2_HB101///dataset: DUKE2///series: DUKE2///id: HB101///er: 0///pgr: 0///her2: NA///size: NA///node: 0///age: 46///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB102///dataset: DUKE2///series: DUKE2///id: HB102A///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: 36///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
   samplename: DUKE2_HB103///dataset: DUKE2///series: DUKE2///id: HB103///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: 52///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
   samplename: DUKE2_HB104///dataset: DUKE2///series: DUKE2///id: HB104///er: 0///pgr: 0///her2: NA///size: NA///node: 0///age: 52///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB105///dataset: DUKE2///series: DUKE2///id: HB105///er: 0///pgr: 0///her2: NA///size: NA///node: NA///age: 56///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
   samplename: DUKE2_HB106///dataset: DUKE2///series: DUKE2///id: HB106///er: 0///pgr: 0///her2: NA///size: NA///node: 0///age: 55///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
   samplename: DUKE2_HB107///dataset: DUKE2///series: DUKE2///id: HB107///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: 38///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
   samplename: DUKE2_HB108///dataset: DUKE2///series: DUKE2///id: HB108///er: 0///pgr: 0///her2: NA///size: NA///node: 0///age: 54///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
   samplename: DUKE2_HB109///dataset: DUKE2///series: DUKE2///id: HB109///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: 35///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
     samplename: DUKE2_HB11///dataset: DUKE2///series: DUKE2///id: HB11///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: 68///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
   samplename: DUKE2_HB110///dataset: DUKE2///series: DUKE2///id: HB110///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: 35///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
   samplename: DUKE2_HB111///dataset: DUKE2///series: DUKE2///id: HB111///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: 35///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
   samplename: DUKE2_HB112///dataset: DUKE2///series: DUKE2///id: HB112///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: 70///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
   samplename: DUKE2_HB113///dataset: DUKE2///series: DUKE2///id: HB113///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: 54///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
   samplename: DUKE2_HB114///dataset: DUKE2///series: DUKE2///id: HB114///er: 0///pgr: 0///her2: NA///size: NA///node: 0///age: 45///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
   samplename: DUKE2_HB115///dataset: DUKE2///series: DUKE2///id: HB115///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: 54///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
   samplename: DUKE2_HB116///dataset: DUKE2///series: DUKE2///id: HB116///er: 0///pgr: 0///her2: NA///size: NA///node: 0///age: 39///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
   samplename: DUKE2_HB117///dataset: DUKE2///series: DUKE2///id: HB117///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: 67///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
   samplename: DUKE2_HB118///dataset: DUKE2///series: DUKE2///id: HB118///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: 43///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB119///dataset: DUKE2///series: DUKE2///id: HB119///er: 0///pgr: 0///her2: NA///size: NA///node: 0///age: 57///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
     samplename: DUKE2_HB12///dataset: DUKE2///series: DUKE2///id: HB12///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: 65///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
   samplename: DUKE2_HB120///dataset: DUKE2///series: DUKE2///id: HB120///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: 42///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
   samplename: DUKE2_HB121///dataset: DUKE2///series: DUKE2///id: HB121///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: 45///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
   samplename: DUKE2_HB122///dataset: DUKE2///series: DUKE2///id: HB122///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: 69///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
   samplename: DUKE2_HB123///dataset: DUKE2///series: DUKE2///id: HB123///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: 36///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
 samplename: DUKE2_HB124///dataset: DUKE2///series: DUKE2///id: HB124///er: 0///pgr: 0///her2: NA///size: NA///node: NA///age: 44///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
   samplename: DUKE2_HB125///dataset: DUKE2///series: DUKE2///id: HB125///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: 68///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
 samplename: DUKE2_HB126///dataset: DUKE2///series: DUKE2///id: HB126///er: 0///pgr: 0///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB127///dataset: DUKE2///series: DUKE2///id: HB127///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB128///dataset: DUKE2///series: DUKE2///id: HB128///er: 0///pgr: 0///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB129///dataset: DUKE2///series: DUKE2///id: HB129///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
     samplename: DUKE2_HB13///dataset: DUKE2///series: DUKE2///id: HB13///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: 56///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB132///dataset: DUKE2///series: DUKE2///id: HB132///er: 0///pgr: 0///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB133///dataset: DUKE2///series: DUKE2///id: HB133///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB134///dataset: DUKE2///series: DUKE2///id: HB134///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB135///dataset: DUKE2///series: DUKE2///id: HB135///er: 0///pgr: 0///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB136///dataset: DUKE2///series: DUKE2///id: HB136///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB137///dataset: DUKE2///series: DUKE2///id: HB137///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB138///dataset: DUKE2///series: DUKE2///id: HB138///er: 0///pgr: 0///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB139///dataset: DUKE2///series: DUKE2///id: HB139///er: 0///pgr: 0///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
    samplename: DUKE2_HB14///dataset: DUKE2///series: DUKE2///id: HB14///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: 50///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB140///dataset: DUKE2///series: DUKE2///id: HB140///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB141///dataset: DUKE2///series: DUKE2///id: HB141///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB142///dataset: DUKE2///series: DUKE2///id: HB142///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB143///dataset: DUKE2///series: DUKE2///id: HB143///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB144///dataset: DUKE2///series: DUKE2///id: HB144///er: 0///pgr: 0///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB145///dataset: DUKE2///series: DUKE2///id: HB145///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB146///dataset: DUKE2///series: DUKE2///id: HB146///er: 0///pgr: 0///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB147///dataset: DUKE2///series: DUKE2///id: HB147///er: 0///pgr: 0///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB148///dataset: DUKE2///series: DUKE2///id: HB148///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB149///dataset: DUKE2///series: DUKE2///id: HB149///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
    samplename: DUKE2_HB15///dataset: DUKE2///series: DUKE2///id: HB15///er: 0///pgr: 0///her2: NA///size: NA///node: 0///age: 41///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB150///dataset: DUKE2///series: DUKE2///id: HB150///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB151///dataset: DUKE2///series: DUKE2///id: HB151///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB152///dataset: DUKE2///series: DUKE2///id: HB152///er: 0///pgr: 0///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB153///dataset: DUKE2///series: DUKE2///id: HB153///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB154///dataset: DUKE2///series: DUKE2///id: HB154///er: 0///pgr: 0///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
 samplename: DUKE2_HB155///dataset: DUKE2///series: DUKE2///id: HB155///er: 0///pgr: 0///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB156///dataset: DUKE2///series: DUKE2///id: HB156///er: 0///pgr: 0///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB157///dataset: DUKE2///series: DUKE2///id: HB157///er: 0///pgr: 0///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB158///dataset: DUKE2///series: DUKE2///id: HB158///er: 0///pgr: 0///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB159///dataset: DUKE2///series: DUKE2///id: HB159///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
     samplename: DUKE2_HB16///dataset: DUKE2///series: DUKE2///id: HB16///er: 1///pgr: 0///her2: NA///size: NA///node: 1///age: 45///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB160///dataset: DUKE2///series: DUKE2///id: HB160///er: 0///pgr: 0///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB161///dataset: DUKE2///series: DUKE2///id: HB161///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
  samplename: DUKE2_HB162///dataset: DUKE2///series: DUKE2///id: HB162///er: 0///pgr: 0///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
    samplename: DUKE2_HB17///dataset: DUKE2///series: DUKE2///id: HB17R///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: 50///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
     samplename: DUKE2_HB18///dataset: DUKE2///series: DUKE2///id: HB18///er: 1///pgr: 1///her2: NA///size: NA///node: 1///age: 49///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
     samplename: DUKE2_HB19///dataset: DUKE2///series: DUKE2///id: HB19///er: 1///pgr: 1///her2: NA///size: NA///node: 0///age: 48///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
     samplename: DUKE2_HB20///dataset: DUKE2///series: DUKE2///id: HB20///er: 0///pgr: 0///her2: NA///size: NA///node: 0///age: 42///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
     samplename: DUKE2_HB21///dataset: DUKE2///series: DUKE2///id: HB21///er: 1///pgr: 1///her2: NA///size: NA///node: 0///age: 46///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
     samplename: DUKE2_HB22///dataset: DUKE2///series: DUKE2///id: HB22///er: 1///pgr: 1///her2: NA///size: NA///node: 1///age: 39///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
   samplename: DUKE2_HB23///dataset: DUKE2///series: DUKE2///id: HB23///er: 1///pgr: NA///her2: NA///size: NA///node: 0///age: 36///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
     samplename: DUKE2_HB24///dataset: DUKE2///series: DUKE2///id: HB24///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: 49///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
     samplename: DUKE2_HB25///dataset: DUKE2///series: DUKE2///id: HB25///er: 1///pgr: 1///her2: NA///size: NA///node: 1///age: 58///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
     samplename: DUKE2_HB26///dataset: DUKE2///series: DUKE2///id: HB26///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: 46///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
     samplename: DUKE2_HB27///dataset: DUKE2///series: DUKE2///id: HB27///er: 1///pgr: 0///her2: NA///size: NA///node: 0///age: 51///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
     samplename: DUKE2_HB28///dataset: DUKE2///series: DUKE2///id: HB28///er: 1///pgr: 1///her2: NA///size: NA///node: 1///age: 26///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
     samplename: DUKE2_HB29///dataset: DUKE2///series: DUKE2///id: HB29///er: 1///pgr: 0///her2: NA///size: NA///node: 1///age: 55///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
     samplename: DUKE2_HB30///dataset: DUKE2///series: DUKE2///id: HB30///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: 60///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
     samplename: DUKE2_HB31///dataset: DUKE2///series: DUKE2///id: HB31///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: 43///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
    samplename: DUKE2_HB32///dataset: DUKE2///series: DUKE2///id: HB32///er: 0///pgr: 1///her2: NA///size: NA///node: 0///age: 43///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
     samplename: DUKE2_HB33///dataset: DUKE2///series: DUKE2///id: HB33///er: 0///pgr: 0///her2: NA///size: NA///node: 0///age: 56///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
     samplename: DUKE2_HB34///dataset: DUKE2///series: DUKE2///id: HB34///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: 55///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
     samplename: DUKE2_HB35///dataset: DUKE2///series: DUKE2///id: HB35///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: 31///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
     samplename: DUKE2_HB36///dataset: DUKE2///series: DUKE2///id: HB36///er: 1///pgr: 1///her2: NA///size: NA///node: 1///age: 40///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
     samplename: DUKE2_HB37///dataset: DUKE2///series: DUKE2///id: HB37///er: 0///pgr: 0///her2: NA///size: NA///node: 1///age: 64///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
     samplename: DUKE2_HB38///dataset: DUKE2///series: DUKE2///id: HB38///er: 0///pgr: 0///her2: NA///size: NA///node: 0///age: 48///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                               1 
                                                                                                                                                                                                                                         (Other) 
                                                                                                                                                                                                                                              61 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse6861 }
\keyword{datasets}
