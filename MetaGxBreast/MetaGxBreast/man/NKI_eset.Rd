\name{ NKI_eset }
\alias{ NKI_eset }
\docType{data}
\title{  }
\usage{data( NKI_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=21910250 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20711 
  PMIDs: GSE20711, PMID 21910250 
  No abstract available.
  notes:
   summary:     
      Breast tumors can be further divided than the currently known expression subtypes based on DNA methylation profiles.
   version:     
      2015-04-27 19:24:07
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: Contig45645_RC Contig44916_RC ... Contig62037_RC (14960 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 14960 features, 337 samples
Platform type: 
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

   42 observations deleted due to missingness 
records   n.max n.start  events  median 0.95LCL 0.95UCL 
    295     295     295      79      NA      NA      NA 

--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
NKI_100 NKI_102 NKI_103 NKI_104 NKI_106 NKI_107 NKI_108 NKI_109  NKI_11 NKI_110 NKI_111 NKI_113 NKI_114 NKI_116 NKI_117 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
NKI_118 NKI_119  NKI_12 NKI_120 NKI_122 NKI_123 NKI_124 NKI_125 NKI_126 NKI_127 NKI_128 NKI_129  NKI_13 NKI_130 NKI_131 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
NKI_132 NKI_133 NKI_134 NKI_135 NKI_136 NKI_137 NKI_138 NKI_139  NKI_14 NKI_140 NKI_141 NKI_142 NKI_144 NKI_145 NKI_146 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
NKI_147 NKI_148 NKI_149 NKI_150 NKI_151 NKI_153 NKI_154 NKI_155 NKI_156 NKI_157 NKI_158 NKI_159 NKI_160 NKI_161 NKI_162 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
NKI_163 NKI_164 NKI_165 NKI_166 NKI_167 NKI_169  NKI_17 NKI_170 NKI_172 NKI_174 NKI_175 NKI_176 NKI_177 NKI_178 NKI_179 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
NKI_180 NKI_181 NKI_182 NKI_183 NKI_184 NKI_185 NKI_186 NKI_187 NKI_188 NKI_189 NKI_190 NKI_191 NKI_192 NKI_193 NKI_194 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
NKI_195 NKI_196 NKI_197 NKI_198 NKI_199 NKI_200 NKI_201 NKI_202 NKI_203 (Other) 
      1       1       1       1       1       1       1       1       1     238 

alt_sample_name: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    4.0   123.0   215.0   214.1   312.0   404.0 

sample_type: 
tumor 
  337 

er: 
negative positive 
      88      249 

tumor_size: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  0.200   1.500   2.000   2.241   2.800   5.500 

N: 
  0   1 
193 144 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   26.0    40.0    45.0    44.2    49.0    62.0 

grade: 
  1   2   3 
 79 109 149 

dmfs_days: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
      9    1252    2414    2546    3602    6699      18 

dmfs_status: 
deceased_or_recurrence    living_norecurrence                   NA's 
                   109                    210                     18 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
      9    1252    2414    2546    3602    6699      18 

recurrence_status: 
deceased_or_recurrence    living_norecurrence                   NA's 
                   109                    210                     18 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
     20    1934    2637    2870    3763    6694      42 

vital_status: 
deceased   living     NA's 
      79      216       42 

treatment: 
  chemotherapy hormonotherapy      untreated 
            90             40            207 

batch: 
 NKI NKI2 
 117  220 

uncurated_author_metadata: 
      samplename: NKI_100///dataset: NKI///series: NKI///id: 100///er: 0///pgr: NA///her2: NA///size: 2///node: 0///age: 48///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
  samplename: NKI_102///dataset: NKI///series: NKI///id: 102///er: 1///pgr: NA///her2: NA///size: 1.5///node: 0///age: 50///grade: 3///t.dmfs: 1201///e.dmfs: 1///t.rfs: 1201///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
   samplename: NKI_103///dataset: NKI///series: NKI///id: 103///er: 0///pgr: NA///her2: NA///size: 3///node: 0///age: 42///grade: 3///t.dmfs: 1809///e.dmfs: 1///t.rfs: 1809///e.rfs: 1///t.os: 2110///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
    samplename: NKI_104///dataset: NKI///series: NKI///id: 104///er: 0///pgr: NA///her2: NA///size: 0.5///node: 0///age: 34///grade: 2///t.dmfs: 672///e.dmfs: 1///t.rfs: 672///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
    samplename: NKI_106///dataset: NKI///series: NKI///id: 106///er: 0///pgr: NA///her2: NA///size: 1.5///node: 0///age: 35///grade: 3///t.dmfs: 701///e.dmfs: 1///t.rfs: 701///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
   samplename: NKI_107///dataset: NKI///series: NKI///id: 107///er: 1///pgr: NA///her2: NA///size: 0.2///node: 0///age: 38///grade: 3///t.dmfs: 929///e.dmfs: 1///t.rfs: 929///e.rfs: 1///t.os: 1261///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
    samplename: NKI_108///dataset: NKI///series: NKI///id: 108///er: 1///pgr: NA///her2: NA///size: 1.7///node: 0///age: 48///grade: 2///t.dmfs: 325///e.dmfs: 1///t.rfs: 325///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
   samplename: NKI_109///dataset: NKI///series: NKI///id: 109///er: 0///pgr: NA///her2: NA///size: 2///node: 0///age: 41///grade: 3///t.dmfs: 1167///e.dmfs: 1///t.rfs: 1167///e.rfs: 1///t.os: 1177///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
   samplename: NKI_11///dataset: NKI///series: NKI///id: 11///er: 1///pgr: NA///her2: NA///size: 2.2///node: 0///age: 37///grade: 3///t.dmfs: 2120///e.dmfs: 0///t.rfs: 2120///e.rfs: 0///t.os: 2119///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
    samplename: NKI_110///dataset: NKI///series: NKI///id: 110///er: 1///pgr: NA///her2: NA///size: 3.5///node: 0///age: 51///grade: 3///t.dmfs: 792///e.dmfs: 1///t.rfs: 792///e.rfs: 1///t.os: 843///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
   samplename: NKI_111///dataset: NKI///series: NKI///id: 111///er: 1///pgr: NA///her2: NA///size: 2.5///node: 0///age: 50///grade: 3///t.dmfs: 464///e.dmfs: 1///t.rfs: 464///e.rfs: 1///t.os: 1187///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
   samplename: NKI_113///dataset: NKI///series: NKI///id: 113///er: 0///pgr: NA///her2: NA///size: 1.2///node: 0///age: 47///grade: 2///t.dmfs: 364///e.dmfs: 1///t.rfs: 364///e.rfs: 1///t.os: 1183///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
  samplename: NKI_114///dataset: NKI///series: NKI///id: 114///er: 1///pgr: NA///her2: NA///size: 1.5///node: 0///age: 54///grade: 1///t.dmfs: 4537///e.dmfs: 0///t.rfs: 4537///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
    samplename: NKI_116///dataset: NKI///series: NKI///id: 116///er: 1///pgr: NA///her2: NA///size: 1///node: 0///age: 54///grade: 1///t.dmfs: 3117///e.dmfs: 0///t.rfs: 3117///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
 samplename: NKI_117///dataset: NKI///series: NKI///id: 117///er: 1///pgr: NA///her2: NA///size: 1.5///node: 1///age: 51///grade: 3///t.dmfs: 1937///e.dmfs: 0///t.rfs: 1937///e.rfs: 0///t.os: 1936///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
 samplename: NKI_118///dataset: NKI///series: NKI///id: 118///er: 1///pgr: NA///her2: NA///size: 1.7///node: 0///age: 47///grade: 2///t.dmfs: 1911///e.dmfs: 0///t.rfs: 1911///e.rfs: 0///t.os: 1910///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
    samplename: NKI_119///dataset: NKI///series: NKI///id: 119///er: 0///pgr: NA///her2: NA///size: 1///node: 0///age: 54///grade: 3///t.dmfs: 3643///e.dmfs: 0///t.rfs: 3643///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
   samplename: NKI_12///dataset: NKI///series: NKI///id: 12///er: 0///pgr: NA///her2: NA///size: 2.2///node: 0///age: 46///grade: 3///t.dmfs: 2870///e.dmfs: 0///t.rfs: 2870///e.rfs: 0///t.os: 2868///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
 samplename: NKI_120///dataset: NKI///series: NKI///id: 120///er: 1///pgr: NA///her2: NA///size: 1.3///node: 0///age: 42///grade: 2///t.dmfs: 3688///e.dmfs: 0///t.rfs: 3688///e.rfs: 0///t.os: 3685///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
 samplename: NKI_122///dataset: NKI///series: NKI///id: 122///er: 1///pgr: NA///her2: NA///size: 2.5///node: 0///age: 43///grade: 2///t.dmfs: 5412///e.dmfs: 0///t.rfs: 5412///e.rfs: 0///t.os: 5408///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
   samplename: NKI_123///dataset: NKI///series: NKI///id: 123///er: 1///pgr: NA///her2: NA///size: 2///node: 0///age: 48///grade: 3///t.dmfs: 5209///e.dmfs: 0///t.rfs: 5209///e.rfs: 0///t.os: 5205///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
 samplename: NKI_124///dataset: NKI///series: NKI///id: 124///er: 1///pgr: NA///her2: NA///size: 1.5///node: 0///age: 38///grade: 2///t.dmfs: 2427///e.dmfs: 0///t.rfs: 2427///e.rfs: 0///t.os: 2425///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
samplename: NKI_125///dataset: NKI///series: NKI2///id: 125///er: 1///pgr: NA///her2: NA///size: 1.5///node: 1///age: 50///grade: 2///t.dmfs: 2830///e.dmfs: 0///t.rfs: 2830///e.rfs: 0///t.os: 2828///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                        1 
samplename: NKI_126///dataset: NKI///series: NKI2///id: 126///er: 1///pgr: NA///her2: NA///size: 1.5///node: 0///age: 38///grade: 2///t.dmfs: 2308///e.dmfs: 1///t.rfs: 2308///e.rfs: 1///t.os: 2349///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
  samplename: NKI_127///dataset: NKI///series: NKI2///id: 127///er: 1///pgr: NA///her2: NA///size: 1///node: 1///age: 42///grade: 1///t.dmfs: 1703///e.dmfs: 1///t.rfs: 1703///e.rfs: 1///t.os: 1839///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                        1 
samplename: NKI_128///dataset: NKI///series: NKI2///id: 128///er: 1///pgr: NA///her2: NA///size: 2.5///node: 1///age: 50///grade: 1///t.dmfs: 3192///e.dmfs: 0///t.rfs: 3192///e.rfs: 0///t.os: 3190///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                        1 
samplename: NKI_129///dataset: NKI///series: NKI2///id: 129///er: 1///pgr: NA///her2: NA///size: 1.5///node: 1///age: 43///grade: 2///t.dmfs: 2764///e.dmfs: 0///t.rfs: 2764///e.rfs: 0///t.os: 2762///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                        1 
   samplename: NKI_13///dataset: NKI///series: NKI///id: 13///er: 1///pgr: NA///her2: NA///size: 2.1///node: 0///age: 41///grade: 3///t.dmfs: 2983///e.dmfs: 0///t.rfs: 2983///e.rfs: 0///t.os: 2981///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
samplename: NKI_130///dataset: NKI///series: NKI2///id: 130///er: 0///pgr: NA///her2: NA///size: 1.8///node: 1///age: 47///grade: 3///t.dmfs: 2665///e.dmfs: 0///t.rfs: 2665///e.rfs: 0///t.os: 2663///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                        1 
  samplename: NKI_131///dataset: NKI///series: NKI2///id: 131///er: 0///pgr: NA///her2: NA///size: 1.7///node: 0///age: 39///grade: 3///t.dmfs: 407///e.dmfs: 0///t.rfs: 407///e.rfs: 0///t.os: 1702///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
samplename: NKI_132///dataset: NKI///series: NKI2///id: 132///er: 1///pgr: NA///her2: NA///size: 1.5///node: 1///age: 47///grade: 2///t.dmfs: 2143///e.dmfs: 0///t.rfs: 2143///e.rfs: 0///t.os: 2452///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                        1 
 samplename: NKI_133///dataset: NKI///series: NKI///id: 133///er: 1///pgr: NA///her2: NA///size: 1.8///node: 0///age: 32///grade: 1///t.dmfs: 3159///e.dmfs: 0///t.rfs: 3159///e.rfs: 0///t.os: 3157///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
samplename: NKI_134///dataset: NKI///series: NKI2///id: 134///er: 1///pgr: NA///her2: NA///size: 1.2///node: 1///age: 38///grade: 3///t.dmfs: 2555///e.dmfs: 1///t.rfs: 2555///e.rfs: 1///t.os: 2589///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                        1 
  samplename: NKI_135///dataset: NKI///series: NKI2///id: 135///er: 0///pgr: NA///her2: NA///size: 2///node: 1///age: 45///grade: 3///t.dmfs: 3408///e.dmfs: 0///t.rfs: 3408///e.rfs: 0///t.os: 3406///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                        1 
  samplename: NKI_136///dataset: NKI///series: NKI2///id: 136///er: 0///pgr: NA///her2: NA///size: 4///node: 1///age: 31///grade: 3///t.dmfs: 1256///e.dmfs: 1///t.rfs: 1256///e.rfs: 1///t.os: 1395///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
samplename: NKI_137///dataset: NKI///series: NKI2///id: 137///er: 1///pgr: NA///her2: NA///size: 1.8///node: 1///age: 41///grade: 1///t.dmfs: 5599///e.dmfs: 0///t.rfs: 5599///e.rfs: 0///t.os: 5595///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                        1 
 samplename: NKI_138///dataset: NKI///series: NKI///id: 138///er: 1///pgr: NA///her2: NA///size: 1.5///node: 0///age: 44///grade: 3///t.dmfs: 1269///e.dmfs: 1///t.rfs: 1269///e.rfs: 1///t.os: 1405///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
 samplename: NKI_139///dataset: NKI///series: NKI///id: 139///er: 1///pgr: NA///her2: NA///size: 4.5///node: 0///age: 41///grade: 3///t.dmfs: 4663///e.dmfs: 0///t.rfs: 4663///e.rfs: 0///t.os: 4660///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
     samplename: NKI_14///dataset: NKI///series: NKI///id: 14///er: 1///pgr: NA///her2: NA///size: 1///node: 0///age: 48///grade: 1///t.dmfs: 3007///e.dmfs: 0///t.rfs: 3007///e.rfs: 0///t.os: 3005///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
samplename: NKI_140///dataset: NKI///series: NKI2///id: 140///er: 1///pgr: NA///her2: NA///size: 1.8///node: 0///age: 46///grade: 1///t.dmfs: 2029///e.dmfs: 0///t.rfs: 2029///e.rfs: 0///t.os: 2028///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
      samplename: NKI_141///dataset: NKI///series: NKI///id: 141///er: 0///pgr: NA///her2: NA///size: 3///node: 0///age: 33///grade: 3///t.dmfs: 512///e.dmfs: 1///t.rfs: 512///e.rfs: 1///t.os: 753///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
  samplename: NKI_142///dataset: NKI///series: NKI2///id: 142///er: 1///pgr: NA///her2: NA///size: 2///node: 0///age: 39///grade: 3///t.dmfs: 5528///e.dmfs: 0///t.rfs: 5528///e.rfs: 0///t.os: 5524///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
  samplename: NKI_144///dataset: NKI///series: NKI2///id: 144///er: 0///pgr: NA///her2: NA///size: 3///node: 0///age: 38///grade: 3///t.dmfs: 5160///e.dmfs: 0///t.rfs: 5160///e.rfs: 0///t.os: 5156///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
samplename: NKI_145///dataset: NKI///series: NKI2///id: 145///er: 1///pgr: NA///her2: NA///size: 1.9///node: 1///age: 48///grade: 2///t.dmfs: 2004///e.dmfs: 0///t.rfs: 2004///e.rfs: 0///t.os: 2003///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                        1 
  samplename: NKI_146///dataset: NKI///series: NKI2///id: 146///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 47///grade: 3///t.dmfs: 1335///e.dmfs: 1///t.rfs: 1335///e.rfs: 1///t.os: 3434///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
   samplename: NKI_147///dataset: NKI///series: NKI2///id: 147///er: 0///pgr: NA///her2: NA///size: 1.7///node: 1///age: 38///grade: 2///t.dmfs: 588///e.dmfs: 1///t.rfs: 588///e.rfs: 1///t.os: 988///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
  samplename: NKI_148///dataset: NKI///series: NKI2///id: 148///er: 1///pgr: NA///her2: NA///size: 5///node: 0///age: 49///grade: 3///t.dmfs: 6699///e.dmfs: 0///t.rfs: 6699///e.rfs: 0///t.os: 6694///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
  samplename: NKI_149///dataset: NKI///series: NKI2///id: 149///er: 0///pgr: NA///her2: NA///size: 3///node: 1///age: 44///grade: 2///t.dmfs: 6297///e.dmfs: 0///t.rfs: 6297///e.rfs: 0///t.os: 6293///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                        1 
     samplename: NKI_150///dataset: NKI///series: NKI2///id: 150///er: 1///pgr: NA///her2: NA///size: 2///node: 1///age: 36///grade: 2///t.dmfs: 351///e.dmfs: 1///t.rfs: 351///e.rfs: 1///t.os: 543///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
  samplename: NKI_151///dataset: NKI///series: NKI2///id: 151///er: 1///pgr: NA///her2: NA///size: 4///node: 1///age: 42///grade: 3///t.dmfs: 5118///e.dmfs: 1///t.rfs: 5118///e.rfs: 1///t.os: 6415///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                        1 
    samplename: NKI_153///dataset: NKI///series: NKI2///id: 153///er: 0///pgr: NA///her2: NA///size: 5///node: 0///age: 37///grade: 3///t.dmfs: 430///e.dmfs: 1///t.rfs: 430///e.rfs: 1///t.os: 1108///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
samplename: NKI_154///dataset: NKI///series: NKI2///id: 154///er: 1///pgr: NA///her2: NA///size: 1.5///node: 0///age: 40///grade: 1///t.dmfs: 5517///e.dmfs: 0///t.rfs: 5517///e.rfs: 0///t.os: 5513///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
     samplename: NKI_155///dataset: NKI///series: NKI2///id: 155///er: 1///pgr: NA///her2: NA///size: 4///node: 0///age: 49///grade: 3///t.dmfs: 340///e.dmfs: 1///t.rfs: 340///e.rfs: 1///t.os: 675///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
samplename: NKI_156///dataset: NKI///series: NKI2///id: 156///er: 1///pgr: NA///her2: NA///size: 2.1///node: 1///age: 42///grade: 3///t.dmfs: 6450///e.dmfs: 0///t.rfs: 6450///e.rfs: 0///t.os: 6446///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                        1 
  samplename: NKI_157///dataset: NKI///series: NKI2///id: 157///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 45///grade: 1///t.dmfs: 2876///e.dmfs: 0///t.rfs: 2876///e.rfs: 0///t.os: 2874///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
samplename: NKI_158///dataset: NKI///series: NKI2///id: 158///er: 0///pgr: NA///her2: NA///size: 2.8///node: 1///age: 49///grade: 3///t.dmfs: 1027///e.dmfs: 1///t.rfs: 1027///e.rfs: 1///t.os: 1426///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
  samplename: NKI_159///dataset: NKI///series: NKI2///id: 159///er: 1///pgr: NA///her2: NA///size: 2///node: 1///age: 44///grade: 2///t.dmfs: 1624///e.dmfs: 1///t.rfs: 1624///e.rfs: 1///t.os: 1977///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                        1 
  samplename: NKI_160///dataset: NKI///series: NKI2///id: 160///er: 1///pgr: NA///her2: NA///size: 4///node: 1///age: 32///grade: 1///t.dmfs: 5898///e.dmfs: 0///t.rfs: 5898///e.rfs: 0///t.os: 5894///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                        1 
  samplename: NKI_161///dataset: NKI///series: NKI2///id: 161///er: 1///pgr: NA///her2: NA///size: 4///node: 1///age: 46///grade: 1///t.dmfs: 2969///e.dmfs: 1///t.rfs: 2969///e.rfs: 1///t.os: 4893///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                        1 
  samplename: NKI_162///dataset: NKI///series: NKI2///id: 162///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 50///grade: 3///t.dmfs: 5593///e.dmfs: 0///t.rfs: 5593///e.rfs: 0///t.os: 5589///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                        1 
samplename: NKI_163///dataset: NKI///series: NKI2///id: 163///er: 1///pgr: NA///her2: NA///size: 1.5///node: 1///age: 38///grade: 3///t.dmfs: 5778///e.dmfs: 0///t.rfs: 5778///e.rfs: 0///t.os: 5774///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                        1 
  samplename: NKI_164///dataset: NKI///series: NKI2///id: 164///er: 0///pgr: NA///her2: NA///size: 3///node: 0///age: 50///grade: 3///t.dmfs: 2069///e.dmfs: 0///t.rfs: 2069///e.rfs: 0///t.os: 2068///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
samplename: NKI_165///dataset: NKI///series: NKI2///id: 165///er: 1///pgr: NA///her2: NA///size: 1.5///node: 0///age: 39///grade: 3///t.dmfs: 3814///e.dmfs: 1///t.rfs: 3814///e.rfs: 1///t.os: 4021///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
  samplename: NKI_166///dataset: NKI///series: NKI2///id: 166///er: 1///pgr: NA///her2: NA///size: 3.5///node: 1///age: 43///grade: 3///t.dmfs: 589///e.dmfs: 1///t.rfs: 589///e.rfs: 1///t.os: 1322///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
samplename: NKI_167///dataset: NKI///series: NKI2///id: 167///er: 1///pgr: NA///her2: NA///size: 1.5///node: 0///age: 44///grade: 1///t.dmfs: 5597///e.dmfs: 0///t.rfs: 5597///e.rfs: 0///t.os: 5593///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
samplename: NKI_169///dataset: NKI///series: NKI2///id: 169///er: 1///pgr: NA///her2: NA///size: 2.5///node: 1///age: 40///grade: 2///t.dmfs: 5437///e.dmfs: 0///t.rfs: 5437///e.rfs: 0///t.os: 5433///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
   samplename: NKI_17///dataset: NKI///series: NKI///id: 17///er: 1///pgr: NA///her2: NA///size: 1.7///node: 0///age: 48///grade: 2///t.dmfs: 2873///e.dmfs: 0///t.rfs: 2873///e.rfs: 0///t.os: 2871///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
samplename: NKI_170///dataset: NKI///series: NKI2///id: 170///er: 1///pgr: NA///her2: NA///size: 2.5///node: 1///age: 42///grade: 1///t.dmfs: 4876///e.dmfs: 0///t.rfs: 4876///e.rfs: 0///t.os: 4873///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                        1 
     samplename: NKI_172///dataset: NKI///series: NKI2///id: 172///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 46///grade: 3///t.dmfs: 507///e.dmfs: 1///t.rfs: 507///e.rfs: 1///t.os: 597///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                        1 
  samplename: NKI_174///dataset: NKI///series: NKI2///id: 174///er: 1///pgr: NA///her2: NA///size: 2///node: 1///age: 41///grade: 3///t.dmfs: 5022///e.dmfs: 0///t.rfs: 5022///e.rfs: 0///t.os: 5019///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
samplename: NKI_175///dataset: NKI///series: NKI2///id: 175///er: 0///pgr: NA///her2: NA///size: 1.5///node: 0///age: 48///grade: 2///t.dmfs: 2774///e.dmfs: 1///t.rfs: 2774///e.rfs: 1///t.os: 2801///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
  samplename: NKI_176///dataset: NKI///series: NKI2///id: 176///er: 1///pgr: NA///her2: NA///size: 5///node: 1///age: 46///grade: 2///t.dmfs: 4592///e.dmfs: 0///t.rfs: 4592///e.rfs: 0///t.os: 4589///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                        1 
samplename: NKI_177///dataset: NKI///series: NKI2///id: 177///er: 0///pgr: NA///her2: NA///size: 2.5///node: 1///age: 48///grade: 3///t.dmfs: 3260///e.dmfs: 1///t.rfs: 3260///e.rfs: 1///t.os: 3545///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                        1 
  samplename: NKI_178///dataset: NKI///series: NKI2///id: 178///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 48///grade: 2///t.dmfs: 4812///e.dmfs: 0///t.rfs: 4812///e.rfs: 0///t.os: 4809///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                        1 
  samplename: NKI_179///dataset: NKI///series: NKI2///id: 179///er: 0///pgr: NA///her2: NA///size: 3///node: 0///age: 40///grade: 3///t.dmfs: 4662///e.dmfs: 0///t.rfs: 4662///e.rfs: 0///t.os: 4659///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
  samplename: NKI_180///dataset: NKI///series: NKI2///id: 180///er: 1///pgr: NA///her2: NA///size: 4.5///node: 1///age: 36///grade: 3///t.dmfs: 955///e.dmfs: 1///t.rfs: 955///e.rfs: 1///t.os: 1930///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                        1 
  samplename: NKI_181///dataset: NKI///series: NKI2///id: 181///er: 0///pgr: NA///her2: NA///size: 2///node: 0///age: 35///grade: 3///t.dmfs: 3951///e.dmfs: 0///t.rfs: 3951///e.rfs: 0///t.os: 4307///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
samplename: NKI_182///dataset: NKI///series: NKI2///id: 182///er: 1///pgr: NA///her2: NA///size: 1.1///node: 1///age: 43///grade: 1///t.dmfs: 4134///e.dmfs: 0///t.rfs: 4134///e.rfs: 0///t.os: 4131///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
samplename: NKI_183///dataset: NKI///series: NKI2///id: 183///er: 1///pgr: NA///her2: NA///size: 2.2///node: 0///age: 42///grade: 1///t.dmfs: 4332///e.dmfs: 0///t.rfs: 4332///e.rfs: 0///t.os: 4329///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
  samplename: NKI_184///dataset: NKI///series: NKI2///id: 184///er: 0///pgr: NA///her2: NA///size: 2.3///node: 1///age: 44///grade: 3///t.dmfs: 442///e.dmfs: 1///t.rfs: 442///e.rfs: 1///t.os: 1608///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                        1 
  samplename: NKI_185///dataset: NKI///series: NKI2///id: 185///er: 1///pgr: NA///her2: NA///size: 2///node: 0///age: 42///grade: 2///t.dmfs: 1950///e.dmfs: 0///t.rfs: 1950///e.rfs: 0///t.os: 2677///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
samplename: NKI_186///dataset: NKI///series: NKI2///id: 186///er: 0///pgr: NA///her2: NA///size: 2.4///node: 1///age: 47///grade: 2///t.dmfs: 4288///e.dmfs: 0///t.rfs: 4288///e.rfs: 0///t.os: 4285///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                        1 
  samplename: NKI_187///dataset: NKI///series: NKI2///id: 187///er: 1///pgr: NA///her2: NA///size: 2///node: 1///age: 40///grade: 2///t.dmfs: 4567///e.dmfs: 0///t.rfs: 4567///e.rfs: 0///t.os: 4564///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                        1 
samplename: NKI_188///dataset: NKI///series: NKI2///id: 188///er: 1///pgr: NA///her2: NA///size: 1.5///node: 1///age: 41///grade: 1///t.dmfs: 4114///e.dmfs: 0///t.rfs: 4114///e.rfs: 0///t.os: 4111///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                        1 
samplename: NKI_189///dataset: NKI///series: NKI2///id: 189///er: 1///pgr: NA///her2: NA///size: 2.5///node: 0///age: 39///grade: 2///t.dmfs: 4410///e.dmfs: 0///t.rfs: 4410///e.rfs: 0///t.os: 4407///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
  samplename: NKI_190///dataset: NKI///series: NKI2///id: 190///er: 1///pgr: NA///her2: NA///size: 2///node: 1///age: 48///grade: 1///t.dmfs: 2721///e.dmfs: 0///t.rfs: 2721///e.rfs: 0///t.os: 4352///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
samplename: NKI_191///dataset: NKI///series: NKI2///id: 191///er: 1///pgr: NA///her2: NA///size: 4.5///node: 0///age: 34///grade: 3///t.dmfs: 4652///e.dmfs: 0///t.rfs: 4652///e.rfs: 0///t.os: 4649///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
    samplename: NKI_192///dataset: NKI///series: NKI2///id: 192///er: 1///pgr: NA///her2: NA///size: 2///node: 1///age: 41///grade: 2///t.dmfs: 985///e.dmfs: 1///t.rfs: 985///e.rfs: 1///t.os: 2298///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                        1 
samplename: NKI_193///dataset: NKI///series: NKI2///id: 193///er: 1///pgr: NA///her2: NA///size: 0.8///node: 0///age: 50///grade: 1///t.dmfs: 4322///e.dmfs: 0///t.rfs: 4322///e.rfs: 0///t.os: 4319///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
samplename: NKI_194///dataset: NKI///series: NKI2///id: 194///er: 1///pgr: NA///her2: NA///size: 1.2///node: 1///age: 45///grade: 2///t.dmfs: 4553///e.dmfs: 1///t.rfs: 4553///e.rfs: 1///t.os: 4770///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                        1 
    samplename: NKI_195///dataset: NKI///series: NKI2///id: 195///er: 0///pgr: NA///her2: NA///size: 2///node: 1///age: 45///grade: 2///t.dmfs: 749///e.dmfs: 0///t.rfs: 749///e.rfs: 0///t.os: 4214///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                        1 
samplename: NKI_196///dataset: NKI///series: NKI2///id: 196///er: 1///pgr: NA///her2: NA///size: 2.5///node: 1///age: 47///grade: 1///t.dmfs: 4089///e.dmfs: 0///t.rfs: 4089///e.rfs: 0///t.os: 4086///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                        1 
  samplename: NKI_197///dataset: NKI///series: NKI2///id: 197///er: 0///pgr: NA///her2: NA///size: 1///node: 1///age: 37///grade: 1///t.dmfs: 4035///e.dmfs: 0///t.rfs: 4035///e.rfs: 0///t.os: 4032///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                        1 
  samplename: NKI_198///dataset: NKI///series: NKI2///id: 198///er: 1///pgr: NA///her2: NA///size: 2///node: 1///age: 43///grade: 2///t.dmfs: 3411///e.dmfs: 0///t.rfs: 3411///e.rfs: 0///t.os: 4067///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                        1 
samplename: NKI_199///dataset: NKI///series: NKI2///id: 199///er: 1///pgr: NA///her2: NA///size: 1.5///node: 0///age: 39///grade: 3///t.dmfs: 3984///e.dmfs: 0///t.rfs: 3984///e.rfs: 0///t.os: 3981///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
  samplename: NKI_200///dataset: NKI///series: NKI2///id: 200///er: 1///pgr: NA///her2: NA///size: 2///node: 1///age: 43///grade: 2///t.dmfs: 3933///e.dmfs: 0///t.rfs: 3933///e.rfs: 0///t.os: 3930///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                        1 
samplename: NKI_201///dataset: NKI///series: NKI2///id: 201///er: 1///pgr: NA///her2: NA///size: 2.5///node: 0///age: 33///grade: 3///t.dmfs: 4091///e.dmfs: 0///t.rfs: 4091///e.rfs: 0///t.os: 4088///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
samplename: NKI_202///dataset: NKI///series: NKI2///id: 202///er: 0///pgr: NA///her2: NA///size: 3.9///node: 0///age: 47///grade: 3///t.dmfs: 1234///e.dmfs: 1///t.rfs: 1234///e.rfs: 1///t.os: 1769///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                        1 
samplename: NKI_203///dataset: NKI///series: NKI2///id: 203///er: 1///pgr: NA///her2: NA///size: 2.3///node: 1///age: 49///grade: 3///t.dmfs: 4031///e.dmfs: 0///t.rfs: 4031///e.rfs: 0///t.os: 4028///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                        1 
                                                                                                                                                                                                                                  (Other) 
                                                                                                                                                                                                                                      238 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20711 }
\keyword{datasets}
