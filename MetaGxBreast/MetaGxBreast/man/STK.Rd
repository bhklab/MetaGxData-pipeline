\name{ STK_eset }
\alias{ STK_eset }
\docType{data}
\title{ STK_eset }
\usage{data( STK_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=12829800 
  Title:  
  URL: http://smd.princeton.edu/cgi-bin/publication/viewPublication.pl?pub_no=248 
  PMIDs: 12829800 
  No abstract available.
  notes:
   summary:     
      Distinct breast cancer subtypes were determined by gene expression profiles and were validated in other published dataset
s.
   version:     
      2015-04-27 19:31:31
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: 1007_s_at 1053_at ... 244889_at (36178 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 36178 features, 159 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
  STK_1  STK_10 STK_101 STK_102 STK_103 STK_106 STK_109  STK_11 STK_113 STK_114 STK_116 STK_117 STK_118 STK_119 STK_121 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
STK_122 STK_123 STK_126 STK_127 STK_128 STK_129  STK_13 STK_130 STK_132 STK_134 STK_135 STK_136 STK_137 STK_138  STK_14 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
STK_140 STK_141 STK_144 STK_145 STK_147 STK_149  STK_15 STK_151 STK_155 STK_156 STK_157 STK_158 STK_159 STK_160 STK_161 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
STK_162 STK_163 STK_165 STK_166 STK_168  STK_17 STK_170 STK_176 STK_177 STK_178 STK_179  STK_18 STK_181 STK_183  STK_19 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
STK_190 STK_191 STK_192 STK_193 STK_194 STK_195 STK_198   STK_2  STK_20 STK_202 STK_203 STK_207 STK_208 STK_209 STK_210 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
STK_211 STK_213 STK_216 STK_217 STK_218  STK_22 STK_222 STK_225 STK_227 STK_228 STK_231 STK_232 STK_235 STK_237 STK_238 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
STK_239 STK_241 STK_243 STK_245 STK_246 STK_249  STK_25 STK_250 STK_251 (Other) 
      1       1       1       1       1       1       1       1       1      60 

alt_sample_name: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    1.0    67.0   136.0   138.3   208.5   277.0 

sample_type: 
tumor 
  159 

er: 
negative positive 
      29      130 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   31.0    48.0    56.0    57.8    68.5    87.0 

grade: 
   1    2    3 NA's 
  28   58   61   12 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
     66    2022    2467    2234    2846    3099 

recurrence_status: 
deceased_or_recurrence    living_norecurrence 
                    46                    113 

treatment: 
  chemotherapy hormonotherapy      untreated 
            89             48             22 

batch: 
STK 
159 

uncurated_author_metadata: 
    samplename: STK_1///dataset: STK///series: STK///id: 1///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 87///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: 3095///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
    samplename: STK_10///dataset: STK///series: STK///id: 10///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: 82///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 832///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                      1 
 samplename: STK_101///dataset: STK///series: STK///id: 101///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: 64///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2183///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                      1 
 samplename: STK_102///dataset: STK///series: STK///id: 102///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: 63///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2303///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                      1 
 samplename: STK_103///dataset: STK///series: STK///id: 103///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 63///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2307///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_106///dataset: STK///series: STK///id: 106///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 63///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2201///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_109///dataset: STK///series: STK///id: 109///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 60///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 3092///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                      1 
   samplename: STK_11///dataset: STK///series: STK///id: 11///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 83///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1018///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_113///dataset: STK///series: STK///id: 113///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 62///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2077///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
  samplename: STK_114///dataset: STK///series: STK///id: 114///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 62///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 475///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                      1 
 samplename: STK_116///dataset: STK///series: STK///id: 116///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 58///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 3037///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                      1 
 samplename: STK_117///dataset: STK///series: STK///id: 117///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: 57///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 3030///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                      1 
 samplename: STK_118///dataset: STK///series: STK///id: 118///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 57///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2427///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                      1 
  samplename: STK_119///dataset: STK///series: STK///id: 119///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 58///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 675///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_121///dataset: STK///series: STK///id: 121///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 57///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2931///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                      1 
 samplename: STK_122///dataset: STK///series: STK///id: 122///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 57///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2770///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                      1 
 samplename: STK_123///dataset: STK///series: STK///id: 123///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: 57///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2694///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                      1 
 samplename: STK_126///dataset: STK///series: STK///id: 126///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 57///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2391///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                      1 
 samplename: STK_127///dataset: STK///series: STK///id: 127///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 56///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2738///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_128///dataset: STK///series: STK///id: 128///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 56///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2949///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_129///dataset: STK///series: STK///id: 129///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 58///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2139///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
  samplename: STK_13///dataset: STK///series: STK///id: 13///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 80///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: 3066///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_130///dataset: STK///series: STK///id: 130///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 58///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2208///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_132///dataset: STK///series: STK///id: 132///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 56///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2657///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_134///dataset: STK///series: STK///id: 134///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: 55///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2278///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                      1 
 samplename: STK_135///dataset: STK///series: STK///id: 135///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 58///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2051///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_136///dataset: STK///series: STK///id: 136///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 56///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1591///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_137///dataset: STK///series: STK///id: 137///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: 56///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2661///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                      1 
 samplename: STK_138///dataset: STK///series: STK///id: 138///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 56///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2632///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
   samplename: STK_14///dataset: STK///series: STK///id: 14///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 80///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 3037///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
samplename: STK_140///dataset: STK///series: STK///id: 140///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: 54///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: 3008///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                      1 
 samplename: STK_141///dataset: STK///series: STK///id: 141///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 56///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2336///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_144///dataset: STK///series: STK///id: 144///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: 54///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 3051///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                      1 
 samplename: STK_145///dataset: STK///series: STK///id: 145///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: 53///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2891///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_147///dataset: STK///series: STK///id: 147///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 54///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2800///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                      1 
 samplename: STK_149///dataset: STK///series: STK///id: 149///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 53///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2960///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
   samplename: STK_15///dataset: STK///series: STK///id: 15///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 81///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2237///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_151///dataset: STK///series: STK///id: 151///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 53///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2741///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                      1 
samplename: STK_155///dataset: STK///series: STK///id: 155///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 55///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: 2121///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                      1 
 samplename: STK_156///dataset: STK///series: STK///id: 156///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 54///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2084///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                      1 
 samplename: STK_157///dataset: STK///series: STK///id: 157///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 54///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2055///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                      1 
 samplename: STK_158///dataset: STK///series: STK///id: 158///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 53///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2657///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                      1 
 samplename: STK_159///dataset: STK///series: STK///id: 159///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 52///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2018///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_160///dataset: STK///series: STK///id: 160///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: 52///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2967///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
  samplename: STK_161///dataset: STK///series: STK///id: 161///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 52///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 380///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_162///dataset: STK///series: STK///id: 162///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 53///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2581///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                      1 
 samplename: STK_163///dataset: STK///series: STK///id: 163///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 52///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2931///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_165///dataset: STK///series: STK///id: 165///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 52///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2814///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_166///dataset: STK///series: STK///id: 166///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 51///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2916///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_168///dataset: STK///series: STK///id: 168///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 51///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2946///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
   samplename: STK_17///dataset: STK///series: STK///id: 17///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 79///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2829///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_170///dataset: STK///series: STK///id: 170///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 50///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 3040///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                      1 
samplename: STK_176///dataset: STK///series: STK///id: 176///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 52///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: 2358///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_177///dataset: STK///series: STK///id: 177///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 52///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2135///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                      1 
 samplename: STK_178///dataset: STK///series: STK///id: 178///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 52///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2077///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                      1 
 samplename: STK_179///dataset: STK///series: STK///id: 179///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 50///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2821///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
   samplename: STK_18///dataset: STK///series: STK///id: 18///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 78///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 1394///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_181///dataset: STK///series: STK///id: 181///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 50///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1150///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_183///dataset: STK///series: STK///id: 183///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 50///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2686///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
    samplename: STK_19///dataset: STK///series: STK///id: 19///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 78///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 810///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
  samplename: STK_190///dataset: STK///series: STK///id: 190///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: 50///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 409///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                      1 
 samplename: STK_191///dataset: STK///series: STK///id: 191///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 50///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2581///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                      1 
 samplename: STK_192///dataset: STK///series: STK///id: 192///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 51///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2172///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                      1 
 samplename: STK_193///dataset: STK///series: STK///id: 193///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 49///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2975///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
  samplename: STK_194///dataset: STK///series: STK///id: 194///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 49///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 442///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                      1 
 samplename: STK_195///dataset: STK///series: STK///id: 195///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 49///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2869///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_198///dataset: STK///series: STK///id: 198///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 49///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2862///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
    samplename: STK_2///dataset: STK///series: STK///id: 2///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 86///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: 2026///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
    samplename: STK_20///dataset: STK///series: STK///id: 20///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: 78///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 628///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                      1 
 samplename: STK_202///dataset: STK///series: STK///id: 202///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 49///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1431///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                      1 
 samplename: STK_203///dataset: STK///series: STK///id: 203///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 51///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2073///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_207///dataset: STK///series: STK///id: 207///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 50///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2139///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                      1 
 samplename: STK_208///dataset: STK///series: STK///id: 208///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 47///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 3099///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_209///dataset: STK///series: STK///id: 209///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 48///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2829///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                      1 
 samplename: STK_210///dataset: STK///series: STK///id: 210///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 48///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2716///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                      1 
 samplename: STK_211///dataset: STK///series: STK///id: 211///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 48///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2840///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                      1 
 samplename: STK_213///dataset: STK///series: STK///id: 213///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 47///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 3037///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
  samplename: STK_216///dataset: STK///series: STK///id: 216///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 47///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 398///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                      1 
 samplename: STK_217///dataset: STK///series: STK///id: 217///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 47///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2821///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_218///dataset: STK///series: STK///id: 218///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 47///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2756///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                      1 
    samplename: STK_22///dataset: STK///series: STK///id: 22///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 77///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 818///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_222///dataset: STK///series: STK///id: 222///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: 48///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2197///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                      1 
 samplename: STK_225///dataset: STK///series: STK///id: 225///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 46///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 3077///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_227///dataset: STK///series: STK///id: 227///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 46///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2902///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                      1 
 samplename: STK_228///dataset: STK///series: STK///id: 228///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 46///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2905///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_231///dataset: STK///series: STK///id: 231///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 45///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2931///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_232///dataset: STK///series: STK///id: 232///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: 46///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2504///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                      1 
 samplename: STK_235///dataset: STK///series: STK///id: 235///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 45///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2785///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_237///dataset: STK///series: STK///id: 237///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: 46///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1778///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                      1 
 samplename: STK_238///dataset: STK///series: STK///id: 238///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 43///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2840///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_239///dataset: STK///series: STK///id: 239///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 44///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2573///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_241///dataset: STK///series: STK///id: 241///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 43///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2756///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_243///dataset: STK///series: STK///id: 243///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: 42///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2814///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_245///dataset: STK///series: STK///id: 245///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: 42///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2624///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_246///dataset: STK///series: STK///id: 246///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: 41///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1559///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_249///dataset: STK///series: STK///id: 249///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 42///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2467///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
   samplename: STK_25///dataset: STK///series: STK///id: 25///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 78///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2566///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
 samplename: STK_250///dataset: STK///series: STK///id: 250///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: 42///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2329///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                      1 
 samplename: STK_251///dataset: STK///series: STK///id: 251///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: 40///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1409///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                      1 
                                                                                                                                                                                                                                (Other) 
                                                                                                                                                                                                                                     60 

}}
\source{ http://smd.princeton.edu/cgi-bin/publication/viewPublication.pl?pub_no=248 }
\keyword{datasets}
