\name{ EMC2_eset }
\alias{ EMC2_eset }
\docType{data}
\title{  }
\usage{data( EMC2_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/19421193 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse12276 
  PMIDs: 19421193 
  No abstract available.
  notes:
   summary:     
      Genes were identified that may increase the ability of breast cancer cells to infiltrate the blood-brain barrier.
   version:     
      2015-04-27 18:59:10
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
assayData: 42447 features, 204 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
  EMC2_1  EMC2_10 EMC2_100 EMC2_101 EMC2_102 EMC2_103 EMC2_104 EMC2_105 EMC2_106 EMC2_107 EMC2_108 EMC2_109  EMC2_11 EMC2_110 
       1        1        1        1        1        1        1        1        1        1        1        1        1        1 
EMC2_111 EMC2_112 EMC2_113 EMC2_114 EMC2_115 EMC2_116 EMC2_117 EMC2_118 EMC2_119  EMC2_12 EMC2_120 EMC2_121 EMC2_122 EMC2_123 
       1        1        1        1        1        1        1        1        1        1        1        1        1        1 
EMC2_124 EMC2_125 EMC2_126 EMC2_127 EMC2_128 EMC2_129  EMC2_13 EMC2_130 EMC2_131 EMC2_132 EMC2_133 EMC2_134 EMC2_135 EMC2_136 
       1        1        1        1        1        1        1        1        1        1        1        1        1        1 
EMC2_137 EMC2_138 EMC2_139  EMC2_14 EMC2_140 EMC2_141 EMC2_142 EMC2_143 EMC2_144 EMC2_145 EMC2_146 EMC2_147 EMC2_148 EMC2_149 
       1        1        1        1        1        1        1        1        1        1        1        1        1        1 
 EMC2_15 EMC2_150 EMC2_151 EMC2_152 EMC2_153 EMC2_154 EMC2_155 EMC2_156 EMC2_157 EMC2_158 EMC2_159  EMC2_16 EMC2_160 EMC2_161 
       1        1        1        1        1        1        1        1        1        1        1        1        1        1 
EMC2_162 EMC2_163 EMC2_164 EMC2_165 EMC2_166 EMC2_167 EMC2_168 EMC2_169  EMC2_17 EMC2_170 EMC2_171 EMC2_172 EMC2_173 EMC2_174 
       1        1        1        1        1        1        1        1        1        1        1        1        1        1 
EMC2_175 EMC2_176 EMC2_177 EMC2_178 EMC2_179  EMC2_18 EMC2_180 EMC2_181 EMC2_182 EMC2_183 EMC2_184 EMC2_185 EMC2_186 EMC2_187 
       1        1        1        1        1        1        1        1        1        1        1        1        1        1 
EMC2_188  (Other) 
       1      105 

alt_sample_name: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   1.00   51.75  102.50  102.50  153.20  204.00 

sample_type: 
tumor 
  204 

N: 
   0 NA's 
  48  156 

dmfs_days: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
      0     335     640     799    1098    3507 

dmfs_status: 
deceased_or_recurrence    living_norecurrence 
                   185                     19 

treatment: 
chemotherapy    untreated 
         156           48 

batch: 
EMC2 
 204 

uncurated_author_metadata: 
      samplename: EMC2_1///dataset: EMC2///series: EMC2///id: 1///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 91///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
   samplename: EMC2_10///dataset: EMC2///series: EMC2///id: 10///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 701///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_100///dataset: EMC2///series: EMC2///id: 100///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 549///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_101///dataset: EMC2///series: EMC2///id: 101///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 488///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
   samplename: EMC2_102///dataset: EMC2///series: EMC2///id: 102///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 0///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
samplename: EMC2_103///dataset: EMC2///series: EMC2///id: 103///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 1098///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
samplename: EMC2_104///dataset: EMC2///series: EMC2///id: 104///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 1189///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_105///dataset: EMC2///series: EMC2///id: 105///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 823///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_106///dataset: EMC2///series: EMC2///id: 106///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 457///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_107///dataset: EMC2///series: EMC2///id: 107///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 671///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_108///dataset: EMC2///series: EMC2///id: 108///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 884///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_109///dataset: EMC2///series: EMC2///id: 109///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 335///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
    samplename: EMC2_11///dataset: EMC2///series: EMC2///id: 11///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 91///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
samplename: EMC2_110///dataset: EMC2///series: EMC2///id: 110///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 1769///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
samplename: EMC2_111///dataset: EMC2///series: EMC2///id: 111///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 1220///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
samplename: EMC2_112///dataset: EMC2///series: EMC2///id: 112///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 1464///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
samplename: EMC2_113///dataset: EMC2///series: EMC2///id: 113///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 1159///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_114///dataset: EMC2///series: EMC2///id: 114///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 366///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_115///dataset: EMC2///series: EMC2///id: 115///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 396///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
samplename: EMC2_116///dataset: EMC2///series: EMC2///id: 116///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 2592///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_117///dataset: EMC2///series: EMC2///id: 117///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 274///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_118///dataset: EMC2///series: EMC2///id: 118///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 183///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_119///dataset: EMC2///series: EMC2///id: 119///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 854///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
   samplename: EMC2_12///dataset: EMC2///series: EMC2///id: 12///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 518///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_120///dataset: EMC2///series: EMC2///id: 120///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 183///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_121///dataset: EMC2///series: EMC2///id: 121///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 579///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_122///dataset: EMC2///series: EMC2///id: 122///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 457///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_123///dataset: EMC2///series: EMC2///id: 123///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 579///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_124///dataset: EMC2///series: EMC2///id: 124///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 213///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
samplename: EMC2_125///dataset: EMC2///series: EMC2///id: 125///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 3507///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
samplename: EMC2_126///dataset: EMC2///series: EMC2///id: 126///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 1403///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_127///dataset: EMC2///series: EMC2///id: 127///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 274///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_128///dataset: EMC2///series: EMC2///id: 128///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 976///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_129///dataset: EMC2///series: EMC2///id: 129///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 823///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
  samplename: EMC2_13///dataset: EMC2///series: EMC2///id: 13///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 1159///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_130///dataset: EMC2///series: EMC2///id: 130///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 152///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_131///dataset: EMC2///series: EMC2///id: 131///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 122///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_132///dataset: EMC2///series: EMC2///id: 132///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 518///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_133///dataset: EMC2///series: EMC2///id: 133///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 457///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
samplename: EMC2_134///dataset: EMC2///series: EMC2///id: 134///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 1860///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_135///dataset: EMC2///series: EMC2///id: 135///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 213///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_136///dataset: EMC2///series: EMC2///id: 136///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 183///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
samplename: EMC2_137///dataset: EMC2///series: EMC2///id: 137///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 1006///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_138///dataset: EMC2///series: EMC2///id: 138///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 671///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
samplename: EMC2_139///dataset: EMC2///series: EMC2///id: 139///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 1037///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
   samplename: EMC2_14///dataset: EMC2///series: EMC2///id: 14///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 244///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_140///dataset: EMC2///series: EMC2///id: 140///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 427///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
samplename: EMC2_141///dataset: EMC2///series: EMC2///id: 141///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 1311///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_142///dataset: EMC2///series: EMC2///id: 142///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 366///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
samplename: EMC2_143///dataset: EMC2///series: EMC2///id: 143///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 1159///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_144///dataset: EMC2///series: EMC2///id: 144///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 335///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_145///dataset: EMC2///series: EMC2///id: 145///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 396///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_146///dataset: EMC2///series: EMC2///id: 146///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 366///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_147///dataset: EMC2///series: EMC2///id: 147///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 335///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_148///dataset: EMC2///series: EMC2///id: 148///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 244///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_149///dataset: EMC2///series: EMC2///id: 149///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 732///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
   samplename: EMC2_15///dataset: EMC2///series: EMC2///id: 15///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 213///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
samplename: EMC2_150///dataset: EMC2///series: EMC2///id: 150///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 1769///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
samplename: EMC2_151///dataset: EMC2///series: EMC2///id: 151///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 2867///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_152///dataset: EMC2///series: EMC2///id: 152///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 122///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
   samplename: EMC2_153///dataset: EMC2///series: EMC2///id: 153///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 0///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
samplename: EMC2_154///dataset: EMC2///series: EMC2///id: 154///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 1159///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_155///dataset: EMC2///series: EMC2///id: 155///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 335///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
samplename: EMC2_156///dataset: EMC2///series: EMC2///id: 156///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 1464///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_157///dataset: EMC2///series: EMC2///id: 157///er: NA///pgr: NA///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: 1494///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_158///dataset: EMC2///series: EMC2///id: 158///er: NA///pgr: NA///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: 2958///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_159///dataset: EMC2///series: EMC2///id: 159///er: NA///pgr: NA///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: 1525///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                          1 
   samplename: EMC2_16///dataset: EMC2///series: EMC2///id: 16///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 884///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_160///dataset: EMC2///series: EMC2///id: 160///er: NA///pgr: NA///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: 1616///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                          1 
  samplename: EMC2_161///dataset: EMC2///series: EMC2///id: 161///er: NA///pgr: NA///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: 701///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                          1 
  samplename: EMC2_162///dataset: EMC2///series: EMC2///id: 162///er: NA///pgr: NA///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: 762///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                          1 
  samplename: EMC2_163///dataset: EMC2///series: EMC2///id: 163///er: NA///pgr: NA///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: 671///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_164///dataset: EMC2///series: EMC2///id: 164///er: NA///pgr: NA///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: 1098///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_165///dataset: EMC2///series: EMC2///id: 165///er: NA///pgr: NA///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: 1067///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                          1 
  samplename: EMC2_166///dataset: EMC2///series: EMC2///id: 166///er: NA///pgr: NA///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: 488///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_167///dataset: EMC2///series: EMC2///id: 167///er: NA///pgr: NA///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: 1067///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                          1 
  samplename: EMC2_168///dataset: EMC2///series: EMC2///id: 168///er: NA///pgr: NA///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: 701///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                          1 
  samplename: EMC2_169///dataset: EMC2///series: EMC2///id: 169///er: NA///pgr: NA///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: 427///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                          1 
   samplename: EMC2_17///dataset: EMC2///series: EMC2///id: 17///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 427///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
  samplename: EMC2_170///dataset: EMC2///series: EMC2///id: 170///er: NA///pgr: NA///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: 274///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_171///dataset: EMC2///series: EMC2///id: 171///er: NA///pgr: NA///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: 1281///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                          1 
  samplename: EMC2_172///dataset: EMC2///series: EMC2///id: 172///er: NA///pgr: NA///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: 518///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                          1 
  samplename: EMC2_173///dataset: EMC2///series: EMC2///id: 173///er: NA///pgr: NA///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: 152///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_174///dataset: EMC2///series: EMC2///id: 174///er: NA///pgr: NA///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: 1677///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                          1 
  samplename: EMC2_175///dataset: EMC2///series: EMC2///id: 175///er: NA///pgr: NA///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: 701///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_176///dataset: EMC2///series: EMC2///id: 176///er: NA///pgr: NA///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: 1220///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_177///dataset: EMC2///series: EMC2///id: 177///er: NA///pgr: NA///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: 1342///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_178///dataset: EMC2///series: EMC2///id: 178///er: NA///pgr: NA///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: 1433///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                          1 
  samplename: EMC2_179///dataset: EMC2///series: EMC2///id: 179///er: NA///pgr: NA///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: 793///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                          1 
   samplename: EMC2_18///dataset: EMC2///series: EMC2///id: 18///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: 579///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                          1 
  samplename: EMC2_180///dataset: EMC2///series: EMC2///id: 180///er: NA///pgr: NA///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: 640///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_181///dataset: EMC2///series: EMC2///id: 181///er: NA///pgr: NA///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: 1067///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_182///dataset: EMC2///series: EMC2///id: 182///er: NA///pgr: NA///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: 1067///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                          1 
  samplename: EMC2_183///dataset: EMC2///series: EMC2///id: 183///er: NA///pgr: NA///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: 122///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                          1 
  samplename: EMC2_184///dataset: EMC2///series: EMC2///id: 184///er: NA///pgr: NA///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: 152///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                          1 
 samplename: EMC2_185///dataset: EMC2///series: EMC2///id: 185///er: NA///pgr: NA///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: 1281///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                          1 
  samplename: EMC2_186///dataset: EMC2///series: EMC2///id: 186///er: NA///pgr: NA///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: 549///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                          1 
  samplename: EMC2_187///dataset: EMC2///series: EMC2///id: 187///er: NA///pgr: NA///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: 427///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                          1 
  samplename: EMC2_188///dataset: EMC2///series: EMC2///id: 188///er: NA///pgr: NA///her2: NA///size: NA///node: 0///age: NA///grade: NA///t.dmfs: 488///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                          1 
                                                                                                                                                                                                                                    (Other) 
                                                                                                                                                                                                                                        105 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse12276 }
\keyword{datasets}
