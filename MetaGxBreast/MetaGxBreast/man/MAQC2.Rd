\name{ MAQC2_eset }
\alias{ MAQC2_eset }
\docType{data}
\title{ MAQC2_eset }
\usage{data( MAQC2_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=20064235 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16716 
  PMIDs: 20064235 
  No abstract available.
  notes:
   summary:     
      It is possible to build multi-gene classifiers of clinical outcome. Prediction accuracy depends on training sample size a
nd classification difficulty. 
   version:     
      2015-04-27 19:17:13
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
assayData: 20967 features, 230 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
MAQC2_106 MAQC2_107 MAQC2_108 MAQC2_111 MAQC2_113 MAQC2_116 MAQC2_117 MAQC2_120 MAQC2_121 MAQC2_123 MAQC2_126 MAQC2_128 
        1         1         1         1         1         1         1         1         1         1         1         1 
MAQC2_129 MAQC2_130 MAQC2_133 MAQC2_135 MAQC2_136 MAQC2_139 MAQC2_141 MAQC2_145 MAQC2_146 MAQC2_153 MAQC2_154 MAQC2_155 
        1         1         1         1         1         1         1         1         1         1         1         1 
MAQC2_156 MAQC2_157 MAQC2_158 MAQC2_159 MAQC2_161 MAQC2_165 MAQC2_176 MAQC2_177 MAQC2_179 MAQC2_180 MAQC2_181 MAQC2_182 
        1         1         1         1         1         1         1         1         1         1         1         1 
MAQC2_186 MAQC2_188 MAQC2_189 MAQC2_195 MAQC2_196 MAQC2_199 MAQC2_201 MAQC2_205 MAQC2_206 MAQC2_211 MAQC2_212 MAQC2_214 
        1         1         1         1         1         1         1         1         1         1         1         1 
MAQC2_215 MAQC2_216 MAQC2_217 MAQC2_220 MAQC2_226 MAQC2_227 MAQC2_228 MAQC2_230 MAQC2_231 MAQC2_232 MAQC2_233 MAQC2_234 
        1         1         1         1         1         1         1         1         1         1         1         1 
MAQC2_235 MAQC2_236 MAQC2_237 MAQC2_238 MAQC2_239 MAQC2_245 MAQC2_246 MAQC2_251 MAQC2_255 MAQC2_256 MAQC2_257 MAQC2_258 
        1         1         1         1         1         1         1         1         1         1         1         1 
MAQC2_259 MAQC2_260 MAQC2_261 MAQC2_263 MAQC2_264 MAQC2_266 MAQC2_270 MAQC2_280 MAQC2_281 MAQC2_282 MAQC2_283 MAQC2_286 
        1         1         1         1         1         1         1         1         1         1         1         1 
MAQC2_287 MAQC2_294 MAQC2_295 MAQC2_297 MAQC2_301 MAQC2_302 MAQC2_304 MAQC2_309 MAQC2_310 MAQC2_315 MAQC2_316 MAQC2_319 
        1         1         1         1         1         1         1         1         1         1         1         1 
MAQC2_322 MAQC2_323 MAQC2_330   (Other) 
        1         1         1       131 

alt_sample_name: 
    106     107     108     111     113     116     117     120     121     123     126     128     129     130     133 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
    135     136     139     141     145     146     153     154     155     156     157     158     159     161     165 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
    176     177     179     180     181     182     186     188     189     195     196     199     201     205     206 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
    211     212     214     215     216     217     220     226     227     228     230     231     232     233     234 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
    235     236     237     238     239     245     246     251     255     256     257     258     259     260     261 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
    263     264     266     270     280     281     282     283     286     287     294     295     297     301     302 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
    304     309     310     315     316     319     322     323     330 (Other) 
      1       1       1       1       1       1       1       1       1     131 

sample_type: 
tumor 
  230 

er: 
negative positive 
      89      141 

pgr: 
negative positive 
     126      104 

her2: 
negative positive 
     190       40 

N: 
  0   1 
 66 164 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  26.00   45.00   51.00   52.02   59.00   79.00 

grade: 
  1   2   3 
 13  94 123 

treatment: 
chemotherapy 
         230 

batch: 
MAQC2 
  230 

uncurated_author_metadata: 
samplename: MAQC2_106///dataset: MAQC2///series: MAQC2///id: 106///er: 1///pgr: 1///her2: 0///size: NA///node: 0///age: 57///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_107///dataset: MAQC2///series: MAQC2///id: 107///er: 1///pgr: 0///her2: 0///size: NA///node: 1///age: 47///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_108///dataset: MAQC2///series: MAQC2///id: 108///er: 1///pgr: 1///her2: 1///size: NA///node: 0///age: 39///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_111///dataset: MAQC2///series: MAQC2///id: 111///er: 0///pgr: 0///her2: 0///size: NA///node: 1///age: 50///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_113///dataset: MAQC2///series: MAQC2///id: 113///er: 0///pgr: 0///her2: 0///size: NA///node: 0///age: 75///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_116///dataset: MAQC2///series: MAQC2///id: 116///er: 0///pgr: 0///her2: 1///size: NA///node: 1///age: 49///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_117///dataset: MAQC2///series: MAQC2///id: 117///er: 1///pgr: 1///her2: 0///size: NA///node: 1///age: 34///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_120///dataset: MAQC2///series: MAQC2///id: 120///er: 0///pgr: 1///her2: 0///size: NA///node: 0///age: 52///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_121///dataset: MAQC2///series: MAQC2///id: 121///er: 0///pgr: 0///her2: 0///size: NA///node: 0///age: 65///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_123///dataset: MAQC2///series: MAQC2///id: 123///er: 1///pgr: 0///her2: 0///size: NA///node: 1///age: 52///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_126///dataset: MAQC2///series: MAQC2///id: 126///er: 0///pgr: 1///her2: 0///size: NA///node: 0///age: 59///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_128///dataset: MAQC2///series: MAQC2///id: 128///er: 1///pgr: 1///her2: 0///size: NA///node: 0///age: 67///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_129///dataset: MAQC2///series: MAQC2///id: 129///er: 0///pgr: 0///her2: 0///size: NA///node: 0///age: 61///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_130///dataset: MAQC2///series: MAQC2///id: 130///er: 0///pgr: 1///her2: 0///size: NA///node: 1///age: 66///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_133///dataset: MAQC2///series: MAQC2///id: 133///er: 0///pgr: 0///her2: 0///size: NA///node: 0///age: 38///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_135///dataset: MAQC2///series: MAQC2///id: 135///er: 1///pgr: 0///her2: 1///size: NA///node: 0///age: 51///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_136///dataset: MAQC2///series: MAQC2///id: 136///er: 0///pgr: 1///her2: 0///size: NA///node: 0///age: 38///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_139///dataset: MAQC2///series: MAQC2///id: 139///er: 1///pgr: 1///her2: 0///size: NA///node: 0///age: 51///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_141///dataset: MAQC2///series: MAQC2///id: 141///er: 1///pgr: 0///her2: 1///size: NA///node: 1///age: 70///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_145///dataset: MAQC2///series: MAQC2///id: 145///er: 1///pgr: 1///her2: 1///size: NA///node: 1///age: 45///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_146///dataset: MAQC2///series: MAQC2///id: 146///er: 0///pgr: 0///her2: 0///size: NA///node: 0///age: 40///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_153///dataset: MAQC2///series: MAQC2///id: 153///er: 0///pgr: 1///her2: 0///size: NA///node: 1///age: 61///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_154///dataset: MAQC2///series: MAQC2///id: 154///er: 0///pgr: 0///her2: 0///size: NA///node: 1///age: 29///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_155///dataset: MAQC2///series: MAQC2///id: 155///er: 1///pgr: 1///her2: 0///size: NA///node: 0///age: 62///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_156///dataset: MAQC2///series: MAQC2///id: 156///er: 0///pgr: 0///her2: 0///size: NA///node: 1///age: 67///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_157///dataset: MAQC2///series: MAQC2///id: 157///er: 1///pgr: 1///her2: 0///size: NA///node: 0///age: 57///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_158///dataset: MAQC2///series: MAQC2///id: 158///er: 1///pgr: 1///her2: 0///size: NA///node: 1///age: 38///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_159///dataset: MAQC2///series: MAQC2///id: 159///er: 0///pgr: 0///her2: 1///size: NA///node: 0///age: 59///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_161///dataset: MAQC2///series: MAQC2///id: 161///er: 0///pgr: 1///her2: 1///size: NA///node: 1///age: 44///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_165///dataset: MAQC2///series: MAQC2///id: 165///er: 1///pgr: 0///her2: 0///size: NA///node: 1///age: 50///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_176///dataset: MAQC2///series: MAQC2///id: 176///er: 1///pgr: 0///her2: 1///size: NA///node: 1///age: 77///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_177///dataset: MAQC2///series: MAQC2///id: 177///er: 1///pgr: 0///her2: 0///size: NA///node: 0///age: 42///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_179///dataset: MAQC2///series: MAQC2///id: 179///er: 0///pgr: 0///her2: 1///size: NA///node: 1///age: 54///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_180///dataset: MAQC2///series: MAQC2///id: 180///er: 0///pgr: 0///her2: 0///size: NA///node: 1///age: 44///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_181///dataset: MAQC2///series: MAQC2///id: 181///er: 1///pgr: 1///her2: 0///size: NA///node: 1///age: 75///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_182///dataset: MAQC2///series: MAQC2///id: 182///er: 1///pgr: 1///her2: 0///size: NA///node: 0///age: 46///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_186///dataset: MAQC2///series: MAQC2///id: 186///er: 0///pgr: 0///her2: 1///size: NA///node: 0///age: 54///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_188///dataset: MAQC2///series: MAQC2///id: 188///er: 1///pgr: 1///her2: 0///size: NA///node: 1///age: 68///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_189///dataset: MAQC2///series: MAQC2///id: 189///er: 1///pgr: 1///her2: 0///size: NA///node: 0///age: 47///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_195///dataset: MAQC2///series: MAQC2///id: 195///er: 1///pgr: 1///her2: 0///size: NA///node: 0///age: 42///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_196///dataset: MAQC2///series: MAQC2///id: 196///er: 1///pgr: 1///her2: 0///size: NA///node: 1///age: 69///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_199///dataset: MAQC2///series: MAQC2///id: 199///er: 0///pgr: 0///her2: 1///size: NA///node: 1///age: 46///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_201///dataset: MAQC2///series: MAQC2///id: 201///er: 1///pgr: 0///her2: 0///size: NA///node: 1///age: 75///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_205///dataset: MAQC2///series: MAQC2///id: 205///er: 0///pgr: 0///her2: 0///size: NA///node: 1///age: 35///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_206///dataset: MAQC2///series: MAQC2///id: 206///er: 0///pgr: 0///her2: 0///size: NA///node: 0///age: 48///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_211///dataset: MAQC2///series: MAQC2///id: 211///er: 0///pgr: 0///her2: 0///size: NA///node: 1///age: 52///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_212///dataset: MAQC2///series: MAQC2///id: 212///er: 0///pgr: 0///her2: 1///size: NA///node: 1///age: 42///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_214///dataset: MAQC2///series: MAQC2///id: 214///er: 1///pgr: 0///her2: 0///size: NA///node: 1///age: 54///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_215///dataset: MAQC2///series: MAQC2///id: 215///er: 0///pgr: 0///her2: 0///size: NA///node: 1///age: 44///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_216///dataset: MAQC2///series: MAQC2///id: 216///er: 1///pgr: 0///her2: 0///size: NA///node: 1///age: 59///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_217///dataset: MAQC2///series: MAQC2///id: 217///er: 1///pgr: 1///her2: 0///size: NA///node: 0///age: 52///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_220///dataset: MAQC2///series: MAQC2///id: 220///er: 1///pgr: 0///her2: 0///size: NA///node: 1///age: 38///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_226///dataset: MAQC2///series: MAQC2///id: 226///er: 1///pgr: 1///her2: 0///size: NA///node: 1///age: 31///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_227///dataset: MAQC2///series: MAQC2///id: 227///er: 1///pgr: 1///her2: 0///size: NA///node: 1///age: 38///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_228///dataset: MAQC2///series: MAQC2///id: 228///er: 1///pgr: 0///her2: 0///size: NA///node: 1///age: 53///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_230///dataset: MAQC2///series: MAQC2///id: 230///er: 0///pgr: 0///her2: 1///size: NA///node: 1///age: 58///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_231///dataset: MAQC2///series: MAQC2///id: 231///er: 1///pgr: 1///her2: 0///size: NA///node: 0///age: 52///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_232///dataset: MAQC2///series: MAQC2///id: 232///er: 0///pgr: 0///her2: 1///size: NA///node: 1///age: 49///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_233///dataset: MAQC2///series: MAQC2///id: 233///er: 0///pgr: 0///her2: 0///size: NA///node: 1///age: 60///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_234///dataset: MAQC2///series: MAQC2///id: 234///er: 1///pgr: 0///her2: 0///size: NA///node: 1///age: 58///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_235///dataset: MAQC2///series: MAQC2///id: 235///er: 0///pgr: 0///her2: 1///size: NA///node: 1///age: 53///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_236///dataset: MAQC2///series: MAQC2///id: 236///er: 1///pgr: 1///her2: 0///size: NA///node: 1///age: 53///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_237///dataset: MAQC2///series: MAQC2///id: 237///er: 1///pgr: 0///her2: 1///size: NA///node: 1///age: 51///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_238///dataset: MAQC2///series: MAQC2///id: 238///er: 1///pgr: 1///her2: 1///size: NA///node: 0///age: 56///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_239///dataset: MAQC2///series: MAQC2///id: 239///er: 0///pgr: 0///her2: 1///size: NA///node: 0///age: 55///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_245///dataset: MAQC2///series: MAQC2///id: 245///er: 1///pgr: 1///her2: 1///size: NA///node: 1///age: 52///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_246///dataset: MAQC2///series: MAQC2///id: 246///er: 1///pgr: 0///her2: 0///size: NA///node: 1///age: 42///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_251///dataset: MAQC2///series: MAQC2///id: 251///er: 1///pgr: 1///her2: 1///size: NA///node: 1///age: 79///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_255///dataset: MAQC2///series: MAQC2///id: 255///er: 0///pgr: 0///her2: 0///size: NA///node: 1///age: 48///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_256///dataset: MAQC2///series: MAQC2///id: 256///er: 1///pgr: 1///her2: 0///size: NA///node: 1///age: 79///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_257///dataset: MAQC2///series: MAQC2///id: 257///er: 1///pgr: 1///her2: 1///size: NA///node: 0///age: 73///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_258///dataset: MAQC2///series: MAQC2///id: 258///er: 0///pgr: 0///her2: 1///size: NA///node: 1///age: 39///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_259///dataset: MAQC2///series: MAQC2///id: 259///er: 1///pgr: 1///her2: 1///size: NA///node: 0///age: 39///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_260///dataset: MAQC2///series: MAQC2///id: 260///er: 1///pgr: 1///her2: 0///size: NA///node: 1///age: 50///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_261///dataset: MAQC2///series: MAQC2///id: 261///er: 1///pgr: 1///her2: 0///size: NA///node: 0///age: 56///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_263///dataset: MAQC2///series: MAQC2///id: 263///er: 1///pgr: 1///her2: 0///size: NA///node: 1///age: 72///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_264///dataset: MAQC2///series: MAQC2///id: 264///er: 1///pgr: 1///her2: 1///size: NA///node: 1///age: 62///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_266///dataset: MAQC2///series: MAQC2///id: 266///er: 1///pgr: 0///her2: 0///size: NA///node: 1///age: 45///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_270///dataset: MAQC2///series: MAQC2///id: 270///er: 1///pgr: 0///her2: 0///size: NA///node: 0///age: 67///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_280///dataset: MAQC2///series: MAQC2///id: 280///er: 1///pgr: 0///her2: 0///size: NA///node: 1///age: 71///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_281///dataset: MAQC2///series: MAQC2///id: 281///er: 1///pgr: 1///her2: 0///size: NA///node: 0///age: 45///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_282///dataset: MAQC2///series: MAQC2///id: 282///er: 1///pgr: 0///her2: 0///size: NA///node: 1///age: 44///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_283///dataset: MAQC2///series: MAQC2///id: 283///er: 1///pgr: 0///her2: 0///size: NA///node: 1///age: 46///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_286///dataset: MAQC2///series: MAQC2///id: 286///er: 0///pgr: 0///her2: 1///size: NA///node: 1///age: 51///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_287///dataset: MAQC2///series: MAQC2///id: 287///er: 0///pgr: 0///her2: 1///size: NA///node: 1///age: 55///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_294///dataset: MAQC2///series: MAQC2///id: 294///er: 0///pgr: 0///her2: 1///size: NA///node: 1///age: 53///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_295///dataset: MAQC2///series: MAQC2///id: 295///er: 0///pgr: 0///her2: 0///size: NA///node: 1///age: 48///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_297///dataset: MAQC2///series: MAQC2///id: 297///er: 1///pgr: 0///her2: 0///size: NA///node: 1///age: 65///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_301///dataset: MAQC2///series: MAQC2///id: 301///er: 0///pgr: 0///her2: 0///size: NA///node: 1///age: 46///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_302///dataset: MAQC2///series: MAQC2///id: 302///er: 0///pgr: 1///her2: 1///size: NA///node: 1///age: 50///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_304///dataset: MAQC2///series: MAQC2///id: 304///er: 1///pgr: 1///her2: 0///size: NA///node: 0///age: 46///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_309///dataset: MAQC2///series: MAQC2///id: 309///er: 0///pgr: 0///her2: 1///size: NA///node: 1///age: 62///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_310///dataset: MAQC2///series: MAQC2///id: 310///er: 0///pgr: 0///her2: 0///size: NA///node: 1///age: 51///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_315///dataset: MAQC2///series: MAQC2///id: 315///er: 1///pgr: 1///her2: 0///size: NA///node: 1///age: 36///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_316///dataset: MAQC2///series: MAQC2///id: 316///er: 0///pgr: 0///her2: 1///size: NA///node: 1///age: 45///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_319///dataset: MAQC2///series: MAQC2///id: 319///er: 1///pgr: 1///her2: 0///size: NA///node: 0///age: 26///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_322///dataset: MAQC2///series: MAQC2///id: 322///er: 1///pgr: 0///her2: 1///size: NA///node: 0///age: 65///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_323///dataset: MAQC2///series: MAQC2///id: 323///er: 1///pgr: 1///her2: 0///size: NA///node: 1///age: 40///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
samplename: MAQC2_330///dataset: MAQC2///series: MAQC2///id: 330///er: 1///pgr: 1///her2: 0///size: NA///node: 1///age: 35///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                       1 
                                                                                                                                                                                                                                 (Other) 
                                                                                                                                                                                                                                     131 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16716 }
\keyword{datasets}
