\name{ TRANSBIG_eset }
\alias{ TRANSBIG_eset }
\docType{data}
\title{ TRANSBIG_eset }
\usage{data( TRANSBIG_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=20813035 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse18229 
  PMIDs: 20813035 
  No abstract available.
  notes:
   summary:     
      Clinically, this subtype is usually triple negative invasive ductal carcinomas with a poor prognosis. Response to standar
d of care preoperative chemotherapy is intermediate between basal-like and luminal tumors. The claudin-low subtype most c
losely resembles the mammary epithelial stem cell.
   version:     
      2015-04-27 19:46:28
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
assayData: 20967 features, 198 samples
Platform type: 
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

records   n.max n.start  events  median 0.95LCL 0.95UCL 
  198.0   198.0   198.0    56.0      NA    17.1      NA 

--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
  VDXGUYU_4002   VDXGUYU_4008   VDXGUYU_4011   VDXGUYU_4014   VDXGUYU_4022   VDXGUYU_4033   VDXGUYU_4034   VDXGUYU_4036 
             1              1              1              1              1              1              1              1 
  VDXGUYU_4040   VDXGUYU_4041   VDXGUYU_4043   VDXGUYU_4044   VDXGUYU_4045   VDXGUYU_4046   VDXGUYU_4049   VDXGUYU_4052 
             1              1              1              1              1              1              1              1 
  VDXGUYU_4054   VDXGUYU_4056   VDXGUYU_4060   VDXGUYU_4061   VDXGUYU_4062   VDXGUYU_4067   VDXGUYU_4072   VDXGUYU_4073 
             1              1              1              1              1              1              1              1 
  VDXGUYU_4074   VDXGUYU_4075   VDXGUYU_4076   VDXGUYU_4080   VDXGUYU_4081   VDXGUYU_4087   VDXGUYU_4088   VDXGUYU_4096 
             1              1              1              1              1              1              1              1 
  VDXGUYU_4097   VDXGUYU_4098   VDXGUYU_4099   VDXGUYU_4105 VDXIGRU_166154 VDXIGRU_171150 VDXIGRU_171260 VDXIGRU_171558 
             1              1              1              1              1              1              1              1 
VDXIGRU_171620 VDXIGRU_172086 VDXIGRU_211980 VDXIGRU_214603 VDXIGRU_216806 VDXIGRU_216939 VDXIGRU_217320 VDXIGRU_219139 
             1              1              1              1              1              1              1              1 
VDXIGRU_219468 VDXIGRU_219490 VDXIGRU_219587 VDXIGRU_219863 VDXIGRU_219873 VDXIGRU_245392 VDXIGRU_245596 VDXIGRU_246047 
             1              1              1              1              1              1              1              1 
VDXIGRU_246166 VDXIGRU_246381 VDXIGRU_246560 VDXIGRU_246821 VDXIGRU_246919 VDXIGRU_247000 VDXIGRU_247243 VDXIGRU_247401 
             1              1              1              1              1              1              1              1 
VDXIGRU_254031 VDXIGRU_254485 VDXIGRU_255160 VDXIGRU_255627 VDXIGRU_256209 VDXIGRU_256734 VDXIGRU_256940 VDXIGRU_256975 
             1              1              1              1              1              1              1              1 
VDXIGRU_271448 VDXIGRU_271982 VDXIGRU_272219 VDXIGRU_272823 VDXIGRU_280328 VDXIGRU_280352 VDXIGRU_280775 VDXIGRU_287897 
             1              1              1              1              1              1              1              1 
VDXIGRU_306818 VDXIGRU_307251 VDXIGRU_307327 VDXIGRU_307506 VDXIGRU_307829 VDXIGRU_314827    VDXKIU_1000  VDXKIU_101B88 
             1              1              1              1              1              1              1              1 
    VDXKIU_121    VDXKIU_1246  VDXKIU_125B43  VDXKIU_136B04  VDXKIU_144B49    VDXKIU_15C5    VDXKIU_15D6    VDXKIU_15D7 
             1              1              1              1              1              1              1              1 
   VDXKIU_15D9    VDXKIU_15E2    VDXKIU_15E6        (Other) 
             1              1              1             99 

sample_type: 
tumor 
  198 

er: 
negative positive 
      64      134 

tumor_size: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  0.600   1.700   2.000   2.181   2.500   5.000 

N: 
  0 
198 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  24.00   42.00   46.00   46.39   51.00   60.00 

grade: 
   1    2    3 NA's 
  30   83   83    2 

dmfs_days: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    125    2375    4384    3954    5566    9108 

dmfs_status: 
deceased_or_recurrence    living_norecurrence 
                    62                    136 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    121    1528    3534    3399    5130    8711 

recurrence_status: 
deceased_or_recurrence    living_norecurrence 
                    91                    107 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    146    2744    4562    4150    5610    9108 

vital_status: 
deceased   living 
      56      142 

treatment: 
untreated 
      198 

batch: 
VDXGUYU VDXIGRU  VDXKIU VDXOXFU  VDXRHU 
     36      50      51      24      37 

uncurated_author_metadata: 
         samplename: VDXGUYU_4002///dataset: TRANSBIG///series: VDXGUYU///id: 4002///er: 0///pgr: NA///her2: NA///size: 3///node: 0///age: 57///grade: 3///t.dmfs: 723///e.dmfs: 1///t.rfs: 723///e.rfs: 1///t.os: 937///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
       samplename: VDXGUYU_4008///dataset: TRANSBIG///series: VDXGUYU///id: 4008///er: 1///pgr: NA///her2: NA///size: 3///node: 0///age: 57///grade: 3///t.dmfs: 6591///e.dmfs: 0///t.rfs: 183///e.rfs: 1///t.os: 6591///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
       samplename: VDXGUYU_4011///dataset: TRANSBIG///series: VDXGUYU///id: 4011///er: 0///pgr: NA///her2: NA///size: 2.5///node: 0///age: 48///grade: 3///t.dmfs: 524///e.dmfs: 1///t.rfs: 524///e.rfs: 1///t.os: 922///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
    samplename: VDXGUYU_4014///dataset: TRANSBIG///series: VDXGUYU///id: 4014///er: 1///pgr: NA///her2: NA///size: 1.8///node: 0///age: 42///grade: 3///t.dmfs: 6255///e.dmfs: 1///t.rfs: 2192///e.rfs: 1///t.os: 6255///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
      samplename: VDXGUYU_4022///dataset: TRANSBIG///series: VDXGUYU///id: 4022///er: 1///pgr: NA///her2: NA///size: 3///node: 0///age: 46///grade: 2///t.dmfs: 3822///e.dmfs: 1///t.rfs: 3822///e.rfs: 1///t.os: 4133///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
      samplename: VDXGUYU_4033///dataset: TRANSBIG///series: VDXGUYU///id: 4033///er: 1///pgr: NA///her2: NA///size: 2///node: 0///age: 58///grade: 2///t.dmfs: 6507///e.dmfs: 0///t.rfs: 6507///e.rfs: 0///t.os: 6507///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
       samplename: VDXGUYU_4034///dataset: TRANSBIG///series: VDXGUYU///id: 4034///er: 0///pgr: NA///her2: NA///size: 2///node: 0///age: 44///grade: 3///t.dmfs: 5947///e.dmfs: 0///t.rfs: 709///e.rfs: 1///t.os: 5947///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
    samplename: VDXGUYU_4036///dataset: TRANSBIG///series: VDXGUYU///id: 4036///er: 1///pgr: NA///her2: NA///size: 2.5///node: 0///age: 58///grade: 1///t.dmfs: 5816///e.dmfs: 1///t.rfs: 5816///e.rfs: 1///t.os: 5816///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
      samplename: VDXGUYU_4040///dataset: TRANSBIG///series: VDXGUYU///id: 4040///er: 1///pgr: NA///her2: NA///size: 3///node: 0///age: 47///grade: 3///t.dmfs: 6007///e.dmfs: 0///t.rfs: 6007///e.rfs: 0///t.os: 6007///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
     samplename: VDXGUYU_4041///dataset: TRANSBIG///series: VDXGUYU///id: 4041///er: 1///pgr: NA///her2: NA///size: 2.5///node: 0///age: 38///grade: 2///t.dmfs: 1233///e.dmfs: 1///t.rfs: 422///e.rfs: 1///t.os: 1484///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
      samplename: VDXGUYU_4043///dataset: TRANSBIG///series: VDXGUYU///id: 4043///er: 1///pgr: NA///her2: NA///size: 2///node: 0///age: 59///grade: 3///t.dmfs: 1136///e.dmfs: 1///t.rfs: 1136///e.rfs: 1///t.os: 1911///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
    samplename: VDXGUYU_4044///dataset: TRANSBIG///series: VDXGUYU///id: 4044///er: 1///pgr: NA///her2: NA///size: 0.8///node: 0///age: 50///grade: 1///t.dmfs: 5984///e.dmfs: 0///t.rfs: 5984///e.rfs: 0///t.os: 5984///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
    samplename: VDXGUYU_4045///dataset: TRANSBIG///series: VDXGUYU///id: 4045///er: 0///pgr: NA///her2: NA///size: 4.5///node: 0///age: 46///grade: 3///t.dmfs: 6249///e.dmfs: 0///t.rfs: 6249///e.rfs: 0///t.os: 6249///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
      samplename: VDXGUYU_4046///dataset: TRANSBIG///series: VDXGUYU///id: 4046///er: 0///pgr: NA///her2: NA///size: 2///node: 0///age: 57///grade: 2///t.dmfs: 4066///e.dmfs: 1///t.rfs: 2841///e.rfs: 1///t.os: 4372///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
       samplename: VDXGUYU_4049///dataset: TRANSBIG///series: VDXGUYU///id: 4049///er: 0///pgr: NA///her2: NA///size: 2.5///node: 0///age: 33///grade: 3///t.dmfs: 690///e.dmfs: 1///t.rfs: 458///e.rfs: 1///t.os: 977///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
        samplename: VDXGUYU_4052///dataset: TRANSBIG///series: VDXGUYU///id: 4052///er: 1///pgr: NA///her2: NA///size: 2///node: 0///age: 42///grade: 3///t.dmfs: 994///e.dmfs: 1///t.rfs: 994///e.rfs: 1///t.os: 3830///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
         samplename: VDXGUYU_4054///dataset: TRANSBIG///series: VDXGUYU///id: 4054///er: 0///pgr: NA///her2: NA///size: 2///node: 0///age: 38///grade: 3///t.dmfs: 434///e.dmfs: 1///t.rfs: 282///e.rfs: 1///t.os: 667///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
     samplename: VDXGUYU_4056///dataset: TRANSBIG///series: VDXGUYU///id: 4056///er: 0///pgr: NA///her2: NA///size: 2.5///node: 0///age: 55///grade: 2///t.dmfs: 1743///e.dmfs: 1///t.rfs: 548///e.rfs: 1///t.os: 2129///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
    samplename: VDXGUYU_4060///dataset: TRANSBIG///series: VDXGUYU///id: 4060///er: 0///pgr: NA///her2: NA///size: 2.5///node: 0///age: 57///grade: 2///t.dmfs: 2863///e.dmfs: 1///t.rfs: 1820///e.rfs: 1///t.os: 3590///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
    samplename: VDXGUYU_4061///dataset: TRANSBIG///series: VDXGUYU///id: 4061///er: 1///pgr: NA///her2: NA///size: 2.5///node: 0///age: 46///grade: 2///t.dmfs: 5693///e.dmfs: 0///t.rfs: 5693///e.rfs: 0///t.os: 5693///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
      samplename: VDXGUYU_4062///dataset: TRANSBIG///series: VDXGUYU///id: 4062///er: 1///pgr: NA///her2: NA///size: 2///node: 0///age: 50///grade: 2///t.dmfs: 5823///e.dmfs: 0///t.rfs: 5823///e.rfs: 0///t.os: 5823///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
      samplename: VDXGUYU_4067///dataset: TRANSBIG///series: VDXGUYU///id: 4067///er: 0///pgr: NA///her2: NA///size: 2///node: 0///age: 46///grade: 3///t.dmfs: 5785///e.dmfs: 0///t.rfs: 5785///e.rfs: 0///t.os: 5785///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
    samplename: VDXGUYU_4072///dataset: TRANSBIG///series: VDXGUYU///id: 4072///er: 1///pgr: NA///her2: NA///size: 2.5///node: 0///age: 46///grade: 2///t.dmfs: 5661///e.dmfs: 0///t.rfs: 5598///e.rfs: 1///t.os: 5661///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
    samplename: VDXGUYU_4073///dataset: TRANSBIG///series: VDXGUYU///id: 4073///er: 1///pgr: NA///her2: NA///size: 1.5///node: 0///age: 51///grade: 2///t.dmfs: 2909///e.dmfs: 1///t.rfs: 2909///e.rfs: 1///t.os: 2909///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
     samplename: VDXGUYU_4074///dataset: TRANSBIG///series: VDXGUYU///id: 4074///er: 0///pgr: NA///her2: NA///size: 2.5///node: 0///age: 58///grade: 3///t.dmfs: 1171///e.dmfs: 1///t.rfs: 121///e.rfs: 1///t.os: 1674///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
    samplename: VDXGUYU_4075///dataset: TRANSBIG///series: VDXGUYU///id: 4075///er: 1///pgr: NA///her2: NA///size: 2.6///node: 0///age: 59///grade: 1///t.dmfs: 5110///e.dmfs: 0///t.rfs: 5110///e.rfs: 0///t.os: 5110///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
    samplename: VDXGUYU_4076///dataset: TRANSBIG///series: VDXGUYU///id: 4076///er: 1///pgr: NA///her2: NA///size: 0.6///node: 0///age: 49///grade: 2///t.dmfs: 5533///e.dmfs: 0///t.rfs: 5533///e.rfs: 0///t.os: 5533///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
    samplename: VDXGUYU_4080///dataset: TRANSBIG///series: VDXGUYU///id: 4080///er: 1///pgr: NA///her2: NA///size: 2.5///node: 0///age: 43///grade: 1///t.dmfs: 5571///e.dmfs: 0///t.rfs: 1429///e.rfs: 1///t.os: 5571///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
      samplename: VDXGUYU_4081///dataset: TRANSBIG///series: VDXGUYU///id: 4081///er: 0///pgr: NA///her2: NA///size: 2///node: 0///age: 47///grade: 3///t.dmfs: 5519///e.dmfs: 0///t.rfs: 5519///e.rfs: 0///t.os: 5519///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
    samplename: VDXGUYU_4087///dataset: TRANSBIG///series: VDXGUYU///id: 4087///er: 1///pgr: NA///her2: NA///size: 1.7///node: 0///age: 37///grade: 3///t.dmfs: 5379///e.dmfs: 0///t.rfs: 5379///e.rfs: 0///t.os: 5379///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
      samplename: VDXGUYU_4088///dataset: TRANSBIG///series: VDXGUYU///id: 4088///er: 1///pgr: NA///her2: NA///size: 4///node: 0///age: 48///grade: 3///t.dmfs: 4675///e.dmfs: 0///t.rfs: 4675///e.rfs: 0///t.os: 4675///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
    samplename: VDXGUYU_4096///dataset: TRANSBIG///series: VDXGUYU///id: 4096///er: 1///pgr: NA///her2: NA///size: 1.9///node: 0///age: 47///grade: 2///t.dmfs: 4652///e.dmfs: 0///t.rfs: 1666///e.rfs: 1///t.os: 4652///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
    samplename: VDXGUYU_4097///dataset: TRANSBIG///series: VDXGUYU///id: 4097///er: 1///pgr: NA///her2: NA///size: 2.5///node: 0///age: 46///grade: 2///t.dmfs: 4227///e.dmfs: 0///t.rfs: 4227///e.rfs: 0///t.os: 4227///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
        samplename: VDXGUYU_4098///dataset: TRANSBIG///series: VDXGUYU///id: 4098///er: 0///pgr: NA///her2: NA///size: 2///node: 0///age: 48///grade: 3///t.dmfs: 910///e.dmfs: 1///t.rfs: 316///e.rfs: 1///t.os: 1530///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
    samplename: VDXGUYU_4099///dataset: TRANSBIG///series: VDXGUYU///id: 4099///er: 1///pgr: NA///her2: NA///size: 1.9///node: 0///age: 58///grade: 2///t.dmfs: 4472///e.dmfs: 0///t.rfs: 4472///e.rfs: 0///t.os: 4472///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
       samplename: VDXGUYU_4105///dataset: TRANSBIG///series: VDXGUYU///id: 4105///er: 1///pgr: NA///her2: NA///size: 2///node: 0///age: 34///grade: 3///t.dmfs: 1206///e.dmfs: 1///t.rfs: 350///e.rfs: 1///t.os: 1275///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
samplename: VDXIGRU_166154///dataset: TRANSBIG///series: VDXIGRU///id: 166154///er: 1///pgr: NA///her2: NA///size: 4.5///node: 0///age: 30///grade: 2///t.dmfs: 7057///e.dmfs: 0///t.rfs: 5632///e.rfs: 1///t.os: 7057///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
     samplename: VDXIGRU_171150///dataset: TRANSBIG///series: VDXIGRU///id: 171150///er: 1///pgr: NA///her2: NA///size: 3///node: 0///age: 35///grade: 2///t.dmfs: 958///e.dmfs: 0///t.rfs: 958///e.rfs: 0///t.os: 958///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
samplename: VDXIGRU_171260///dataset: TRANSBIG///series: VDXIGRU///id: 171260///er: 0///pgr: NA///her2: NA///size: 2.6///node: 0///age: 56///grade: 3///t.dmfs: 5359///e.dmfs: 0///t.rfs: 5359///e.rfs: 0///t.os: 5359///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
samplename: VDXIGRU_171558///dataset: TRANSBIG///series: VDXIGRU///id: 171558///er: 0///pgr: NA///her2: NA///size: 3.2///node: 0///age: 35///grade: 2///t.dmfs: 5945///e.dmfs: 0///t.rfs: 5945///e.rfs: 0///t.os: 5945///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
samplename: VDXIGRU_171620///dataset: TRANSBIG///series: VDXIGRU///id: 171620///er: 0///pgr: NA///her2: NA///size: 2.7///node: 0///age: 52///grade: 3///t.dmfs: 5651///e.dmfs: 0///t.rfs: 5651///e.rfs: 0///t.os: 5651///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
    samplename: VDXIGRU_172086///dataset: TRANSBIG///series: VDXIGRU///id: 172086///er: 1///pgr: NA///her2: NA///size: 3///node: 0///age: 47///grade: 3///t.dmfs: 528///e.dmfs: 1///t.rfs: 528///e.rfs: 1///t.os: 6036///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
  samplename: VDXIGRU_211980///dataset: TRANSBIG///series: VDXIGRU///id: 211980///er: 0///pgr: NA///her2: NA///size: 1.6///node: 0///age: 42///grade: 3///t.dmfs: 796///e.dmfs: 1///t.rfs: 408///e.rfs: 1///t.os: 1436///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
samplename: VDXIGRU_214603///dataset: TRANSBIG///series: VDXIGRU///id: 214603///er: 1///pgr: NA///her2: NA///size: 1.7///node: 0///age: 54///grade: 2///t.dmfs: 5685///e.dmfs: 0///t.rfs: 5685///e.rfs: 0///t.os: 5685///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
samplename: VDXIGRU_216806///dataset: TRANSBIG///series: VDXIGRU///id: 216806///er: 1///pgr: NA///her2: NA///size: 2.5///node: 0///age: 42///grade: 3///t.dmfs: 1023///e.dmfs: 1///t.rfs: 1023///e.rfs: 1///t.os: 1754///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
samplename: VDXIGRU_216939///dataset: TRANSBIG///series: VDXIGRU///id: 216939///er: 1///pgr: NA///her2: NA///size: 2.5///node: 0///age: 50///grade: 2///t.dmfs: 4791///e.dmfs: 0///t.rfs: 4791///e.rfs: 0///t.os: 4791///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
  samplename: VDXIGRU_217320///dataset: TRANSBIG///series: VDXIGRU///id: 217320///er: 1///pgr: NA///her2: NA///size: 3///node: 0///age: 42///grade: 1///t.dmfs: 1215///e.dmfs: 0///t.rfs: 1215///e.rfs: 0///t.os: 1215///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
samplename: VDXIGRU_219139///dataset: TRANSBIG///series: VDXIGRU///id: 219139///er: 1///pgr: NA///her2: NA///size: 3.5///node: 0///age: 38///grade: 1///t.dmfs: 3779///e.dmfs: 0///t.rfs: 3779///e.rfs: 0///t.os: 3779///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
samplename: VDXIGRU_219468///dataset: TRANSBIG///series: VDXIGRU///id: 219468///er: 1///pgr: NA///her2: NA///size: 2.5///node: 0///age: 37///grade: 3///t.dmfs: 5306///e.dmfs: 0///t.rfs: 5306///e.rfs: 0///t.os: 5306///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
  samplename: VDXIGRU_219490///dataset: TRANSBIG///series: VDXIGRU///id: 219490///er: 1///pgr: NA///her2: NA///size: 3///node: 0///age: 44///grade: 3///t.dmfs: 5678///e.dmfs: 0///t.rfs: 5678///e.rfs: 0///t.os: 5678///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
   samplename: VDXIGRU_219587///dataset: TRANSBIG///series: VDXIGRU///id: 219587///er: 0///pgr: NA///her2: NA///size: 3.1///node: 0///age: 37///grade: 3///t.dmfs: 421///e.dmfs: 1///t.rfs: 421///e.rfs: 1///t.os: 794///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
 samplename: VDXIGRU_219863///dataset: TRANSBIG///series: VDXIGRU///id: 219863///er: 1///pgr: NA///her2: NA///size: 1.5///node: 0///age: 51///grade: 2///t.dmfs: 5612///e.dmfs: 0///t.rfs: 840///e.rfs: 1///t.os: 5612///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
  samplename: VDXIGRU_219873///dataset: TRANSBIG///series: VDXIGRU///id: 219873///er: 0///pgr: NA///her2: NA///size: 4///node: 0///age: 53///grade: 3///t.dmfs: 3121///e.dmfs: 1///t.rfs: 3121///e.rfs: 1///t.os: 3326///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
samplename: VDXIGRU_245392///dataset: TRANSBIG///series: VDXIGRU///id: 245392///er: 1///pgr: NA///her2: NA///size: 2.5///node: 0///age: 49///grade: 1///t.dmfs: 1104///e.dmfs: 0///t.rfs: 1104///e.rfs: 0///t.os: 1104///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
  samplename: VDXIGRU_245596///dataset: TRANSBIG///series: VDXIGRU///id: 245596///er: 0///pgr: NA///her2: NA///size: 3///node: 0///age: 43///grade: 3///t.dmfs: 4691///e.dmfs: 0///t.rfs: 4691///e.rfs: 0///t.os: 4691///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
samplename: VDXIGRU_246047///dataset: TRANSBIG///series: VDXIGRU///id: 246047///er: 1///pgr: NA///her2: NA///size: 2.2///node: 0///age: 39///grade: 2///t.dmfs: 5238///e.dmfs: 0///t.rfs: 1609///e.rfs: 1///t.os: 5238///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
  samplename: VDXIGRU_246166///dataset: TRANSBIG///series: VDXIGRU///id: 246166///er: 1///pgr: NA///her2: NA///size: 2///node: 0///age: 47///grade: 2///t.dmfs: 5101///e.dmfs: 0///t.rfs: 5101///e.rfs: 0///t.os: 5101///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
samplename: VDXIGRU_246381///dataset: TRANSBIG///series: VDXIGRU///id: 246381///er: 1///pgr: NA///her2: NA///size: 1.9///node: 0///age: 42///grade: 3///t.dmfs: 2809///e.dmfs: 0///t.rfs: 2809///e.rfs: 0///t.os: 2809///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
samplename: VDXIGRU_246560///dataset: TRANSBIG///series: VDXIGRU///id: 246560///er: 0///pgr: NA///her2: NA///size: 1.2///node: 0///age: 39///grade: 2///t.dmfs: 3296///e.dmfs: 0///t.rfs: 3296///e.rfs: 0///t.os: 3296///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
   samplename: VDXIGRU_246821///dataset: TRANSBIG///series: VDXIGRU///id: 246821///er: 1///pgr: NA///her2: NA///size: 3///node: 0///age: 52///grade: 2///t.dmfs: 3555///e.dmfs: 0///t.rfs: 546///e.rfs: 1///t.os: 3555///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
samplename: VDXIGRU_246919///dataset: TRANSBIG///series: VDXIGRU///id: 246919///er: 1///pgr: NA///her2: NA///size: 1.5///node: 0///age: 50///grade: 2///t.dmfs: 5433///e.dmfs: 0///t.rfs: 3525///e.rfs: 1///t.os: 5433///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
 samplename: VDXIGRU_247000///dataset: TRANSBIG///series: VDXIGRU///id: 247000///er: 1///pgr: NA///her2: NA///size: 2.2///node: 0///age: 36///grade: 3///t.dmfs: 5228///e.dmfs: 0///t.rfs: 854///e.rfs: 1///t.os: 5228///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
samplename: VDXIGRU_247243///dataset: TRANSBIG///series: VDXIGRU///id: 247243///er: 1///pgr: NA///her2: NA///size: 2.5///node: 0///age: 45///grade: 2///t.dmfs: 5316///e.dmfs: 0///t.rfs: 1813///e.rfs: 1///t.os: 5316///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
samplename: VDXIGRU_247401///dataset: TRANSBIG///series: VDXIGRU///id: 247401///er: 1///pgr: NA///her2: NA///size: 1.7///node: 0///age: 50///grade: 1///t.dmfs: 4925///e.dmfs: 0///t.rfs: 4925///e.rfs: 0///t.os: 4925///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
samplename: VDXIGRU_254031///dataset: TRANSBIG///series: VDXIGRU///id: 254031///er: 1///pgr: NA///her2: NA///size: 2.3///node: 0///age: 39///grade: 3///t.dmfs: 5012///e.dmfs: 0///t.rfs: 3853///e.rfs: 1///t.os: 5012///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
  samplename: VDXIGRU_254485///dataset: TRANSBIG///series: VDXIGRU///id: 254485///er: 1///pgr: NA///her2: NA///size: 3///node: 0///age: 44///grade: 1///t.dmfs: 4642///e.dmfs: 0///t.rfs: 4642///e.rfs: 0///t.os: 4642///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
samplename: VDXIGRU_255160///dataset: TRANSBIG///series: VDXIGRU///id: 255160///er: 1///pgr: NA///her2: NA///size: 1.3///node: 0///age: 50///grade: 1///t.dmfs: 3311///e.dmfs: 0///t.rfs: 3311///e.rfs: 0///t.os: 3311///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
samplename: VDXIGRU_255627///dataset: TRANSBIG///series: VDXIGRU///id: 255627///er: 1///pgr: NA///her2: NA///size: 1.8///node: 0///age: 51///grade: 2///t.dmfs: 4863///e.dmfs: 0///t.rfs: 4863///e.rfs: 0///t.os: 4863///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
  samplename: VDXIGRU_256209///dataset: TRANSBIG///series: VDXIGRU///id: 256209///er: 1///pgr: NA///her2: NA///size: 2///node: 0///age: 49///grade: 1///t.dmfs: 1815///e.dmfs: 1///t.rfs: 1780///e.rfs: 1///t.os: 3347///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
    samplename: VDXIGRU_256734///dataset: TRANSBIG///series: VDXIGRU///id: 256734///er: 0///pgr: NA///her2: NA///size: 2///node: 0///age: 37///grade: 2///t.dmfs: 803///e.dmfs: 1///t.rfs: 803///e.rfs: 1///t.os: 1049///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
samplename: VDXIGRU_256940///dataset: TRANSBIG///series: VDXIGRU///id: 256940///er: 1///pgr: NA///her2: NA///size: 1.8///node: 0///age: 51///grade: 2///t.dmfs: 2105///e.dmfs: 0///t.rfs: 2105///e.rfs: 0///t.os: 2105///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
     samplename: VDXIGRU_256975///dataset: TRANSBIG///series: VDXIGRU///id: 256975///er: 0///pgr: NA///her2: NA///size: 2///node: 0///age: 40///grade: 3///t.dmfs: 530///e.dmfs: 1///t.rfs: 530///e.rfs: 1///t.os: 929///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
  samplename: VDXIGRU_271448///dataset: TRANSBIG///series: VDXIGRU///id: 271448///er: 1///pgr: NA///her2: NA///size: 2///node: 0///age: 44///grade: 1///t.dmfs: 4830///e.dmfs: 0///t.rfs: 4830///e.rfs: 0///t.os: 4830///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
  samplename: VDXIGRU_271982///dataset: TRANSBIG///series: VDXIGRU///id: 271982///er: 1///pgr: NA///her2: NA///size: 2///node: 0///age: 49///grade: 3///t.dmfs: 1721///e.dmfs: 1///t.rfs: 1721///e.rfs: 1///t.os: 2148///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
samplename: VDXIGRU_272219///dataset: TRANSBIG///series: VDXIGRU///id: 272219///er: 1///pgr: NA///her2: NA///size: 1.5///node: 0///age: 53///grade: 1///t.dmfs: 4720///e.dmfs: 0///t.rfs: 4720///e.rfs: 0///t.os: 4720///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
samplename: VDXIGRU_272823///dataset: TRANSBIG///series: VDXIGRU///id: 272823///er: 1///pgr: NA///her2: NA///size: 2.9///node: 0///age: 44///grade: 2///t.dmfs: 1259///e.dmfs: 1///t.rfs: 1259///e.rfs: 1///t.os: 2453///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
samplename: VDXIGRU_280328///dataset: TRANSBIG///series: VDXIGRU///id: 280328///er: 1///pgr: NA///her2: NA///size: 2.3///node: 0///age: 48///grade: 2///t.dmfs: 3976///e.dmfs: 1///t.rfs: 3976///e.rfs: 1///t.os: 4438///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
samplename: VDXIGRU_280352///dataset: TRANSBIG///series: VDXIGRU///id: 280352///er: 1///pgr: NA///her2: NA///size: 2.5///node: 0///age: 57///grade: 2///t.dmfs: 4050///e.dmfs: 1///t.rfs: 3331///e.rfs: 1///t.os: 4469///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
  samplename: VDXIGRU_280775///dataset: TRANSBIG///series: VDXIGRU///id: 280775///er: 1///pgr: NA///her2: NA///size: 3///node: 0///age: 45///grade: 3///t.dmfs: 1740///e.dmfs: 0///t.rfs: 1740///e.rfs: 0///t.os: 1740///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
  samplename: VDXIGRU_287897///dataset: TRANSBIG///series: VDXIGRU///id: 287897///er: 1///pgr: NA///her2: NA///size: 2///node: 0///age: 44///grade: 1///t.dmfs: 3608///e.dmfs: 1///t.rfs: 3218///e.rfs: 1///t.os: 3608///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
samplename: VDXIGRU_306818///dataset: TRANSBIG///series: VDXIGRU///id: 306818///er: 1///pgr: NA///her2: NA///size: 1.5///node: 0///age: 39///grade: 2///t.dmfs: 3900///e.dmfs: 0///t.rfs: 2536///e.rfs: 1///t.os: 3900///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
  samplename: VDXIGRU_307251///dataset: TRANSBIG///series: VDXIGRU///id: 307251///er: 0///pgr: NA///her2: NA///size: 3///node: 0///age: 40///grade: 3///t.dmfs: 1554///e.dmfs: 0///t.rfs: 1554///e.rfs: 0///t.os: 1554///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
  samplename: VDXIGRU_307327///dataset: TRANSBIG///series: VDXIGRU///id: 307327///er: 1///pgr: NA///her2: NA///size: 2///node: 0///age: 38///grade: 3///t.dmfs: 3544///e.dmfs: 0///t.rfs: 3544///e.rfs: 0///t.os: 3544///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
samplename: VDXIGRU_307506///dataset: TRANSBIG///series: VDXIGRU///id: 307506///er: 1///pgr: NA///her2: NA///size: 2.6///node: 0///age: 45///grade: 1///t.dmfs: 3721///e.dmfs: 0///t.rfs: 3721///e.rfs: 0///t.os: 3721///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
samplename: VDXIGRU_307829///dataset: TRANSBIG///series: VDXIGRU///id: 307829///er: 1///pgr: NA///her2: NA///size: 2.6///node: 0///age: 43///grade: 2///t.dmfs: 3771///e.dmfs: 0///t.rfs: 3771///e.rfs: 0///t.os: 3771///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
samplename: VDXIGRU_314827///dataset: TRANSBIG///series: VDXIGRU///id: 314827///er: 0///pgr: NA///her2: NA///size: 1.5///node: 0///age: 47///grade: 3///t.dmfs: 3293///e.dmfs: 0///t.rfs: 3293///e.rfs: 0///t.os: 3293///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
           samplename: VDXKIU_1000///dataset: TRANSBIG///series: VDXKIU///id: 1000///er: 0///pgr: NA///her2: NA///size: 4///node: 0///age: 43///grade: 2///t.dmfs: 394///e.dmfs: 1///t.rfs: 394///e.rfs: 1///t.os: 491///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
  samplename: VDXKIU_101B88///dataset: TRANSBIG///series: VDXKIU///id: 101B88///er: 0///pgr: NA///her2: NA///size: 1.2///node: 0///age: 40///grade: 3///t.dmfs: 6240///e.dmfs: 0///t.rfs: 6240///e.rfs: 0///t.os: 6240///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
          samplename: VDXKIU_121///dataset: TRANSBIG///series: VDXKIU///id: 121///er: 1///pgr: NA///her2: NA///size: 1///node: 0///age: 54///grade: 2///t.dmfs: 4568///e.dmfs: 0///t.rfs: 4568///e.rfs: 0///t.os: 4568///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
      samplename: VDXKIU_1246///dataset: TRANSBIG///series: VDXKIU///id: 1246///er: 1///pgr: NA///her2: NA///size: 0.9///node: 0///age: 24///grade: 2///t.dmfs: 6777///e.dmfs: 1///t.rfs: 3512///e.rfs: 1///t.os: 8251///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
  samplename: VDXKIU_125B43///dataset: TRANSBIG///series: VDXKIU///id: 125B43///er: 1///pgr: NA///her2: NA///size: 1.5///node: 0///age: 38///grade: 2///t.dmfs: 6161///e.dmfs: 0///t.rfs: 6161///e.rfs: 0///t.os: 6161///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
   samplename: VDXKIU_136B04///dataset: TRANSBIG///series: VDXKIU///id: 136B04///er: 1///pgr: NA///her2: NA///size: 1.7///node: 0///age: 54///grade: 2///t.dmfs: 6132///e.dmfs: 0///t.rfs: 892///e.rfs: 1///t.os: 6132///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
  samplename: VDXKIU_144B49///dataset: TRANSBIG///series: VDXKIU///id: 144B49///er: 1///pgr: NA///her2: NA///size: 2.1///node: 0///age: 40///grade: 2///t.dmfs: 4212///e.dmfs: 1///t.rfs: 4212///e.rfs: 1///t.os: 5337///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
        samplename: VDXKIU_15C5///dataset: TRANSBIG///series: VDXKIU///id: 15C5///er: 0///pgr: NA///her2: NA///size: 1///node: 0///age: 56///grade: 3///t.dmfs: 4980///e.dmfs: 0///t.rfs: 4980///e.rfs: 0///t.os: 4980///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
        samplename: VDXKIU_15D6///dataset: TRANSBIG///series: VDXKIU///id: 15D6///er: 1///pgr: NA///her2: NA///size: 2///node: 0///age: 45///grade: 2///t.dmfs: 1703///e.dmfs: 1///t.rfs: 1703///e.rfs: 1///t.os: 1703///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
        samplename: VDXKIU_15D7///dataset: TRANSBIG///series: VDXKIU///id: 15D7///er: 1///pgr: NA///her2: NA///size: 2///node: 0///age: 46///grade: 2///t.dmfs: 5014///e.dmfs: 0///t.rfs: 5014///e.rfs: 0///t.os: 5014///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
        samplename: VDXKIU_15D9///dataset: TRANSBIG///series: VDXKIU///id: 15D9///er: 1///pgr: NA///her2: NA///size: 1.6///node: 0///age: 45///grade: 3///t.dmfs: 880///e.dmfs: 1///t.rfs: 880///e.rfs: 1///t.os: 1924///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
  samplename: VDXKIU_15E2///dataset: TRANSBIG///series: VDXKIU///id: 1.50E+03///er: 1///pgr: NA///her2: NA///size: 1.4///node: 0///age: 40///grade: 3///t.dmfs: 2549///e.dmfs: 0///t.rfs: 2549///e.rfs: 0///t.os: 2549///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
  samplename: VDXKIU_15E6///dataset: TRANSBIG///series: VDXKIU///id: 1.50E+07///er: 0///pgr: NA///her2: NA///size: 1.1///node: 0///age: 39///grade: 3///t.dmfs: 6298///e.dmfs: 0///t.rfs: 6298///e.rfs: 0///t.os: 6298///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                          1 
                                                                                                                                                                                                                                                    (Other) 
                                                                                                                                                                                                                                                         99 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse18229 }
\keyword{datasets}
