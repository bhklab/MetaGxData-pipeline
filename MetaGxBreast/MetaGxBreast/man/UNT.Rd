\name{ UNT_eset }
\alias{ UNT_eset }
\docType{data}
\title{ UNT_eset }
\usage{data( UNT_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information:  http://www.ncbi.nlm.nih.gov/pubmed/?term=15721472; http://www.ncbi.nlm.nih.gov/pubmed/?term=17420468 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse2034; http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse5327 
  PMIDs: 15721472, 17420468 
  No abstract available.
  notes:
   summary:     
      15721472: A gene signature was identified that can accurately predict distant metastases in node-negative cases. 17420468
: Tumors with a lung metastatic gene signature were shown to be larger.  
   version:     
      2015-04-27 19:50:26
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: 1007_s_at 1053_at ... 244889_at (36084 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 36084 features, 133 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
KIU_101B88 KIU_105B13 KIU_106B55 KIU_111B51 KIU_113B11 KIU_120B73 KIU_124B25 KIU_127B00 KIU_134B33 KIU_136B04 KIU_140B91 
         1          1          1          1          1          1          1          1          1          1          1 
KIU_144B49 KIU_151B84 KIU_155B52 KIU_163B27 KIU_164B81 KIU_172B19 KIU_177B67 KIU_184B38 KIU_188B13 KIU_196B81 KIU_197B95 
         1          1          1          1          1          1          1          1          1          1          1 
KIU_199B55 KIU_205B99 KIU_220C70 KIU_227C50 KIU_229C44 KIU_231C80 KIU_233C91 KIU_242C21 KIU_243C70 KIU_247C76 KIU_248C91 
         1          1          1          1          1          1          1          1          1          1          1 
 KIU_24C30 KIU_259C74 KIU_260C91 KIU_266C51 KIU_268C87 KIU_272C88 KIU_278C80 KIU_279C61 KIU_280C43 KIU_282C51 KIU_284C63 
         1          1          1          1          1          1          1          1          1          1          1 
KIU_286C91  KIU_28C76 KIU_292C66 KIU_303C36 KIU_304C89 KIU_308C93 KIU_309C49 KIU_314B55 KIU_316C64  KIU_36C17  KIU_42C67 
         1          1          1          1          1          1          1          1          1          1          1 
 KIU_43C47  KIU_52A90   KIU_5B97  KIU_65A68  KIU_74A63  KIU_86A40  KIU_87A79  KIU_88A67  KIU_89A64   OXFU_104  OXFU_1065 
         1          1          1          1          1          1          1          1          1          1          1 
 OXFU_1183    OXFU_12  OXFU_1210  OXFU_1248   OXFU_126   OXFU_127  OXFU_1286  OXFU_1328  OXFU_1373   OXFU_138  OXFU_1415 
         1          1          1          1          1          1          1          1          1          1          1 
  OXFU_145   OXFU_157    OXFU_16  OXFU_1605   OXFU_181   OXFU_217   OXFU_220   OXFU_223    OXFU_23   OXFU_245   OXFU_247 
         1          1          1          1          1          1          1          1          1          1          1 
  OXFU_254   OXFU_281   OXFU_316   OXFU_320   OXFU_337   OXFU_348   OXFU_359   OXFU_360   OXFU_366    OXFU_37   OXFU_373 
         1          1          1          1          1          1          1          1          1          1          1 
   (Other) 
        34 

alt_sample_name: 
 101B88     104  105B13    1065  106B55  111B51  113B11    1183      12  120B73    1210    1248  124B25     126     127 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
 127B00    1286    1328  134B33  136B04    1373     138  140B91    1415  144B49     145  151B84  155B52     157      16 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
   1605  163B27  164B81  172B19  177B67     181  184B38  188B13  196B81  197B95  199B55  205B99     217     220  220C70 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
    223  227C50  229C44      23  231C80  233C91  242C21  243C70     245     247  247C76  248C91   24C30     254  259C74 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
 260C91  266C51  268C87  272C88  278C80  279C61  280C43     281  282C51  284C63  286C91   28C76  292C66  303C36  304C89 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
 308C93  309C49  314B55     316  316C64     320     337     348     359     360     366   36C17      37     373     382 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
    397     419   42C67   43C47     446     449     471     476     484 (Other) 
      1       1       1       1       1       1       1       1       1      34 

sample_type: 
tumor 
  133 

er: 
negative positive     NA's 
      40       86        7 

pgr: 
negative positive     NA's 
       6       56       71 

tumor_size: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  0.000   1.200   1.900   1.892   2.300   6.000 

N: 
  0 
133 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  24.00   44.00   53.00   51.79   60.00   73.00 

grade: 
   1    2    3 NA's 
  32   51   29   21 

dmfs_days: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
     61    1338    2809    2724    4078    5305 

dmfs_status: 
deceased_or_recurrence    living_norecurrence                   NA's 
                    28                     97                      8 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
     61    1338    2675    2687    3912    5305 

recurrence_status: 
deceased_or_recurrence    living_norecurrence                   NA's 
                    49                     76                      8 

treatment: 
untreated 
      133 

batch: 
 KIU OXFU 
  64   69 

uncurated_author_metadata: 
samplename: KIU_101B88///dataset: UNT///series: KIU///id: 101B88///er: 0///pgr: 0///her2: NA///size: 1.2///node: 0///age: 40///grade: 3///t.dmfs: 2280///e.dmfs: 0///t.rfs: 2280///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_105B13///dataset: UNT///series: KIU///id: 105B13///er: 1///pgr: 1///her2: NA///size: 1.3///node: 0///age: 46///grade: 1///t.dmfs: 2675///e.dmfs: 0///t.rfs: 2675///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
    samplename: KIU_106B55///dataset: UNT///series: KIU///id: 106B55///er: 1///pgr: 1///her2: NA///size: 6///node: 0///age: 37///grade: 1///t.dmfs: 426///e.dmfs: 0///t.rfs: 426///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
  samplename: KIU_111B51///dataset: UNT///series: KIU///id: 111B51///er: 1///pgr: 1///her2: NA///size: 3.3///node: 0///age: 41///grade: 3///t.dmfs: 182///e.dmfs: 1///t.rfs: 182///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_113B11///dataset: UNT///series: KIU///id: 113B11///er: 1///pgr: 1///her2: NA///size: 3.2///node: 0///age: 38///grade: 3///t.dmfs: 1125///e.dmfs: 1///t.rfs: 1125///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_120B73///dataset: UNT///series: KIU///id: 120B73///er: 1///pgr: 1///her2: NA///size: 1.6///node: 0///age: 34///grade: 2///t.dmfs: 3952///e.dmfs: 0///t.rfs: 3952///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_124B25///dataset: UNT///series: KIU///id: 124B25///er: 1///pgr: 1///her2: NA///size: 2.1///node: 0///age: 46///grade: 2///t.dmfs: 1824///e.dmfs: 1///t.rfs: 1824///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
  samplename: KIU_127B00///dataset: UNT///series: KIU///id: 127B00///er: 1///pgr: 1///her2: NA///size: 2.2///node: 0///age: 57///grade: 3///t.dmfs: 699///e.dmfs: 1///t.rfs: 699///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
  samplename: KIU_134B33///dataset: UNT///series: KIU///id: 134B33///er: 1///pgr: 1///her2: NA///size: 2.8///node: 0///age: 63///grade: 2///t.dmfs: 730///e.dmfs: 1///t.rfs: 730///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
  samplename: KIU_136B04///dataset: UNT///series: KIU///id: 136B04///er: 1///pgr: 1///her2: NA///size: 1.7///node: 0///age: 54///grade: 2///t.dmfs: 882///e.dmfs: 0///t.rfs: 882///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_140B91///dataset: UNT///series: KIU///id: 140B91///er: 1///pgr: 1///her2: NA///size: 1.2///node: 0///age: 61///grade: 2///t.dmfs: 2827///e.dmfs: 0///t.rfs: 2827///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_144B49///dataset: UNT///series: KIU///id: 144B49///er: 1///pgr: 1///her2: NA///size: 2.1///node: 0///age: 40///grade: 2///t.dmfs: 2006///e.dmfs: 0///t.rfs: 2006///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_151B84///dataset: UNT///series: KIU///id: 151B84///er: 1///pgr: 1///her2: NA///size: 1.5///node: 0///age: 57///grade: 2///t.dmfs: 2523///e.dmfs: 0///t.rfs: 2523///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_155B52///dataset: UNT///series: KIU///id: 155B52///er: 0///pgr: 1///her2: NA///size: 1.3///node: 0///age: 57///grade: 1///t.dmfs: 3648///e.dmfs: 0///t.rfs: 3648///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_163B27///dataset: UNT///series: KIU///id: 163B27///er: 1///pgr: 1///her2: NA///size: 0.8///node: 0///age: 49///grade: 1///t.dmfs: 2250///e.dmfs: 0///t.rfs: 2250///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_164B81///dataset: UNT///series: KIU///id: 164B81///er: 0///pgr: 0///her2: NA///size: 2.3///node: 0///age: 62///grade: 2///t.dmfs: 3587///e.dmfs: 0///t.rfs: 3587///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_172B19///dataset: UNT///series: KIU///id: 172B19///er: 1///pgr: 1///her2: NA///size: 2.3///node: 0///age: 42///grade: 3///t.dmfs: 3131///e.dmfs: 0///t.rfs: 3131///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_177B67///dataset: UNT///series: KIU///id: 177B67///er: 1///pgr: 1///her2: NA///size: 1.8///node: 0///age: 41///grade: 1///t.dmfs: 2493///e.dmfs: 1///t.rfs: 2493///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
  samplename: KIU_184B38///dataset: UNT///series: KIU///id: 184B38///er: 1///pgr: 1///her2: NA///size: 1///node: 0///age: 63///grade: 1///t.dmfs: 3162///e.dmfs: 0///t.rfs: 3162///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_188B13///dataset: UNT///series: KIU///id: 188B13///er: 0///pgr: 0///her2: NA///size: 1.4///node: 0///age: 60///grade: 2///t.dmfs: 3496///e.dmfs: 0///t.rfs: 3496///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_196B81///dataset: UNT///series: KIU///id: 196B81///er: 1///pgr: 1///her2: NA///size: 1.4///node: 0///age: 65///grade: 1///t.dmfs: 2158///e.dmfs: 0///t.rfs: 2158///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_197B95///dataset: UNT///series: KIU///id: 197B95///er: 1///pgr: 1///her2: NA///size: 1.6///node: 0///age: 44///grade: 2///t.dmfs: 1885///e.dmfs: 0///t.rfs: 1885///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_199B55///dataset: UNT///series: KIU///id: 199B55///er: 1///pgr: 1///her2: NA///size: 2.3///node: 0///age: 54///grade: 2///t.dmfs: 3678///e.dmfs: 0///t.rfs: 3678///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_205B99///dataset: UNT///series: KIU///id: 205B99///er: 1///pgr: 1///her2: NA///size: 2.2///node: 0///age: 59///grade: 1///t.dmfs: 1611///e.dmfs: 1///t.rfs: 1611///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
 samplename: KIU_220C70///dataset: UNT///series: KIU///id: 220C70///er: 1///pgr: NA///her2: NA///size: 2///node: 0///age: 42///grade: 1///t.dmfs: 2797///e.dmfs: 0///t.rfs: 2797///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_227C50///dataset: UNT///series: KIU///id: 227C50///er: 1///pgr: 0///her2: NA///size: 1.2///node: 0///age: 57///grade: 1///t.dmfs: 3314///e.dmfs: 0///t.rfs: 3314///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_229C44///dataset: UNT///series: KIU///id: 229C44///er: 1///pgr: 1///her2: NA///size: 1.3///node: 0///age: 52///grade: 1///t.dmfs: 3466///e.dmfs: 0///t.rfs: 3466///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_231C80///dataset: UNT///series: KIU///id: 231C80///er: 1///pgr: 1///her2: NA///size: 2.2///node: 0///age: 56///grade: 1///t.dmfs: 2341///e.dmfs: 1///t.rfs: 2341///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_233C91///dataset: UNT///series: KIU///id: 233C91///er: 1///pgr: 1///her2: NA///size: 1.1///node: 0///age: 49///grade: 1///t.dmfs: 3344///e.dmfs: 0///t.rfs: 3344///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
  samplename: KIU_242C21///dataset: UNT///series: KIU///id: 242C21///er: 1///pgr: 1///her2: NA///size: 1.6///node: 0///age: 64///grade: 2///t.dmfs: 790///e.dmfs: 0///t.rfs: 790///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_243C70///dataset: UNT///series: KIU///id: 243C70///er: 1///pgr: 1///her2: NA///size: 1.8///node: 0///age: 50///grade: 1///t.dmfs: 2189///e.dmfs: 0///t.rfs: 2189///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
  samplename: KIU_247C76///dataset: UNT///series: KIU///id: 247C76///er: 1///pgr: 1///her2: NA///size: 1///node: 0///age: 56///grade: 2///t.dmfs: 1520///e.dmfs: 0///t.rfs: 1520///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_248C91///dataset: UNT///series: KIU///id: 248C91///er: 1///pgr: 1///her2: NA///size: 2.5///node: 0///age: 57///grade: 1///t.dmfs: 1064///e.dmfs: 0///t.rfs: 1064///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
  samplename: KIU_24C30///dataset: UNT///series: KIU///id: 24C30///er: 1///pgr: 1///her2: NA///size: 2.3///node: 0///age: 55///grade: 2///t.dmfs: 2189///e.dmfs: 0///t.rfs: 2189///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
  samplename: KIU_259C74///dataset: UNT///series: KIU///id: 259C74///er: 1///pgr: 1///her2: NA///size: 1///node: 0///age: 65///grade: 1///t.dmfs: 2462///e.dmfs: 0///t.rfs: 2462///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_260C91///dataset: UNT///series: KIU///id: 260C91///er: 1///pgr: 1///her2: NA///size: 2.1///node: 0///age: 58///grade: 2///t.dmfs: 1642///e.dmfs: 0///t.rfs: 1642///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_266C51///dataset: UNT///series: KIU///id: 266C51///er: 1///pgr: 1///her2: NA///size: 2.2///node: 0///age: 58///grade: 1///t.dmfs: 3222///e.dmfs: 0///t.rfs: 3222///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_268C87///dataset: UNT///series: KIU///id: 268C87///er: 1///pgr: 1///her2: NA///size: 1.5///node: 0///age: 32///grade: 2///t.dmfs: 3253///e.dmfs: 0///t.rfs: 3253///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_272C88///dataset: UNT///series: KIU///id: 272C88///er: 1///pgr: 1///her2: NA///size: 1.7///node: 0///age: 45///grade: 2///t.dmfs: 1338///e.dmfs: 0///t.rfs: 1338///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_278C80///dataset: UNT///series: KIU///id: 278C80///er: 1///pgr: 1///her2: NA///size: 1.1///node: 0///age: 56///grade: 2///t.dmfs: 3162///e.dmfs: 0///t.rfs: 3162///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_279C61///dataset: UNT///series: KIU///id: 279C61///er: 1///pgr: 1///her2: NA///size: 1.9///node: 0///age: 50///grade: 3///t.dmfs: 3192///e.dmfs: 0///t.rfs: 3192///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
  samplename: KIU_280C43///dataset: UNT///series: KIU///id: 280C43///er: 1///pgr: 1///her2: NA///size: 0.9///node: 0///age: 45///grade: 2///t.dmfs: 365///e.dmfs: 0///t.rfs: 365///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
  samplename: KIU_282C51///dataset: UNT///series: KIU///id: 282C51///er: 1///pgr: 1///her2: NA///size: 1.1///node: 0///age: 55///grade: 1///t.dmfs: 851///e.dmfs: 0///t.rfs: 851///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
  samplename: KIU_284C63///dataset: UNT///series: KIU///id: 284C63///er: 1///pgr: 1///her2: NA///size: 1///node: 0///age: 48///grade: 1///t.dmfs: 3435///e.dmfs: 0///t.rfs: 3435///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_286C91///dataset: UNT///series: KIU///id: 286C91///er: 1///pgr: 1///her2: NA///size: 1.8///node: 0///age: 62///grade: 2///t.dmfs: 2675///e.dmfs: 0///t.rfs: 2675///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
    samplename: KIU_28C76///dataset: UNT///series: KIU///id: 28C76///er: 1///pgr: 1///her2: NA///size: 2///node: 0///age: 56///grade: 1///t.dmfs: 2280///e.dmfs: 0///t.rfs: 2280///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
  samplename: KIU_292C66///dataset: UNT///series: KIU///id: 292C66///er: 1///pgr: 1///her2: NA///size: 2///node: 0///age: 51///grade: 2///t.dmfs: 3283///e.dmfs: 0///t.rfs: 3283///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_303C36///dataset: UNT///series: KIU///id: 303C36///er: 0///pgr: 1///her2: NA///size: 2.3///node: 0///age: 37///grade: 3///t.dmfs: 2462///e.dmfs: 0///t.rfs: 2462///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
  samplename: KIU_304C89///dataset: UNT///series: KIU///id: 304C89///er: 0///pgr: 1///her2: NA///size: 1.5///node: 0///age: 54///grade: 3///t.dmfs: 942///e.dmfs: 1///t.rfs: 942///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
  samplename: KIU_308C93///dataset: UNT///series: KIU///id: 308C93///er: 0///pgr: 1///her2: NA///size: 2.1///node: 0///age: 38///grade: 2///t.dmfs: 821///e.dmfs: 1///t.rfs: 821///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_309C49///dataset: UNT///series: KIU///id: 309C49///er: 1///pgr: 1///her2: NA///size: 1.2///node: 0///age: 44///grade: 1///t.dmfs: 3070///e.dmfs: 0///t.rfs: 3070///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
      samplename: KIU_314B55///dataset: UNT///series: KIU///id: 314B55///er: 0///pgr: 0///her2: NA///size: 3///node: 0///age: 38///grade: 3///t.dmfs: 61///e.dmfs: 0///t.rfs: 61///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_316C64///dataset: UNT///series: KIU///id: 316C64///er: 1///pgr: 1///her2: NA///size: 1.3///node: 0///age: 51///grade: 1///t.dmfs: 3374///e.dmfs: 0///t.rfs: 3374///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
  samplename: KIU_36C17///dataset: UNT///series: KIU///id: 36C17///er: 1///pgr: 1///her2: NA///size: 2.2///node: 0///age: 46///grade: 2///t.dmfs: 3010///e.dmfs: 0///t.rfs: 3010///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: KIU_42C67///dataset: UNT///series: KIU///id: 42C67///er: NA///pgr: NA///her2: NA///size: 2.6///node: 0///age: 59///grade: 1///t.dmfs: 3222///e.dmfs: 0///t.rfs: 3222///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
    samplename: KIU_43C47///dataset: UNT///series: KIU///id: 43C47///er: 0///pgr: 1///her2: NA///size: 1.2///node: 0///age: 46///grade: 2///t.dmfs: 213///e.dmfs: 0///t.rfs: 213///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
  samplename: KIU_52A90///dataset: UNT///series: KIU///id: 52A90///er: 1///pgr: 1///her2: NA///size: 2.6///node: 0///age: 53///grade: 1///t.dmfs: 4347///e.dmfs: 0///t.rfs: 4347///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
      samplename: KIU_5B97///dataset: UNT///series: KIU///id: 5B97///er: 1///pgr: 1///her2: NA///size: 2.4///node: 0///age: 37///grade: 2///t.dmfs: 274///e.dmfs: 0///t.rfs: 274///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
  samplename: KIU_65A68///dataset: UNT///series: KIU///id: 65A68///er: 1///pgr: 1///her2: NA///size: 1.8///node: 0///age: 49///grade: 1///t.dmfs: 1277///e.dmfs: 0///t.rfs: 1277///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
  samplename: KIU_74A63///dataset: UNT///series: KIU///id: 74A63///er: 1///pgr: 1///her2: NA///size: 2.2///node: 0///age: 56///grade: 1///t.dmfs: 2158///e.dmfs: 1///t.rfs: 2158///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
  samplename: KIU_86A40///dataset: UNT///series: KIU///id: 86A40///er: 0///pgr: 0///her2: NA///size: 2.4///node: 0///age: 61///grade: 2///t.dmfs: 3648///e.dmfs: 0///t.rfs: 3648///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
  samplename: KIU_87A79///dataset: UNT///series: KIU///id: 87A79///er: 1///pgr: 1///her2: NA///size: 1.2///node: 0///age: 36///grade: 2///t.dmfs: 3709///e.dmfs: 0///t.rfs: 3709///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
  samplename: KIU_88A67///dataset: UNT///series: KIU///id: 88A67///er: 1///pgr: 1///her2: NA///size: 2.4///node: 0///age: 63///grade: 2///t.dmfs: 1550///e.dmfs: 0///t.rfs: 1550///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
  samplename: KIU_89A64///dataset: UNT///series: KIU///id: 89A64///er: 1///pgr: 1///her2: NA///size: 2.3///node: 0///age: 60///grade: 3///t.dmfs: 4165///e.dmfs: 0///t.rfs: 4165///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
    samplename: OXFU_104///dataset: UNT///series: OXFU///id: 104///er: 1///pgr: NA///her2: NA///size: 3.1///node: 0///age: 60///grade: NA///t.dmfs: 648///e.dmfs: 1///t.rfs: 648///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
     samplename: OXFU_1065///dataset: UNT///series: OXFU///id: 1065///er: 0///pgr: NA///her2: NA///size: 2///node: 0///age: 43///grade: 2///t.dmfs: 729///e.dmfs: 1///t.rfs: 729///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
 samplename: OXFU_1183///dataset: UNT///series: OXFU///id: 1183///er: 1///pgr: NA///her2: NA///size: 0.8///node: 0///age: 50///grade: 1///t.dmfs: 1591///e.dmfs: 0///t.rfs: 1591///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
      samplename: OXFU_12///dataset: UNT///series: OXFU///id: 12///er: 1///pgr: NA///her2: NA///size: 0///node: 0///age: 44///grade: NA///t.dmfs: 5305///e.dmfs: 0///t.rfs: 5305///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
 samplename: OXFU_1210///dataset: UNT///series: OXFU///id: 1210///er: 1///pgr: NA///her2: NA///size: 0.8///node: 0///age: 43///grade: 1///t.dmfs: 1920///e.dmfs: 0///t.rfs: 1920///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
   samplename: OXFU_1248///dataset: UNT///series: OXFU///id: 1248///er: 1///pgr: NA///her2: NA///size: 2///node: 0///age: 70///grade: 1///t.dmfs: 3262///e.dmfs: 0///t.rfs: 3262///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
    samplename: OXFU_126///dataset: UNT///series: OXFU///id: 126///er: 0///pgr: NA///her2: NA///size: 1///node: 0///age: 45///grade: NA///t.dmfs: 4163///e.dmfs: 0///t.rfs: 4163///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
   samplename: OXFU_127///dataset: UNT///series: OXFU///id: 127///er: 1///pgr: NA///her2: NA///size: 1.9///node: 0///age: 42///grade: 1///t.dmfs: 4869///e.dmfs: 0///t.rfs: 4869///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
  samplename: OXFU_1286///dataset: UNT///series: OXFU///id: 1286///er: NA///pgr: NA///her2: NA///size: 2///node: 0///age: 52///grade: 1///t.dmfs: 1913///e.dmfs: 0///t.rfs: 1913///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: OXFU_1328///dataset: UNT///series: OXFU///id: 1328///er: 1///pgr: NA///her2: NA///size: 1.3///node: 0///age: 49///grade: NA///t.dmfs: 2648///e.dmfs: 0///t.rfs: 2648///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
     samplename: OXFU_1373///dataset: UNT///series: OXFU///id: 1373///er: 0///pgr: NA///her2: NA///size: 2///node: 0///age: 38///grade: 2///t.dmfs: 268///e.dmfs: 1///t.rfs: 268///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
     samplename: OXFU_138///dataset: UNT///series: OXFU///id: 138///er: 1///pgr: NA///her2: NA///size: 3///node: 0///age: 55///grade: 3///t.dmfs: 5029///e.dmfs: 0///t.rfs: 5029///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
samplename: OXFU_1415///dataset: UNT///series: OXFU///id: 1415///er: 0///pgr: NA///her2: NA///size: 0.9///node: 0///age: 47///grade: NA///t.dmfs: 1290///e.dmfs: 0///t.rfs: 1290///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
   samplename: OXFU_145///dataset: UNT///series: OXFU///id: 145///er: 0///pgr: NA///her2: NA///size: 2.5///node: 0///age: 45///grade: 2///t.dmfs: 4714///e.dmfs: 0///t.rfs: 4714///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
     samplename: OXFU_157///dataset: UNT///series: OXFU///id: 157///er: 0///pgr: NA///her2: NA///size: 2///node: 0///age: 42///grade: 3///t.dmfs: 4906///e.dmfs: 0///t.rfs: 4906///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
     samplename: OXFU_16///dataset: UNT///series: OXFU///id: 16///er: 1///pgr: NA///her2: NA///size: 2.6///node: 0///age: 46///grade: 2///t.dmfs: 4081///e.dmfs: 0///t.rfs: 4081///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
   samplename: OXFU_1605///dataset: UNT///series: OXFU///id: 1605///er: 0///pgr: NA///her2: NA///size: 1///node: 0///age: 39///grade: 2///t.dmfs: 2809///e.dmfs: 0///t.rfs: 2297///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
   samplename: OXFU_181///dataset: UNT///series: OXFU///id: 181///er: 1///pgr: NA///her2: NA///size: 1.5///node: 0///age: 64///grade: 2///t.dmfs: 4579///e.dmfs: 1///t.rfs: 4579///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
       samplename: OXFU_217///dataset: UNT///series: OXFU///id: 217///er: 0///pgr: NA///her2: NA///size: 1///node: 0///age: 53///grade: 2///t.dmfs: 562///e.dmfs: 1///t.rfs: 562///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
     samplename: OXFU_220///dataset: UNT///series: OXFU///id: 220///er: 0///pgr: NA///her2: NA///size: 1///node: 0///age: 47///grade: 2///t.dmfs: 4605///e.dmfs: 0///t.rfs: 4605///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
   samplename: OXFU_223///dataset: UNT///series: OXFU///id: 223///er: 1///pgr: NA///her2: NA///size: 2.1///node: 0///age: 64///grade: 3///t.dmfs: 1875///e.dmfs: 1///t.rfs: 1875///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
  samplename: OXFU_23///dataset: UNT///series: OXFU///id: 23///er: NA///pgr: NA///her2: NA///size: 1.2///node: 0///age: 46///grade: 2///t.dmfs: 4772///e.dmfs: NA///t.rfs: 4772///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
   samplename: OXFU_245///dataset: UNT///series: OXFU///id: 245///er: NA///pgr: NA///her2: NA///size: 1///node: 0///age: 54///grade: NA///t.dmfs: 4859///e.dmfs: 0///t.rfs: 4859///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
     samplename: OXFU_247///dataset: UNT///series: OXFU///id: 247///er: 0///pgr: NA///her2: NA///size: 4.5///node: 0///age: 73///grade: 3///t.dmfs: 562///e.dmfs: 0///t.rfs: 252///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
    samplename: OXFU_254///dataset: UNT///series: OXFU///id: 254///er: 1///pgr: NA///her2: NA///size: 2///node: 0///age: 48///grade: NA///t.dmfs: 4614///e.dmfs: 0///t.rfs: 4614///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
   samplename: OXFU_281///dataset: UNT///series: OXFU///id: 281///er: 0///pgr: NA///her2: NA///size: 1.6///node: 0///age: 64///grade: 2///t.dmfs: 3854///e.dmfs: 0///t.rfs: 3854///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
   samplename: OXFU_316///dataset: UNT///series: OXFU///id: 316///er: 0///pgr: NA///her2: NA///size: 2.2///node: 0///age: 47///grade: 3///t.dmfs: 4683///e.dmfs: 0///t.rfs: 4683///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
     samplename: OXFU_320///dataset: UNT///series: OXFU///id: 320///er: 0///pgr: NA///her2: NA///size: 5///node: 0///age: 39///grade: 3///t.dmfs: 2562///e.dmfs: 0///t.rfs: 2562///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
    samplename: OXFU_337///dataset: UNT///series: OXFU///id: 337///er: 0///pgr: NA///her2: NA///size: 1///node: 0///age: 51///grade: NA///t.dmfs: 769///e.dmfs: NA///t.rfs: 769///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
     samplename: OXFU_348///dataset: UNT///series: OXFU///id: 348///er: 1///pgr: NA///her2: NA///size: 4.5///node: 0///age: 65///grade: 2///t.dmfs: 221///e.dmfs: 1///t.rfs: 221///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
   samplename: OXFU_359///dataset: UNT///series: OXFU///id: 359///er: 0///pgr: NA///her2: NA///size: 0///node: 0///age: 48///grade: 2///t.dmfs: 3939///e.dmfs: NA///t.rfs: 3939///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
     samplename: OXFU_360///dataset: UNT///series: OXFU///id: 360///er: 0///pgr: NA///her2: NA///size: 3///node: 0///age: 32///grade: 3///t.dmfs: 4598///e.dmfs: 0///t.rfs: 4598///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
   samplename: OXFU_366///dataset: UNT///series: OXFU///id: 366///er: 0///pgr: NA///her2: NA///size: 2.5///node: 0///age: 57///grade: 3///t.dmfs: 1112///e.dmfs: 1///t.rfs: 1038///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
     samplename: OXFU_37///dataset: UNT///series: OXFU///id: 37///er: 0///pgr: NA///her2: NA///size: 1.8///node: 0///age: 38///grade: 2///t.dmfs: 3094///e.dmfs: 0///t.rfs: 3094///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
  samplename: OXFU_373///dataset: UNT///series: OXFU///id: 373///er: 1///pgr: NA///her2: NA///size: 1.8///node: 0///age: 64///grade: NA///t.dmfs: 1065///e.dmfs: 1///t.rfs: 1065///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                           1 
                                                                                                                                                                                                                                     (Other) 
                                                                                                                                                                                                                                          34 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse2034; http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse5327 }
\keyword{datasets}
