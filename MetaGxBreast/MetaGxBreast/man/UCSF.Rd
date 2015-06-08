\name{ UCSF_eset }
\alias{ UCSF_eset }
\docType{data}
\title{ UCSF_eset }
\usage{data( UCSF_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=16478745; http://www.ncbi.nlm.nih.gov/pubmed/?term=17401012 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse2990 
  PMIDs: 16478745, 17401012  
  No abstract available.
  notes:
   summary:     
      A gene expression grading index (GGI) was developed. The GGI reclassified grade 2 patients into two groups with low and h
igh risks of recurrence.
   version:     
      2015-04-27 19:46:42
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: probe_1 probe_3 ... probe_10365 (8015 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 8015 features, 162 samples
Platform type: 
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

   29 observations deleted due to missingness 
records   n.max n.start  events  median 0.95LCL 0.95UCL 
 133.00  133.00  133.00   44.00   11.56    9.25      NA 

--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
UCSF_B0241 UCSF_B0244 UCSF_B0251 UCSF_B0296 UCSF_B0299 UCSF_B0305 UCSF_B0317 UCSF_B0318 UCSF_B0326 UCSF_B0327 UCSF_B0334 
         1          1          1          1          1          1          1          1          1          1          1 
UCSF_B0335 UCSF_B0336 UCSF_B0339 UCSF_B0340 UCSF_B0341 UCSF_B0343 UCSF_B0344 UCSF_B0346 UCSF_B0348 UCSF_B0349 UCSF_B0354 
         1          1          1          1          1          1          1          1          1          1          1 
UCSF_B0357 UCSF_B0359 UCSF_B0365 UCSF_B0369 UCSF_B0370 UCSF_B0371 UCSF_B0372 UCSF_B0374 UCSF_B0382 UCSF_B0383 UCSF_B0387 
         1          1          1          1          1          1          1          1          1          1          1 
UCSF_B0394 UCSF_B0410 UCSF_B0413 UCSF_B0414 UCSF_B0419 UCSF_B0420 UCSF_B0427 UCSF_B0428 UCSF_B0433 UCSF_B0447 UCSF_B0454 
         1          1          1          1          1          1          1          1          1          1          1 
UCSF_B0481 UCSF_B0494 UCSF_B0499 UCSF_B0500 UCSF_B0502 UCSF_B0504 UCSF_B0505 UCSF_B0507 UCSF_B0512 UCSF_B0513 UCSF_B0515 
         1          1          1          1          1          1          1          1          1          1          1 
UCSF_B0516 UCSF_B0519 UCSF_B0521 UCSF_B0530 UCSF_B0533 UCSF_B0537 UCSF_B0542 UCSF_B0544 UCSF_B0551 UCSF_B0554 UCSF_B0556 
         1          1          1          1          1          1          1          1          1          1          1 
UCSF_B0564 UCSF_B0565 UCSF_B0566 UCSF_B0567 UCSF_B0574 UCSF_B0624 UCSF_B0626 UCSF_B0627 UCSF_B0629 UCSF_B0630 UCSF_B0638 
         1          1          1          1          1          1          1          1          1          1          1 
UCSF_B0647 UCSF_B0648 UCSF_B0657 UCSF_B0664 UCSF_B0667 UCSF_B0668 UCSF_B0671 UCSF_B0672 UCSF_B0673 UCSF_B0676 UCSF_B0678 
         1          1          1          1          1          1          1          1          1          1          1 
UCSF_B0679 UCSF_B0680 UCSF_B0686 UCSF_B0689 UCSF_B0694 UCSF_B0745 UCSF_B0749 UCSF_B0751 UCSF_B0752 UCSF_B0772 UCSF_B0785 
         1          1          1          1          1          1          1          1          1          1          1 
   (Other) 
        63 

alt_sample_name: 
  B0241   B0244   B0251   B0296   B0299   B0305   B0317   B0318   B0326   B0327   B0334   B0335   B0336   B0339   B0340 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
  B0341   B0343   B0344   B0346   B0348   B0349   B0354   B0357   B0359   B0365   B0369   B0370   B0371   B0372   B0374 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
  B0382   B0383   B0387   B0394   B0410   B0413   B0414   B0419   B0420   B0427   B0428   B0433   B0447   B0454   B0481 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
  B0494   B0499   B0500   B0502   B0504   B0505   B0507   B0512   B0513   B0515   B0516   B0519   B0521   B0530   B0533 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
  B0537   B0542   B0544   B0551   B0554   B0556   B0564   B0565   B0566   B0567   B0574   B0624   B0626   B0627   B0629 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
  B0630   B0638   B0647   B0648   B0657   B0664   B0667   B0668   B0671   B0672   B0673   B0676   B0678   B0679   B0680 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
  B0686   B0689   B0694   B0745   B0749   B0751   B0752   B0772   B0785 (Other) 
      1       1       1       1       1       1       1       1       1      63 

sample_type: 
tumor 
  162 

er: 
negative positive     NA's 
      41      101       20 

pgr: 
negative positive     NA's 
      46       94       22 

her2: 
negative positive     NA's 
      35       19      108 

tumor_size: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  0.000   1.800   2.000   2.682   3.200  11.000       7 

N: 
   0    1 NA's 
  67   82   13 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  28.00   44.00   53.00   56.61   70.00   88.00       9 

grade: 
   1    2    3 NA's 
  14   62   74   12 

dmfs_days: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
     47     897    2040    2084    2992    8267      29 

dmfs_status: 
deceased_or_recurrence    living_norecurrence 
                    22                    140 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
     47     861    1865    1985    2847    8267      29 

recurrence_status: 
deceased_or_recurrence    living_norecurrence 
                    37                    125 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
     47    1087    2054    2140    3087    8267      29 

vital_status: 
deceased   living     NA's 
      54       99        9 

treatment: 
chemo.plus.hormono       chemotherapy     hormonotherapy          untreated               NA's 
                31                 38                 61                 22                 10 

batch: 
UCSF 
 162 

uncurated_author_metadata: 
  samplename: UCSF_B0241///dataset: UCSF///series: UCSF///id: B0241///er: 1///pgr: 1///her2: 0///size: 2///node: 0///age: 68///grade: 3///t.dmfs: 3573///e.dmfs: 0///t.rfs: 3573///e.rfs: 0///t.os: 3573///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0244///dataset: UCSF///series: UCSF///id: B0244///er: 0///pgr: 0///her2: 0///size: 3///node: 0///age: 49///grade: 3///t.dmfs: 4591///e.dmfs: 0///t.rfs: 4591///e.rfs: 0///t.os: 4591///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0251///dataset: UCSF///series: UCSF///id: B0251///er: 1///pgr: 0///her2: 0///size: 3///node: 0///age: 78///grade: 3///t.dmfs: 1474///e.dmfs: 0///t.rfs: 1474///e.rfs: 0///t.os: 1474///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
 samplename: UCSF_B0296///dataset: UCSF///series: UCSF///id: B0296///er: 1///pgr: 1///her2: NA///size: 2///node: 0///age: 82///grade: 2///t.dmfs: 2646///e.dmfs: 0///t.rfs: 2646///e.rfs: 0///t.os: 2646///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
   samplename: UCSF_B0299///dataset: UCSF///series: UCSF///id: B0299///er: 1///pgr: NA///her2: NA///size: 4///node: 1///age: 77///grade: 2///t.dmfs: 861///e.dmfs: 1///t.rfs: 861///e.rfs: 1///t.os: 861///e.os: 1///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0305///dataset: UCSF///series: UCSF///id: B0305///er: 0///pgr: 0///her2: 0///size: 2///node: 1///age: 72///grade: 3///t.dmfs: 1054///e.dmfs: 1///t.rfs: 1054///e.rfs: 1///t.os: 1222///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0317///dataset: UCSF///series: UCSF///id: B0317///er: 1///pgr: 1///her2: 0///size: 2///node: 0///age: 73///grade: 1///t.dmfs: 4219///e.dmfs: 0///t.rfs: 4219///e.rfs: 0///t.os: 4219///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0318///dataset: UCSF///series: UCSF///id: B0318///er: 0///pgr: 1///her2: 1///size: 4///node: 1///age: 40///grade: 2///t.dmfs: 1193///e.dmfs: 0///t.rfs: 1193///e.rfs: 0///t.os: 1193///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0326///dataset: UCSF///series: UCSF///id: B0326///er: 1///pgr: 0///her2: 0///size: 3///node: 1///age: 46///grade: 3///t.dmfs: 2197///e.dmfs: 0///t.rfs: 2197///e.rfs: 0///t.os: 2197///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0327///dataset: UCSF///series: UCSF///id: B0327///er: 0///pgr: 0///her2: 0///size: 2///node: 0///age: 57///grade: 3///t.dmfs: 2292///e.dmfs: 0///t.rfs: 2292///e.rfs: 0///t.os: 2292///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0334///dataset: UCSF///series: UCSF///id: B0334///er: 0///pgr: 1///her2: 1///size: 6///node: 0///age: 54///grade: 2///t.dmfs: 1445///e.dmfs: 1///t.rfs: 1445///e.rfs: 1///t.os: 1445///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0335///dataset: UCSF///series: UCSF///id: B0335///er: 1///pgr: 1///her2: 0///size: 1///node: 0///age: 40///grade: 1///t.dmfs: 2080///e.dmfs: 0///t.rfs: 2080///e.rfs: 0///t.os: 2080///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0336///dataset: UCSF///series: UCSF///id: B0336///er: 1///pgr: 1///her2: 1///size: 4///node: 1///age: 73///grade: 3///t.dmfs: 2456///e.dmfs: 0///t.rfs: 2456///e.rfs: 0///t.os: 2456///e.os: 0///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0339///dataset: UCSF///series: UCSF///id: B0339///er: 1///pgr: 1///her2: 0///size: 3///node: 1///age: 80///grade: 3///t.dmfs: 2467///e.dmfs: 0///t.rfs: 2467///e.rfs: 0///t.os: 2467///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0340///dataset: UCSF///series: UCSF///id: B0340///er: 1///pgr: 1///her2: 0///size: 5///node: 1///age: 39///grade: 3///t.dmfs: 4033///e.dmfs: 0///t.rfs: 4033///e.rfs: 0///t.os: 4033///e.os: 0///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0341///dataset: UCSF///series: UCSF///id: B0341///er: 1///pgr: 1///her2: 0///size: 3///node: 0///age: 37///grade: 2///t.dmfs: 4142///e.dmfs: 0///t.rfs: 1212///e.rfs: 1///t.os: 4142///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                            1 
 samplename: UCSF_B0343///dataset: UCSF///series: UCSF///id: B0343///er: 1///pgr: 0///her2: NA///size: 2///node: 1///age: 62///grade: 3///t.dmfs: 2434///e.dmfs: 0///t.rfs: 2434///e.rfs: 0///t.os: 2434///e.os: 0///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0344///dataset: UCSF///series: UCSF///id: B0344///er: 0///pgr: 0///her2: 0///size: 2///node: 1///age: 51///grade: 3///t.dmfs: 2069///e.dmfs: 1///t.rfs: 2069///e.rfs: 1///t.os: 2069///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0346///dataset: UCSF///series: UCSF///id: B0346///er: 0///pgr: 0///her2: 0///size: 3///node: 0///age: 47///grade: 3///t.dmfs: 1766///e.dmfs: 0///t.rfs: 1766///e.rfs: 0///t.os: 1766///e.os: 0///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0348///dataset: UCSF///series: UCSF///id: B0348///er: 1///pgr: 1///her2: 1///size: 2///node: 0///age: 53///grade: 2///t.dmfs: 2076///e.dmfs: 1///t.rfs: 2076///e.rfs: 1///t.os: 2076///e.os: 1///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                            1 
 samplename: UCSF_B0349///dataset: UCSF///series: UCSF///id: B0349///er: 1///pgr: 1///her2: 1///size: 3///node: NA///age: 83///grade: 2///t.dmfs: 2018///e.dmfs: 0///t.rfs: 2018///e.rfs: 0///t.os: 2018///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
     samplename: UCSF_B0354///dataset: UCSF///series: UCSF///id: B0354///er: 1///pgr: 0///her2: 0///size: 5///node: 1///age: 51///grade: 2///t.dmfs: 748///e.dmfs: 0///t.rfs: 748///e.rfs: 0///t.os: 748///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0357///dataset: UCSF///series: UCSF///id: B0357///er: 1///pgr: 1///her2: 0///size: 2///node: 0///age: 45///grade: 3///t.dmfs: 3912///e.dmfs: 0///t.rfs: 3912///e.rfs: 0///t.os: 3912///e.os: 0///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0359///dataset: UCSF///series: UCSF///id: B0359///er: 0///pgr: 1///her2: 0///size: 2///node: 0///age: 50///grade: 3///t.dmfs: 2054///e.dmfs: 0///t.rfs: 2054///e.rfs: 0///t.os: 2054///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0365///dataset: UCSF///series: UCSF///id: B0365///er: 1///pgr: 1///her2: 0///size: 1///node: 0///age: 48///grade: 3///t.dmfs: 3898///e.dmfs: 0///t.rfs: 3898///e.rfs: 0///t.os: 3898///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
 samplename: UCSF_B0369///dataset: UCSF///series: UCSF///id: B0369///er: 1///pgr: 1///her2: NA///size: 1///node: 0///age: 30///grade: 2///t.dmfs: 1865///e.dmfs: 0///t.rfs: 1865///e.rfs: 0///t.os: 1865///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0370///dataset: UCSF///series: UCSF///id: B0370///er: 1///pgr: 1///her2: 0///size: 4///node: 0///age: 43///grade: 2///t.dmfs: 4007///e.dmfs: 0///t.rfs: 4007///e.rfs: 0///t.os: 4007///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0371///dataset: UCSF///series: UCSF///id: B0371///er: 1///pgr: 1///her2: 0///size: 1///node: 1///age: 63///grade: 2///t.dmfs: 3839///e.dmfs: 0///t.rfs: 3839///e.rfs: 0///t.os: 3839///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
     samplename: UCSF_B0372///dataset: UCSF///series: UCSF///id: B0372///er: 0///pgr: 0///her2: 1///size: 5///node: 1///age: 31///grade: 2///t.dmfs: 244///e.dmfs: 1///t.rfs: 244///e.rfs: 1///t.os: 244///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                            1 
     samplename: UCSF_B0374///dataset: UCSF///series: UCSF///id: B0374///er: 1///pgr: 1///her2: 1///size: 3///node: 1///age: 88///grade: 2///t.dmfs: 281///e.dmfs: 0///t.rfs: 281///e.rfs: 0///t.os: 281///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0382///dataset: UCSF///series: UCSF///id: B0382///er: 0///pgr: 0///her2: 0///size: 2///node: 0///age: 45///grade: 3///t.dmfs: 2233///e.dmfs: 0///t.rfs: 2233///e.rfs: 0///t.os: 2233///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0383///dataset: UCSF///series: UCSF///id: B0383///er: 0///pgr: 1///her2: 0///size: 2///node: 1///age: 36///grade: 3///t.dmfs: 1835///e.dmfs: 0///t.rfs: 1835///e.rfs: 0///t.os: 1835///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0387///dataset: UCSF///series: UCSF///id: B0387///er: 0///pgr: 1///her2: 0///size: 2///node: 0///age: 47///grade: 3///t.dmfs: 1631///e.dmfs: 0///t.rfs: 1631///e.rfs: 0///t.os: 1631///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0394///dataset: UCSF///series: UCSF///id: B0394///er: 0///pgr: 0///her2: 0///size: 2///node: 1///age: 39///grade: 3///t.dmfs: 3712///e.dmfs: 0///t.rfs: 3712///e.rfs: 0///t.os: 3712///e.os: 0///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                            1 
       samplename: UCSF_B0410///dataset: UCSF///series: UCSF///id: B0410///er: 1///pgr: 1///her2: 0///size: 0///node: NA///age: 42///grade: 2///t.dmfs: NA///e.dmfs: 0///t.rfs: NA///e.rfs: 0///t.os: NA///e.os: 1///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                            1 
        samplename: UCSF_B0413///dataset: UCSF///series: UCSF///id: B0413///er: 1///pgr: 1///her2: 1///size: 8///node: 1///age: 47///grade: 2///t.dmfs: NA///e.dmfs: 1///t.rfs: NA///e.rfs: 1///t.os: NA///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0414///dataset: UCSF///series: UCSF///id: B0414///er: 1///pgr: 1///her2: 1///size: 1///node: 0///age: 43///grade: 2///t.dmfs: 3642///e.dmfs: 0///t.rfs: 3642///e.rfs: 0///t.os: 3642///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0419///dataset: UCSF///series: UCSF///id: B0419///er: 1///pgr: 1///her2: 0///size: 5///node: 1///age: 58///grade: 2///t.dmfs: 1945///e.dmfs: 0///t.rfs: 1945///e.rfs: 0///t.os: 1945///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0420///dataset: UCSF///series: UCSF///id: B0420///er: 1///pgr: 1///her2: 0///size: 2///node: 0///age: 37///grade: 1///t.dmfs: 1566///e.dmfs: 1///t.rfs: 1566///e.rfs: 1///t.os: 1580///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
    samplename: UCSF_B0427///dataset: UCSF///series: UCSF///id: B0427///er: 1///pgr: 1///her2: 1///size: 2///node: 1///age: 41///grade: 3///t.dmfs: 680///e.dmfs: 1///t.rfs: 680///e.rfs: 1///t.os: 1087///e.os: 1///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0428///dataset: UCSF///series: UCSF///id: B0428///er: 1///pgr: 1///her2: 0///size: 4///node: 1///age: 38///grade: 3///t.dmfs: 3606///e.dmfs: 0///t.rfs: 3606///e.rfs: 0///t.os: 3606///e.os: 0///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                            1 
     samplename: UCSF_B0433///dataset: UCSF///series: UCSF///id: B0433///er: 0///pgr: 0///her2: 0///size: 4///node: 0///age: 60///grade: 3///t.dmfs: 483///e.dmfs: 1///t.rfs: 483///e.rfs: 1///t.os: 817///e.os: 1///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                            1 
    samplename: UCSF_B0447///dataset: UCSF///series: UCSF///id: B0447///er: 0///pgr: 0///her2: NA///size: 2///node: 1///age: 62///grade: 3///t.dmfs: 897///e.dmfs: 0///t.rfs: 897///e.rfs: 0///t.os: 897///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                            1 
    samplename: UCSF_B0454///dataset: UCSF///series: UCSF///id: B0454///er: NA///pgr: NA///her2: NA///size: 1///node: 0///age: 80///grade: NA///t.dmfs: NA///e.dmfs: 0///t.rfs: NA///e.rfs: 1///t.os: NA///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                            1 
      samplename: UCSF_B0481///dataset: UCSF///series: UCSF///id: B0481///er: NA///pgr: NA///her2: 1///size: 0///node: 0///age: 33///grade: 3///t.dmfs: NA///e.dmfs: 0///t.rfs: NA///e.rfs: 0///t.os: NA///e.os: 0///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                            1 
 samplename: UCSF_B0494///dataset: UCSF///series: UCSF///id: B0494///er: 0///pgr: 0///her2: NA///size: 4///node: 0///age: 44///grade: 3///t.dmfs: 3113///e.dmfs: 0///t.rfs: 3113///e.rfs: 0///t.os: 3113///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0499///dataset: UCSF///series: UCSF///id: B0499///er: 1///pgr: 1///her2: 1///size: 2///node: 1///age: 71///grade: 2///t.dmfs: 2182///e.dmfs: 0///t.rfs: 2182///e.rfs: 0///t.os: 2182///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
     samplename: UCSF_B0500///dataset: UCSF///series: UCSF///id: B0500///er: 0///pgr: 1///her2: 0///size: 0///node: 0///age: 28///grade: 2///t.dmfs: 624///e.dmfs: 1///t.rfs: 624///e.rfs: 1///t.os: 766///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0502///dataset: UCSF///series: UCSF///id: B0502///er: 0///pgr: 0///her2: 1///size: 2///node: 0///age: 55///grade: 3///t.dmfs: 2883///e.dmfs: 0///t.rfs: 2883///e.rfs: 0///t.os: 2883///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
 samplename: UCSF_B0504///dataset: UCSF///series: UCSF///id: B0504///er: 0///pgr: 0///her2: NA///size: 3///node: 1///age: 49///grade: 3///t.dmfs: 1116///e.dmfs: 0///t.rfs: 1116///e.rfs: 0///t.os: 1116///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                            1 
        samplename: UCSF_B0505///dataset: UCSF///series: UCSF///id: B0505///er: 1///pgr: 1///her2: 1///size: 2///node: 1///age: 69///grade: 2///t.dmfs: NA///e.dmfs: 0///t.rfs: NA///e.rfs: 0///t.os: NA///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
   samplename: UCSF_B0507///dataset: UCSF///series: UCSF///id: B0507///er: 0///pgr: 1///her2: NA///size: 11///node: 0///age: 28///grade: 3///t.dmfs: 624///e.dmfs: 1///t.rfs: 624///e.rfs: 1///t.os: 766///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
 samplename: UCSF_B0512///dataset: UCSF///series: UCSF///id: B0512///er: 1///pgr: 1///her2: NA///size: 3///node: 1///age: 47///grade: 3///t.dmfs: 1372///e.dmfs: 0///t.rfs: 1372///e.rfs: 1///t.os: 1372///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0513///dataset: UCSF///series: UCSF///id: B0513///er: 1///pgr: 1///her2: 0///size: 2///node: 0///age: 37///grade: 2///t.dmfs: 4142///e.dmfs: 0///t.rfs: 1212///e.rfs: 1///t.os: 4142///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0515///dataset: UCSF///series: UCSF///id: B0515///er: 0///pgr: 0///her2: 1///size: 2///node: 0///age: 56///grade: 3///t.dmfs: 1058///e.dmfs: 0///t.rfs: 1058///e.rfs: 0///t.os: 1058///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
 samplename: UCSF_B0516///dataset: UCSF///series: UCSF///id: B0516///er: 1///pgr: 1///her2: 0///size: 1///node: NA///age: 86///grade: 1///t.dmfs: 1609///e.dmfs: 0///t.rfs: 1609///e.rfs: 0///t.os: 1609///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0519///dataset: UCSF///series: UCSF///id: B0519///er: 1///pgr: 1///her2: 1///size: 2///node: 0///age: 56///grade: 3///t.dmfs: 1598///e.dmfs: 0///t.rfs: 1598///e.rfs: 0///t.os: 1598///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0521///dataset: UCSF///series: UCSF///id: B0521///er: 1///pgr: 1///her2: 1///size: 2///node: 1///age: 56///grade: 2///t.dmfs: 1602///e.dmfs: 0///t.rfs: 1602///e.rfs: 0///t.os: 1602///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
    samplename: UCSF_B0530///dataset: UCSF///series: UCSF///id: B0530///er: 1///pgr: 1///her2: 1///size: 2///node: 1///age: 50///grade: 3///t.dmfs: 857///e.dmfs: 1///t.rfs: 857///e.rfs: 1///t.os: 1503///e.os: 1///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                            1 
 samplename: UCSF_B0533///dataset: UCSF///series: UCSF///id: B0533///er: 1///pgr: 1///her2: NA///size: 1///node: 1///age: 75///grade: 2///t.dmfs: 1069///e.dmfs: 0///t.rfs: 1069///e.rfs: 0///t.os: 1069///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
    samplename: UCSF_B0537///dataset: UCSF///series: UCSF///id: B0537///er: 0///pgr: 1///her2: NA///size: 2///node: 1///age: 36///grade: 3///t.dmfs: 653///e.dmfs: 0///t.rfs: 653///e.rfs: 0///t.os: 653///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                            1 
     samplename: UCSF_B0542///dataset: UCSF///series: UCSF///id: B0542///er: 1///pgr: 0///her2: 1///size: 8///node: 1///age: 33///grade: 3///t.dmfs: 288///e.dmfs: 1///t.rfs: 288///e.rfs: 1///t.os: 288///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
     samplename: UCSF_B0544///dataset: UCSF///series: UCSF///id: B0544///er: NA///pgr: NA///her2: NA///size: 0///node: 0///age: 49///grade: 1///t.dmfs: NA///e.dmfs: 0///t.rfs: NA///e.rfs: 0///t.os: NA///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                            1 
       samplename: UCSF_B0551///dataset: UCSF///series: UCSF///id: B0551///er: 0///pgr: 0///her2: NA///size: 5///node: 1///age: 37///grade: 3///t.dmfs: NA///e.dmfs: 1///t.rfs: NA///e.rfs: 1///t.os: NA///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                            1 
      samplename: UCSF_B0554///dataset: UCSF///series: UCSF///id: B0554///er: NA///pgr: NA///her2: 0///size: 4///node: 1///age: 43///grade: 3///t.dmfs: NA///e.dmfs: 0///t.rfs: NA///e.rfs: 1///t.os: NA///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                            1 
    samplename: UCSF_B0556///dataset: UCSF///series: UCSF///id: B0556///er: 1///pgr: 0///her2: NA///size: 2///node: 1///age: 37///grade: 3///t.dmfs: 496///e.dmfs: 0///t.rfs: 496///e.rfs: 0///t.os: 496///e.os: 0///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                            1 
       samplename: UCSF_B0564///dataset: UCSF///series: UCSF///id: B0564///er: 1///pgr: 1///her2: NA///size: 2///node: 1///age: 78///grade: 3///t.dmfs: NA///e.dmfs: 0///t.rfs: NA///e.rfs: 0///t.os: NA///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
       samplename: UCSF_B0565///dataset: UCSF///series: UCSF///id: B0565///er: 1///pgr: 1///her2: 0///size: 3///node: NA///age: 86///grade: 1///t.dmfs: NA///e.dmfs: 0///t.rfs: NA///e.rfs: 0///t.os: NA///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
        samplename: UCSF_B0566///dataset: UCSF///series: UCSF///id: B0566///er: 1///pgr: 1///her2: 0///size: 5///node: 1///age: 55///grade: 2///t.dmfs: NA///e.dmfs: 0///t.rfs: NA///e.rfs: 1///t.os: NA///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                            1 
        samplename: UCSF_B0567///dataset: UCSF///series: UCSF///id: B0567///er: 1///pgr: 0///her2: 0///size: 2///node: 0///age: 70///grade: 2///t.dmfs: NA///e.dmfs: 0///t.rfs: NA///e.rfs: 0///t.os: NA///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
       samplename: UCSF_B0574///dataset: UCSF///series: UCSF///id: B0574///er: 0///pgr: 0///her2: NA///size: 1///node: 0///age: 69///grade: 3///t.dmfs: NA///e.dmfs: 0///t.rfs: NA///e.rfs: 0///t.os: NA///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
       samplename: UCSF_B0624///dataset: UCSF///series: UCSF///id: B0624///er: 1///pgr: 1///her2: NA///size: 8///node: 1///age: 45///grade: 3///t.dmfs: NA///e.dmfs: 0///t.rfs: NA///e.rfs: 1///t.os: NA///e.os: 1///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                            1 
   samplename: UCSF_B0626///dataset: UCSF///series: UCSF///id: B0626///er: 1///pgr: NA///her2: NA///size: 5///node: 0///age: 86///grade: 1///t.dmfs: 766///e.dmfs: 0///t.rfs: 766///e.rfs: 0///t.os: 766///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                            1 
    samplename: UCSF_B0627///dataset: UCSF///series: UCSF///id: B0627///er: NA///pgr: NA///her2: NA///size: 3///node: NA///age: 66///grade: 3///t.dmfs: NA///e.dmfs: 0///t.rfs: NA///e.rfs: 0///t.os: NA///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
     samplename: UCSF_B0629///dataset: UCSF///series: UCSF///id: B0629///er: NA///pgr: NA///her2: NA///size: 2///node: 1///age: 88///grade: 3///t.dmfs: NA///e.dmfs: 0///t.rfs: NA///e.rfs: 0///t.os: NA///e.os: 1///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
 samplename: UCSF_B0630///dataset: UCSF///series: UCSF///id: B0630///er: 1///pgr: 1///her2: NA///size: 1///node: 0///age: 44///grade: 3///t.dmfs: 1923///e.dmfs: 0///t.rfs: 1923///e.rfs: 0///t.os: 1923///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                            1 
 samplename: UCSF_B0638///dataset: UCSF///series: UCSF///id: B0638///er: 1///pgr: 1///her2: NA///size: 3///node: 1///age: 66///grade: 2///t.dmfs: 1766///e.dmfs: 0///t.rfs: 1766///e.rfs: 0///t.os: 1766///e.os: 1///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                            1 
   samplename: UCSF_B0647///dataset: UCSF///series: UCSF///id: B0647///er: NA///pgr: NA///her2: NA///size: 1///node: 0///age: NA///grade: 1///t.dmfs: NA///e.dmfs: 0///t.rfs: NA///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0648///dataset: UCSF///series: UCSF///id: B0648///er: NA///pgr: NA///her2: NA///size: 0///node: 1///age: 69///grade: 2///t.dmfs: 620///e.dmfs: 0///t.rfs: 620///e.rfs: 0///t.os: 620///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
 samplename: UCSF_B0657///dataset: UCSF///series: UCSF///id: B0657///er: 1///pgr: 1///her2: NA///size: 3///node: 1///age: 78///grade: 1///t.dmfs: 1573///e.dmfs: 0///t.rfs: 1573///e.rfs: 0///t.os: 1573///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
       samplename: UCSF_B0664///dataset: UCSF///series: UCSF///id: B0664///er: 1///pgr: 1///her2: NA///size: 1///node: 0///age: 47///grade: 3///t.dmfs: 47///e.dmfs: 0///t.rfs: 47///e.rfs: 0///t.os: 47///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0667///dataset: UCSF///series: UCSF///id: B0667///er: NA///pgr: NA///her2: NA///size: 3///node: 0///age: 44///grade: 3///t.dmfs: 412///e.dmfs: 0///t.rfs: 412///e.rfs: 0///t.os: 412///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                            1 
    samplename: UCSF_B0668///dataset: UCSF///series: UCSF///id: B0668///er: 1///pgr: 1///her2: NA///size: 2///node: 0///age: 42///grade: 1///t.dmfs: 156///e.dmfs: 0///t.rfs: 156///e.rfs: 0///t.os: 156///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                            1 
    samplename: UCSF_B0671///dataset: UCSF///series: UCSF///id: B0671///er: 1///pgr: 1///her2: NA///size: 2///node: 0///age: 57///grade: 2///t.dmfs: 149///e.dmfs: 0///t.rfs: 149///e.rfs: 0///t.os: 149///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                            1 
    samplename: UCSF_B0672///dataset: UCSF///series: UCSF///id: B0672///er: 0///pgr: 0///her2: NA///size: 3///node: 1///age: 60///grade: 3///t.dmfs: 671///e.dmfs: 0///t.rfs: 671///e.rfs: 0///t.os: 671///e.os: 1///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                            1 
   samplename: UCSF_B0673///dataset: UCSF///series: UCSF///id: B0673///er: NA///pgr: NA///her2: NA///size: 1///node: 0///age: NA///grade: 2///t.dmfs: NA///e.dmfs: 0///t.rfs: NA///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                            1 
 samplename: UCSF_B0676///dataset: UCSF///series: UCSF///id: B0676///er: NA///pgr: NA///her2: NA///size: 0///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: 0///t.rfs: NA///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                            1 
    samplename: UCSF_B0678///dataset: UCSF///series: UCSF///id: B0678///er: 1///pgr: 1///her2: NA///size: 2///node: 1///age: 40///grade: 1///t.dmfs: 127///e.dmfs: 0///t.rfs: 127///e.rfs: 0///t.os: 127///e.os: 0///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                            1 
 samplename: UCSF_B0679///dataset: UCSF///series: UCSF///id: B0679///er: 1///pgr: 1///her2: NA///size: 2///node: 1///age: 53///grade: 2///t.dmfs: 1722///e.dmfs: 0///t.rfs: 1722///e.rfs: 0///t.os: 1722///e.os: 0///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                            1 
       samplename: UCSF_B0680///dataset: UCSF///series: UCSF///id: B0680///er: 1///pgr: 1///her2: NA///size: 1///node: 1///age: 48///grade: 1///t.dmfs: NA///e.dmfs: 0///t.rfs: NA///e.rfs: 0///t.os: NA///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0686///dataset: UCSF///series: UCSF///id: B0686///er: 1///pgr: 1///her2: NA///size: 2///node: 0///age: 78///grade: 2///t.dmfs: 1759///e.dmfs: 0///t.rfs: 350///e.rfs: 1///t.os: 1759///e.os: 0///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                            1 
       samplename: UCSF_B0689///dataset: UCSF///series: UCSF///id: B0689///er: 0///pgr: 0///her2: NA///size: 1///node: 0///age: 54///grade: 3///t.dmfs: 91///e.dmfs: 0///t.rfs: 91///e.rfs: 0///t.os: 91///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                            1 
   samplename: UCSF_B0694///dataset: UCSF///series: UCSF///id: B0694///er: NA///pgr: NA///her2: NA///size: 5///node: 0///age: NA///grade: 3///t.dmfs: NA///e.dmfs: 0///t.rfs: NA///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                            1 
     samplename: UCSF_B0745///dataset: UCSF///series: UCSF///id: B0745///er: 0///pgr: 0///her2: NA///size: 1.4///node: 0///age: 38///grade: 3///t.dmfs: 62///e.dmfs: 0///t.rfs: 62///e.rfs: 0///t.os: 62///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                            1 
     samplename: UCSF_B0749///dataset: UCSF///series: UCSF///id: B0749///er: 1///pgr: 1///her2: NA///size: 3.5///node: 1///age: 70///grade: 3///t.dmfs: NA///e.dmfs: 0///t.rfs: NA///e.rfs: 1///t.os: NA///e.os: 0///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                            1 
samplename: UCSF_B0751///dataset: UCSF///series: UCSF///id: B0751///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: 0///t.rfs: NA///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                            1 
      samplename: UCSF_B0752///dataset: UCSF///series: UCSF///id: B0752///er: 1///pgr: 1///her2: NA///size: 4///node: 0///age: 72///grade: 3///t.dmfs: NA///e.dmfs: 0///t.rfs: NA///e.rfs: 1///t.os: NA///e.os: 0///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                            1 
  samplename: UCSF_B0772///dataset: UCSF///series: UCSF///id: B0772///er: 1///pgr: 0///her2: NA///size: 6.2///node: 1///age: 52///grade: 2///t.dmfs: 244///e.dmfs: 0///t.rfs: 244///e.rfs: 0///t.os: 244///e.os: 0///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                            1 
samplename: UCSF_B0785///dataset: UCSF///series: UCSF///id: B0785///er: NA///pgr: NA///her2: NA///size: 3.2///node: NA///age: NA///grade: 2///t.dmfs: NA///e.dmfs: 0///t.rfs: NA///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                            1 
                                                                                                                                                                                                                                      (Other) 
                                                                                                                                                                                                                                           63 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse2990 }
\keyword{datasets}
