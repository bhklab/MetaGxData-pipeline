\name{ MDA4_eset }
\alias{ MDA4_eset }
\docType{data}
\title{ MDA4_eset }
\usage{data( MDA4_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=16896004 
  Title:  
  URL: http://bioinformatics.mdanderson.org/pubdata.html 
  PMIDs: 16896004 
  No abstract available.
  notes:
   summary:     
      The developed 30-probe set has high sensitivity and negative predictive value, accurately identifying 12 out of 13 patien
ts with pCR and 27 out of 28 patients with residual disease.
   version:     
      2015-04-27 19:19:46
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: 1007_s_at 1053_at ... AFFX-HUMISGF3A/M97935_MB_at (21169 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 21169 features, 129 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
    MDA4_M141     MDA4_M150     MDA4_M172     MDA4_M194     MDA4_M207     MDA4_M208     MDA4_M248     MDA4_M294     MDA4_M313 
            1             1             1             1             1             1             1             1             1 
    MDA4_M315     MDA4_M322     MDA4_M323 MDA4_M323_bis     MDA4_M330     MDA4_M331     MDA4_M332     MDA4_M333     MDA4_M334 
            1             1             1             1             1             1             1             1             1 
    MDA4_M337     MDA4_M338     MDA4_M339     MDA4_M349     MDA4_M353     MDA4_M357     MDA4_M367     MDA4_M369     MDA4_M371 
            1             1             1             1             1             1             1             1             1 
    MDA4_M373     MDA4_M375     MDA4_M378     MDA4_M380     MDA4_M382     MDA4_M384     MDA4_M385     MDA4_M386     MDA4_M387 
            1             1             1             1             1             1             1             1             1 
    MDA4_M395     MDA4_M396     MDA4_M399     MDA4_M400     MDA4_M402     MDA4_M406     MDA4_M414     MDA4_M421     MDA4_M430 
            1             1             1             1             1             1             1             1             1 
    MDA4_M431     MDA4_M433     MDA4_M439     MDA4_M444     MDA4_M447     MDA4_M452     MDA4_M458     MDA4_M461     MDA4_M463 
            1             1             1             1             1             1             1             1             1 
    MDA4_M464     MDA4_M466     MDA4_M473     MDA4_M484     MDA4_M485     MDA4_M497     MDA4_M503     MDA4_M506     MDA4_M508 
            1             1             1             1             1             1             1             1             1 
    MDA4_M521     MDA4_M523     MDA4_M524     MDA4_M536     MDA4_M557     MDA4_M559     MDA4_M564     MDA4_M570     MDA4_M572 
            1             1             1             1             1             1             1             1             1 
    MDA4_M578     MDA4_M595     MDA4_M604     MDA4_M610     MDA4_M631     MDA4_M632     MDA4_M633     MDA4_M644     MDA4_M651 
            1             1             1             1             1             1             1             1             1 
    MDA4_M652     MDA4_M658     MDA4_M662     MDA4_M665     MDA4_M675     MDA4_M679     MDA4_M680     MDA4_M683     MDA4_M690 
            1             1             1             1             1             1             1             1             1 
    MDA4_M694     MDA4_M696     MDA4_M702     MDA4_M703     MDA4_M709     MDA4_M717     MDA4_M718     MDA4_M728     MDA4_M729 
            1             1             1             1             1             1             1             1             1 
      (Other) 
           30 

unique_patient_ID: 
    323     141     150     172     194     207     208     248     294     313     315     322     330     331     332 
      2       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
    333     334     337     338     339     349     353     357     367     369     371     373     375     378     380 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
    382     384     385     386     387     395     396     399     400     402     406     414     421     430     431 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
    433     439     444     447     452     458     461     463     464     466     473     484     485     497     503 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
    506     508     521     523     524     536     557     559     564     570     572     578     595     604     610 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
    631     632     633     644     651     652     658     662     665     675     679     680     683     690     694 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
    696     702     703     709     717     718     728     729     732 (Other) 
      1       1       1       1       1       1       1       1       1      29 

sample_type: 
tumor 
  129 

er: 
negative positive     NA's 
      48       79        2 

pgr: 
negative positive     NA's 
      73       54        2 

her2: 
negative positive 
     114       15 

tumor_size: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  0.000   0.500   1.800   2.162   3.000  10.000       8 

N: 
   0    1 NA's 
  59   62    8 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  28.00   43.00   51.00   51.43   61.00   73.00 

treatment: 
chemotherapy 
         129 

batch: 
MDA4 
 129 

uncurated_author_metadata: 
      samplename: MDA4_M141///dataset: MDA4///series: MDA4///id: 141///er: 1///pgr: 0///her2: 1///size: 0///node: 0///age: 70///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M150///dataset: MDA4///series: MDA4///id: 150///er: 0///pgr: 0///her2: 1///size: 1.7///node: 1///age: 30///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M172///dataset: MDA4///series: MDA4///id: 172///er: 0///pgr: 0///her2: 0///size: 1.5///node: 0///age: 33///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M194///dataset: MDA4///series: MDA4///id: 194///er: 1///pgr: 1///her2: 1///size: 0///node: 0///age: 62///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M207///dataset: MDA4///series: MDA4///id: 207///er: 0///pgr: 1///her2: 0///size: 2///node: 0///age: 51///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M208///dataset: MDA4///series: MDA4///id: 208///er: 1///pgr: 1///her2: 0///size: 2.5///node: 0///age: 49///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M248///dataset: MDA4///series: MDA4///id: 248///er: 1///pgr: 1///her2: 1///size: 1.8///node: 0///age: 68///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M294///dataset: MDA4///series: MDA4///id: 294///er: 0///pgr: 0///her2: 1///size: 2///node: 0///age: 53///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M313///dataset: MDA4///series: MDA4///id: 313///er: 1///pgr: 1///her2: 0///size: 1.3///node: 0///age: 73///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M315///dataset: MDA4///series: MDA4///id: 315///er: 1///pgr: 1///her2: 0///size: 3///node: 1///age: 36///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M322///dataset: MDA4///series: MDA4///id: 322///er: 1///pgr: 0///her2: 1///size: 4.5///node: 0///age: 65///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
samplename: MDA4_M323_bis///dataset: MDA4///series: MDA4///id: 323///er: 1///pgr: 1///her2: 0///size: 1.6///node: 1///age: 40///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M323///dataset: MDA4///series: MDA4///id: 323///er: 1///pgr: 1///her2: 0///size: 1.6///node: 1///age: 40///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M330///dataset: MDA4///series: MDA4///id: 330///er: 1///pgr: 1///her2: 0///size: 0///node: 1///age: 35///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M331///dataset: MDA4///series: MDA4///id: 331///er: 0///pgr: 0///her2: 1///size: 0///node: 0///age: 71///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M332///dataset: MDA4///series: MDA4///id: 332///er: 1///pgr: 1///her2: 0///size: 0.2///node: 0///age: 42///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M333///dataset: MDA4///series: MDA4///id: 333///er: 1///pgr: 1///her2: 0///size: 1.3///node: 1///age: 67///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M334///dataset: MDA4///series: MDA4///id: 334///er: 0///pgr: 0///her2: 0///size: 1.3///node: 1///age: 63///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M337///dataset: MDA4///series: MDA4///id: 337///er: 1///pgr: 1///her2: 0///size: 4///node: 1///age: 52///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M338///dataset: MDA4///series: MDA4///id: 338///er: 1///pgr: 1///her2: 0///size: 3///node: 1///age: 45///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M339///dataset: MDA4///series: MDA4///id: 339///er: 0///pgr: 0///her2: 0///size: 0.5///node: 0///age: 44///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M349///dataset: MDA4///series: MDA4///id: 349///er: 1///pgr: 1///her2: 0///size: 1.1///node: 1///age: 40///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M353///dataset: MDA4///series: MDA4///id: 353///er: 1///pgr: 0///her2: 0///size: 0///node: 0///age: 51///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M357///dataset: MDA4///series: MDA4///id: 357///er: 1///pgr: 1///her2: 0///size: 4///node: 1///age: 62///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M367///dataset: MDA4///series: MDA4///id: 367///er: 1///pgr: 0///her2: 0///size: 0///node: 0///age: 61///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M369///dataset: MDA4///series: MDA4///id: 369///er: 0///pgr: 0///her2: 0///size: 0.5///node: 0///age: 52///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M371///dataset: MDA4///series: MDA4///id: 371///er: 0///pgr: 0///her2: 1///size: 0///node: 0///age: 51///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M373///dataset: MDA4///series: MDA4///id: 373///er: 0///pgr: 0///her2: 0///size: 0///node: 0///age: 45///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M375///dataset: MDA4///series: MDA4///id: 375///er: 0///pgr: 0///her2: 0///size: 0.7///node: 1///age: 63///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M378///dataset: MDA4///series: MDA4///id: 378///er: 1///pgr: 1///her2: 0///size: 3.5///node: 1///age: 50///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M380///dataset: MDA4///series: MDA4///id: 380///er: 1///pgr: 1///her2: 0///size: 3.8///node: 1///age: 50///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M382///dataset: MDA4///series: MDA4///id: 382///er: 0///pgr: 0///her2: 1///size: 0///node: 0///age: 55///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M384///dataset: MDA4///series: MDA4///id: 384///er: 1///pgr: 1///her2: 0///size: 8///node: 0///age: 60///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M385///dataset: MDA4///series: MDA4///id: 385///er: 1///pgr: 1///her2: 1///size: 2///node: 1///age: 65///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M386///dataset: MDA4///series: MDA4///id: 386///er: 1///pgr: 0///her2: 0///size: 0.9///node: 0///age: 42///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M387///dataset: MDA4///series: MDA4///id: 387///er: 1///pgr: 1///her2: 1///size: 4///node: 1///age: 28///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M395///dataset: MDA4///series: MDA4///id: 395///er: 1///pgr: 1///her2: 0///size: NA///node: NA///age: 58///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M396///dataset: MDA4///series: MDA4///id: 396///er: 1///pgr: 0///her2: 0///size: 3.5///node: 1///age: 35///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M399///dataset: MDA4///series: MDA4///id: 399///er: 1///pgr: 1///her2: 0///size: 3///node: 1///age: 40///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M400///dataset: MDA4///series: MDA4///id: 400///er: 1///pgr: 1///her2: 0///size: 2///node: 1///age: 55///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M402///dataset: MDA4///series: MDA4///id: 402///er: 1///pgr: 0///her2: 0///size: 0///node: 0///age: 38///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M406///dataset: MDA4///series: MDA4///id: 406///er: 0///pgr: 0///her2: 0///size: 1.7///node: 1///age: 58///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M414///dataset: MDA4///series: MDA4///id: 414///er: 1///pgr: 0///her2: 0///size: 2///node: 0///age: 66///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M421///dataset: MDA4///series: MDA4///id: 421///er: 1///pgr: 1///her2: 0///size: 3///node: 1///age: 48///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M430///dataset: MDA4///series: MDA4///id: 430///er: 0///pgr: 0///her2: 0///size: 1.1///node: 0///age: 52///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M431///dataset: MDA4///series: MDA4///id: 431///er: 1///pgr: 1///her2: 0///size: 2.5///node: 1///age: 50///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M433///dataset: MDA4///series: MDA4///id: 433///er: 1///pgr: 1///her2: 0///size: 2.8///node: 1///age: 47///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M439///dataset: MDA4///series: MDA4///id: 439///er: 0///pgr: 0///her2: 1///size: 3.3///node: 0///age: 73///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M444///dataset: MDA4///series: MDA4///id: 444///er: 1///pgr: 0///her2: 0///size: 0.2///node: 1///age: 53///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M447///dataset: MDA4///series: MDA4///id: 447///er: 1///pgr: 0///her2: 0///size: 2.5///node: 0///age: 44///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M452///dataset: MDA4///series: MDA4///id: 452///er: 1///pgr: 0///her2: 0///size: NA///node: NA///age: 42///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M458///dataset: MDA4///series: MDA4///id: 458///er: 0///pgr: 0///her2: 0///size: 0.9///node: 0///age: 65///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M461///dataset: MDA4///series: MDA4///id: 461///er: 1///pgr: 1///her2: 0///size: 0.5///node: 1///age: 44///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M463///dataset: MDA4///series: MDA4///id: 463///er: 1///pgr: 0///her2: 0///size: 2.2///node: 1///age: 60///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M464///dataset: MDA4///series: MDA4///id: 464///er: 0///pgr: 0///her2: 0///size: NA///node: NA///age: 73///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M466///dataset: MDA4///series: MDA4///id: 466///er: 1///pgr: 1///her2: 0///size: 4///node: 1///age: 64///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M473///dataset: MDA4///series: MDA4///id: 473///er: 0///pgr: 0///her2: 0///size: 1.5///node: 0///age: 54///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M484///dataset: MDA4///series: MDA4///id: 484///er: 1///pgr: 0///her2: 0///size: 4///node: 1///age: 55///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M485///dataset: MDA4///series: MDA4///id: 485///er: 1///pgr: 0///her2: 0///size: 0///node: 0///age: 50///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M497///dataset: MDA4///series: MDA4///id: 497///er: 1///pgr: 0///her2: 0///size: 4///node: 0///age: 62///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M503///dataset: MDA4///series: MDA4///id: 503///er: 0///pgr: 0///her2: 0///size: 1.1///node: 0///age: 38///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M506///dataset: MDA4///series: MDA4///id: 506///er: 1///pgr: 1///her2: 0///size: 1.1///node: 1///age: 46///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M508///dataset: MDA4///series: MDA4///id: 508///er: 1///pgr: 0///her2: 0///size: NA///node: NA///age: 58///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M521///dataset: MDA4///series: MDA4///id: 521///er: 0///pgr: 0///her2: 0///size: 3///node: 0///age: 31///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M523///dataset: MDA4///series: MDA4///id: 523///er: 1///pgr: 1///her2: 0///size: 2.3///node: 0///age: 47///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M524///dataset: MDA4///series: MDA4///id: 524///er: 0///pgr: 0///her2: 0///size: 2///node: 0///age: 31///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M536///dataset: MDA4///series: MDA4///id: 536///er: 1///pgr: 0///her2: 0///size: 1.2///node: 1///age: 43///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M557///dataset: MDA4///series: MDA4///id: 557///er: 1///pgr: 0///her2: 0///size: 0///node: 0///age: 57///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M559///dataset: MDA4///series: MDA4///id: 559///er: 1///pgr: 0///her2: 0///size: 2.2///node: 0///age: 55///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M564///dataset: MDA4///series: MDA4///id: 564///er: 0///pgr: 0///her2: 0///size: 1.4///node: 1///age: 68///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M570///dataset: MDA4///series: MDA4///id: 570///er: 0///pgr: 0///her2: 0///size: 3///node: 1///age: 46///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M572///dataset: MDA4///series: MDA4///id: 572///er: 0///pgr: 0///her2: 0///size: 5///node: 1///age: 58///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M578///dataset: MDA4///series: MDA4///id: 578///er: 1///pgr: 1///her2: 0///size: 1.5///node: 1///age: 45///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M595///dataset: MDA4///series: MDA4///id: 595///er: 0///pgr: 0///her2: 1///size: 1.9///node: 1///age: 54///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M604///dataset: MDA4///series: MDA4///id: 604///er: 1///pgr: 1///her2: 0///size: 2.5///node: 1///age: 62///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
     samplename: MDA4_M610///dataset: MDA4///series: MDA4///id: 610///er: 1///pgr: 1///her2: 0///size: 10///node: 1///age: 47///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M631///dataset: MDA4///series: MDA4///id: 631///er: 1///pgr: 0///her2: 0///size: 1///node: 0///age: 32///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M632///dataset: MDA4///series: MDA4///id: 632///er: 1///pgr: 1///her2: 0///size: 2.1///node: 1///age: 48///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M633///dataset: MDA4///series: MDA4///id: 633///er: 0///pgr: 0///her2: 0///size: 0.8///node: 0///age: 46///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M644///dataset: MDA4///series: MDA4///id: 644///er: 1///pgr: 1///her2: 0///size: 3.5///node: 0///age: 53///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M651///dataset: MDA4///series: MDA4///id: 651///er: 1///pgr: 1///her2: 0///size: 0///node: 0///age: 61///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
   samplename: MDA4_M652///dataset: MDA4///series: MDA4///id: 652///er: 1///pgr: 0///her2: 0///size: 0.05///node: 1///age: 66///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M658///dataset: MDA4///series: MDA4///id: 658///er: 0///pgr: 0///her2: 0///size: 0///node: 0///age: 61///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M662///dataset: MDA4///series: MDA4///id: 662///er: 1///pgr: 1///her2: 0///size: 1.3///node: 1///age: 46///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M665///dataset: MDA4///series: MDA4///id: 665///er: 0///pgr: 0///her2: 0///size: 0///node: 0///age: 57///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M675///dataset: MDA4///series: MDA4///id: 675///er: 1///pgr: 0///her2: 0///size: 2///node: 1///age: 64///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M679///dataset: MDA4///series: MDA4///id: 679///er: 0///pgr: 0///her2: 0///size: 4.7///node: 0///age: 72///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M680///dataset: MDA4///series: MDA4///id: 680///er: 1///pgr: 1///her2: 0///size: 2.2///node: 1///age: 34///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M683///dataset: MDA4///series: MDA4///id: 683///er: 1///pgr: 1///her2: 0///size: 0.2///node: 1///age: 59///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M690///dataset: MDA4///series: MDA4///id: 690///er: 1///pgr: 1///her2: 0///size: 0.7///node: 1///age: 73///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M694///dataset: MDA4///series: MDA4///id: 694///er: 0///pgr: 0///her2: 0///size: 1.3///node: 1///age: 38///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M696///dataset: MDA4///series: MDA4///id: 696///er: 0///pgr: 0///her2: 0///size: 0///node: 0///age: 53///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M702///dataset: MDA4///series: MDA4///id: 702///er: 1///pgr: 1///her2: 0///size: 1///node: 0///age: 49///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M703///dataset: MDA4///series: MDA4///id: 703///er: 0///pgr: 0///her2: 0///size: 2.2///node: 0///age: 57///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
    samplename: MDA4_M709///dataset: MDA4///series: MDA4///id: 709///er: NA///pgr: NA///her2: 0///size: 0///node: 0///age: 51///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
   samplename: MDA4_M717///dataset: MDA4///series: MDA4///id: 717///er: 1///pgr: 1///her2: 0///size: 0.21///node: 0///age: 52///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
   samplename: MDA4_M718///dataset: MDA4///series: MDA4///id: 718///er: 0///pgr: 0///her2: 0///size: 0.03///node: 0///age: 36///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
      samplename: MDA4_M728///dataset: MDA4///series: MDA4///id: 728///er: 1///pgr: 1///her2: 0///size: 5///node: 1///age: 49///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
 samplename: MDA4_M729///dataset: MDA4///series: MDA4///id: 729///er: NA///pgr: NA///her2: 0///size: 0.45///node: 1///age: 28///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                           1 
                                                                                                                                                                                                                                     (Other) 
                                                                                                                                                                                                                                          30 

duplicates: 
MDA4.MDA4_M207 MDA4.MDA4_M400           NA's 
             1              1            127 

}}
\source{ http://bioinformatics.mdanderson.org/pubdata.html }
\keyword{datasets}
