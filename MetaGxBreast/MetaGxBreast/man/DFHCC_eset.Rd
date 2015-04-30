\name{ DFHCC_eset }
\alias{ DFHCC_eset }
\docType{data}
\title{  }
\usage{data( DFHCC_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2826790/ 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19615 
  PMIDs: 20098429 
  No abstract available.
  notes:
   summary:     
      A small number of over-expressed and over-amplified genes were significantly associated with early recurrence despite adj
uvant therapy.  This was verified in independent cohorts. 
   version:     
      2015-04-27 18:43:59
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
assayData: 42447 features, 115 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
DFHCC_AR2004020603 DFHCC_AR2004020604 DFHCC_AR2004020606 DFHCC_AR2004020607 DFHCC_AR2004020608 DFHCC_AR2004020609 
                 1                  1                  1                  1                  1                  1 
DFHCC_AR2004020610 DFHCC_AR2004020611 DFHCC_AR2004020612 DFHCC_AR2004020615 DFHCC_AR2004020618 DFHCC_AR2004020619 
                 1                  1                  1                  1                  1                  1 
DFHCC_AR2004020621 DFHCC_AR2004020622 DFHCC_AR2004020623 DFHCC_AR2004020624 DFHCC_AR2004020625 DFHCC_AR2004020626 
                 1                  1                  1                  1                  1                  1 
DFHCC_AR2004020628 DFHCC_AR2004020629 DFHCC_AR2004020631 DFHCC_AR2004020632 DFHCC_AR2004020633 DFHCC_AR2004020634 
                 1                  1                  1                  1                  1                  1 
DFHCC_AR2004020635 DFHCC_AR2004020636 DFHCC_AR2004020639 DFHCC_AR2004020640 DFHCC_AR2004020642 DFHCC_AR2004021101 
                 1                  1                  1                  1                  1                  1 
DFHCC_AR2004021102 DFHCC_AR2004021103 DFHCC_AR2004021104 DFHCC_AR2004021105 DFHCC_AR2004021107 DFHCC_AR2004021108 
                 1                  1                  1                  1                  1                  1 
DFHCC_AR2004021110 DFHCC_AR2004021111 DFHCC_AR2004021112 DFHCC_AR2004021113 DFHCC_AR2004021114 DFHCC_AR2004021115 
                 1                  1                  1                  1                  1                  1 
DFHCC_AR2004021117 DFHCC_AR2004021118 DFHCC_AR2004021121 DFHCC_AR2004021122 DFHCC_AR2004021123 DFHCC_AR2004021124 
                 1                  1                  1                  1                  1                  1 
DFHCC_AR2004021125 DFHCC_AR2004021126 DFHCC_AR2004021129 DFHCC_AR2004021130 DFHCC_AR2004021131 DFHCC_AR2004021132 
                 1                  1                  1                  1                  1                  1 
DFHCC_AR2004021133 DFHCC_AR2004021134 DFHCC_AR2004021135 DFHCC_AR2004021137 DFHCC_AR2004021138 DFHCC_AR2004021139 
                 1                  1                  1                  1                  1                  1 
DFHCC_AR2004021141 DFHCC_AR2004021143 DFHCC_AR2004021144 DFHCC_AR2004021145 DFHCC_AR2004021147 DFHCC_AR2004021252 
                 1                  1                  1                  1                  1                  1 
DFHCC_AR2004021253 DFHCC_AR2004021254 DFHCC_AR2004021255 DFHCC_AR2004021257 DFHCC_AR2004021263 DFHCC_AR2004021265 
                 1                  1                  1                  1                  1                  1 
DFHCC_AR2004021267 DFHCC_AR2004021271 DFHCC_AR2004021272 DFHCC_AR2004021273 DFHCC_AR2004021274 DFHCC_AR2004021275 
                 1                  1                  1                  1                  1                  1 
DFHCC_AR2004021276 DFHCC_AR2004021277 DFHCC_AR2004021278 DFHCC_AR2004021279 DFHCC_AR2004021281 DFHCC_AR2004021285 
                 1                  1                  1                  1                  1                  1 
DFHCC_AR2004021286 DFHCC_AR2004021288 DFHCC_AR2004021289 DFHCC_AR2004021290 DFHCC_AR2004021291 DFHCC_AR2004021292 
                 1                  1                  1                  1                  1                  1 
DFHCC_AR2004021293 DFHCC_AR2004021294 DFHCC_AR2004021295 DFHCC_AR2004021296 DFHCC_AR2004021297 DFHCC_AR2004031202 
                 1                  1                  1                  1                  1                  1 
DFHCC_AR2004031204 DFHCC_AR2004031208 DFHCC_AR2004031209            (Other) 
                 1                  1                  1                 16 

alt_sample_name: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    6.0   155.0   230.0   293.3   398.5   828.0 

sample_type: 
tumor 
  115 

er: 
negative positive 
      45       70 

pgr: 
negative positive 
      51       64 

her2: 
negative positive 
      79       36 

tumor_size: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  0.800   1.350   2.100   2.312   2.850   6.500 

N: 
 0  1 
62 53 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  32.00   45.00   53.00   53.89   60.00   85.00 

grade: 
 1  2  3 
23 28 64 

dmfs_days: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
     30    1500    1920    1799    2325    2640 

dmfs_status: 
deceased_or_recurrence    living_norecurrence 
                    14                    101 

treatment: 
chemo.plus.hormono       chemotherapy     hormonotherapy          untreated               NA's 
                42                 38                 22                  7                  6 

batch: 
DFHCC 
  115 

uncurated_author_metadata: 
  samplename: DFHCC_AR2004020603///dataset: DFHCC///series: DFHCC///id: 292///er: 0///pgr: 0///her2: 1///size: 2.3///node: 1///age: 53///grade: 3///t.dmfs: 720///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004020604///dataset: DFHCC///series: DFHCC///id: 310///er: 0///pgr: 0///her2: 0///size: 0.9///node: 0///age: 60///grade: 3///t.dmfs: 1140///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004020606///dataset: DFHCC///series: DFHCC///id: 316///er: 0///pgr: 0///her2: 0///size: 3.5///node: 0///age: 65///grade: 3///t.dmfs: 1800///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
  samplename: DFHCC_AR2004020607///dataset: DFHCC///series: DFHCC///id: 325///er: 0///pgr: 0///her2: 1///size: 1.8///node: 1///age: 48///grade: 3///t.dmfs: 270///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004020608///dataset: DFHCC///series: DFHCC///id: 330///er: 0///pgr: 0///her2: 0///size: 2.1///node: 1///age: 58///grade: 3///t.dmfs: 1860///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
  samplename: DFHCC_AR2004020609///dataset: DFHCC///series: DFHCC///id: 343///er: 0///pgr: 0///her2: 1///size: 2.5///node: 1///age: 45///grade: 3///t.dmfs: 480///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004020610///dataset: DFHCC///series: DFHCC///id: 390///er: 0///pgr: 0///her2: 0///size: 1.2///node: 1///age: 49///grade: 3///t.dmfs: 1920///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004020611///dataset: DFHCC///series: DFHCC///id: 397///er: 0///pgr: 0///her2: 1///size: 2.1///node: 0///age: 56///grade: 3///t.dmfs: 2010///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
  samplename: DFHCC_AR2004020612///dataset: DFHCC///series: DFHCC///id: 400///er: 0///pgr: 0///her2: 1///size: 2.5///node: 1///age: 85///grade: 2///t.dmfs: 810///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004020615///dataset: DFHCC///series: DFHCC///id: 451///er: 0///pgr: 0///her2: 0///size: 2.5///node: 0///age: 66///grade: 3///t.dmfs: 1860///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
samplename: DFHCC_AR2004020618///dataset: DFHCC///series: DFHCC///id: 488///er: 0///pgr: 0///her2: 0///size: 1.8///node: 0///age: 58///grade: 2///t.dmfs: 1560///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                   1 
     samplename: DFHCC_AR2004020619///dataset: DFHCC///series: DFHCC///id: 498///er: 0///pgr: 0///her2: 0///size: 5///node: 0///age: 79///grade: 3///t.dmfs: 30///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004020621///dataset: DFHCC///series: DFHCC///id: 534///er: 0///pgr: 0///her2: 0///size: 2.1///node: 0///age: 58///grade: 3///t.dmfs: 1590///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
  samplename: DFHCC_AR2004020622///dataset: DFHCC///series: DFHCC///id: 538///er: 0///pgr: 0///her2: 0///size: 1///node: 0///age: 48///grade: 3///t.dmfs: 1620///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004020623///dataset: DFHCC///series: DFHCC///id: 554///er: 0///pgr: 0///her2: 0///size: 0.9///node: 0///age: 50///grade: 3///t.dmfs: 1530///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004020624///dataset: DFHCC///series: DFHCC///id: 558///er: 0///pgr: 0///her2: 0///size: 2.2///node: 0///age: 40///grade: 3///t.dmfs: 1440///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
   samplename: DFHCC_AR2004020625///dataset: DFHCC///series: DFHCC///id: 582///er: 0///pgr: 0///her2: 0///size: 1///node: 0///age: 58///grade: 3///t.dmfs: 1560///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004020626///dataset: DFHCC///series: DFHCC///id: 588///er: 0///pgr: 0///her2: 1///size: 3.5///node: 1///age: 63///grade: 3///t.dmfs: 1500///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004020628///dataset: DFHCC///series: DFHCC///id: 591///er: 0///pgr: 0///her2: 0///size: 2.7///node: 0///age: 41///grade: 3///t.dmfs: 1500///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004020629///dataset: DFHCC///series: DFHCC///id: 598///er: 0///pgr: 0///her2: 0///size: 2.5///node: 0///age: 74///grade: 2///t.dmfs: 1200///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                   1 
   samplename: DFHCC_AR2004020631///dataset: DFHCC///series: DFHCC///id: 609///er: 0///pgr: 0///her2: 0///size: 3///node: 0///age: 48///grade: 3///t.dmfs: 1470///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004020632///dataset: DFHCC///series: DFHCC///id: 613///er: 0///pgr: 0///her2: 0///size: 1.3///node: 0///age: 78///grade: 3///t.dmfs: 1590///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004020633///dataset: DFHCC///series: DFHCC///id: 653///er: 0///pgr: 0///her2: 0///size: 1.4///node: 0///age: 55///grade: 3///t.dmfs: 1500///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004020634///dataset: DFHCC///series: DFHCC///id: 664///er: 0///pgr: 0///her2: 0///size: 1.7///node: 0///age: 43///grade: 3///t.dmfs: 1050///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004020635///dataset: DFHCC///series: DFHCC///id: 665///er: 0///pgr: 0///her2: 0///size: 2.3///node: 0///age: 72///grade: 3///t.dmfs: 1500///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004020636///dataset: DFHCC///series: DFHCC///id: 670///er: 0///pgr: 0///her2: 1///size: 3.4///node: 1///age: 49///grade: 3///t.dmfs: 1080///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004020639///dataset: DFHCC///series: DFHCC///id: 799///er: 0///pgr: 0///her2: 0///size: 2.8///node: 0///age: 43///grade: 3///t.dmfs: 1080///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004020640///dataset: DFHCC///series: DFHCC///id: 828///er: 0///pgr: 0///her2: 0///size: 2.2///node: 0///age: 56///grade: 3///t.dmfs: 1170///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004020642///dataset: DFHCC///series: DFHCC///id: 396///er: 1///pgr: 1///her2: 1///size: 1.1///node: 1///age: 41///grade: 3///t.dmfs: 1950///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                   1 
  samplename: DFHCC_AR2004021101///dataset: DFHCC///series: DFHCC///id: 92///er: 1///pgr: 1///her2: 0///size: 3.5///node: 0///age: 62///grade: 1///t.dmfs: 2400///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                   1 
  samplename: DFHCC_AR2004021102///dataset: DFHCC///series: DFHCC///id: 129///er: 1///pgr: 1///her2: 1///size: 1///node: 0///age: 43///grade: 3///t.dmfs: 2340///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004021103///dataset: DFHCC///series: DFHCC///id: 131///er: 1///pgr: 1///her2: 1///size: 1.7///node: 0///age: 54///grade: 1///t.dmfs: 2340///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                   1 
  samplename: DFHCC_AR2004021104///dataset: DFHCC///series: DFHCC///id: 96///er: 1///pgr: 1///her2: 1///size: 1.4///node: 1///age: 55///grade: 3///t.dmfs: 1530///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                   1 
   samplename: DFHCC_AR2004021105///dataset: DFHCC///series: DFHCC///id: 86///er: 1///pgr: 1///her2: 0///size: 2.6///node: 1///age: 46///grade: 2///t.dmfs: 180///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                   1 
  samplename: DFHCC_AR2004021107///dataset: DFHCC///series: DFHCC///id: 73///er: 1///pgr: 1///her2: 0///size: 0.9///node: 0///age: 46///grade: 2///t.dmfs: 2550///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004021108///dataset: DFHCC///series: DFHCC///id: 212///er: 1///pgr: 1///her2: 0///size: 0.8///node: 0///age: 44///grade: 1///t.dmfs: 2190///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                   1 
  samplename: DFHCC_AR2004021110///dataset: DFHCC///series: DFHCC///id: 74///er: 1///pgr: 1///her2: 1///size: 2.2///node: 0///age: 54///grade: 2///t.dmfs: 2460///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004021111///dataset: DFHCC///series: DFHCC///id: 196///er: 0///pgr: 0///her2: 1///size: 1.8///node: 0///age: 49///grade: 3///t.dmfs: 2250///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004021112///dataset: DFHCC///series: DFHCC///id: 123///er: 1///pgr: 1///her2: 0///size: 2.2///node: 1///age: 84///grade: 2///t.dmfs: 2370///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004021113///dataset: DFHCC///series: DFHCC///id: 180///er: 1///pgr: 0///her2: 1///size: 4.2///node: 0///age: 43///grade: 3///t.dmfs: 1560///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004021114///dataset: DFHCC///series: DFHCC///id: 249///er: 0///pgr: 0///her2: 1///size: 2.5///node: 0///age: 54///grade: 2///t.dmfs: 2430///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004021115///dataset: DFHCC///series: DFHCC///id: 207///er: 1///pgr: 1///her2: 0///size: 2.3///node: 0///age: 45///grade: 2///t.dmfs: 1770///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                   1 
  samplename: DFHCC_AR2004021117///dataset: DFHCC///series: DFHCC///id: 65///er: 1///pgr: 1///her2: 1///size: 1.7///node: 1///age: 45///grade: 2///t.dmfs: 2490///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004021118///dataset: DFHCC///series: DFHCC///id: 157///er: 1///pgr: 0///her2: 0///size: 1.1///node: 1///age: 35///grade: 2///t.dmfs: 2550///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                   1 
    samplename: DFHCC_AR2004021121///dataset: DFHCC///series: DFHCC///id: 148///er: 0///pgr: 0///her2: 0///size: 3///node: 0///age: 53///grade: 3///t.dmfs: 540///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004021122///dataset: DFHCC///series: DFHCC///id: 228///er: 1///pgr: 1///her2: 0///size: 2.3///node: 0///age: 60///grade: 1///t.dmfs: 2100///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004021123///dataset: DFHCC///series: DFHCC///id: 162///er: 1///pgr: 1///her2: 0///size: 1.5///node: 1///age: 32///grade: 2///t.dmfs: 2490///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                   1 
   samplename: DFHCC_AR2004021124///dataset: DFHCC///series: DFHCC///id: 175///er: 1///pgr: 0///her2: 0///size: 4///node: 1///age: 40///grade: 1///t.dmfs: 2220///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                   1 
  samplename: DFHCC_AR2004021125///dataset: DFHCC///series: DFHCC///id: 14///er: 0///pgr: 0///her2: 1///size: 1.5///node: 0///age: 45///grade: 3///t.dmfs: 2610///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004021126///dataset: DFHCC///series: DFHCC///id: 246///er: 1///pgr: 1///her2: 0///size: 1.1///node: 0///age: 46///grade: 1///t.dmfs: 2580///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
    samplename: DFHCC_AR2004021129///dataset: DFHCC///series: DFHCC///id: 60///er: 1///pgr: 1///her2: 0///size: 3///node: 1///age: 41///grade: 1///t.dmfs: 2580///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004021130///dataset: DFHCC///series: DFHCC///id: 138///er: 1///pgr: 1///her2: 1///size: 2.3///node: 1///age: 46///grade: 3///t.dmfs: 2400///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                   1 
   samplename: DFHCC_AR2004021131///dataset: DFHCC///series: DFHCC///id: 204///er: 1///pgr: 1///her2: 0///size: 1///node: 0///age: 36///grade: 2///t.dmfs: 2340///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                   1 
samplename: DFHCC_AR2004021132///dataset: DFHCC///series: DFHCC///id: 239///er: 0///pgr: 0///her2: 0///size: 1.1///node: 0///age: 68///grade: 3///t.dmfs: 2190///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004021133///dataset: DFHCC///series: DFHCC///id: 253///er: 0///pgr: 0///her2: 0///size: 1.1///node: 0///age: 54///grade: 3///t.dmfs: 1680///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
  samplename: DFHCC_AR2004021134///dataset: DFHCC///series: DFHCC///id: 11///er: 1///pgr: 1///her2: 0///size: 2.5///node: 1///age: 81///grade: 1///t.dmfs: 2100///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                   1 
    samplename: DFHCC_AR2004021135///dataset: DFHCC///series: DFHCC///id: 21///er: 0///pgr: 0///her2: 1///size: 3///node: 1///age: 45///grade: 3///t.dmfs: 2490///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
  samplename: DFHCC_AR2004021137///dataset: DFHCC///series: DFHCC///id: 214///er: 0///pgr: 0///her2: 0///size: 5.2///node: 1///age: 53///grade: 3///t.dmfs: 420///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
  samplename: DFHCC_AR2004021138///dataset: DFHCC///series: DFHCC///id: 12///er: 0///pgr: 0///her2: 0///size: 2.5///node: 0///age: 48///grade: 3///t.dmfs: 2250///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004021139///dataset: DFHCC///series: DFHCC///id: 165///er: 1///pgr: 0///her2: 0///size: 2.5///node: 0///age: 53///grade: 2///t.dmfs: 2310///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004021141///dataset: DFHCC///series: DFHCC///id: 126///er: 1///pgr: 1///her2: 1///size: 1.7///node: 1///age: 57///grade: 3///t.dmfs: 2190///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                   1 
   samplename: DFHCC_AR2004021143///dataset: DFHCC///series: DFHCC///id: 6///er: 1///pgr: 1///her2: 0///size: 2.1///node: 0///age: 57///grade: 2///t.dmfs: 1770///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                   1 
  samplename: DFHCC_AR2004021144///dataset: DFHCC///series: DFHCC///id: 58///er: 1///pgr: 1///her2: 0///size: 1.5///node: 1///age: 46///grade: 3///t.dmfs: 2550///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004021145///dataset: DFHCC///series: DFHCC///id: 164///er: 0///pgr: 0///her2: 1///size: 1.9///node: 0///age: 50///grade: 3///t.dmfs: 2310///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004021147///dataset: DFHCC///series: DFHCC///id: 119///er: 1///pgr: 1///her2: 0///size: 3.5///node: 1///age: 58///grade: 2///t.dmfs: 2430///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004021252///dataset: DFHCC///series: DFHCC///id: 154///er: 1///pgr: 1///her2: 0///size: 1.5///node: 0///age: 62///grade: 1///t.dmfs: 2430///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                   1 
   samplename: DFHCC_AR2004021253///dataset: DFHCC///series: DFHCC///id: 182///er: 0///pgr: 0///her2: 0///size: 2///node: 0///age: 44///grade: 3///t.dmfs: 1710///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
  samplename: DFHCC_AR2004021254///dataset: DFHCC///series: DFHCC///id: 35///er: 1///pgr: 1///her2: 1///size: 1.3///node: 1///age: 37///grade: 2///t.dmfs: 1440///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004021255///dataset: DFHCC///series: DFHCC///id: 151///er: 0///pgr: 0///her2: 1///size: 2.5///node: 1///age: 55///grade: 3///t.dmfs: 1740///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
  samplename: DFHCC_AR2004021257///dataset: DFHCC///series: DFHCC///id: 203///er: 0///pgr: 0///her2: 0///size: 1.7///node: 0///age: 64///grade: 3///t.dmfs: 630///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                   1 
   samplename: DFHCC_AR2004021263///dataset: DFHCC///series: DFHCC///id: 187///er: 0///pgr: 0///her2: 1///size: 3///node: 0///age: 53///grade: 3///t.dmfs: 2310///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004021265///dataset: DFHCC///series: DFHCC///id: 237///er: 1///pgr: 1///her2: 0///size: 1.1///node: 1///age: 66///grade: 1///t.dmfs: 2460///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004021267///dataset: DFHCC///series: DFHCC///id: 238///er: 1///pgr: 1///her2: 1///size: 1.2///node: 0///age: 63///grade: 3///t.dmfs: 2130///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004021271///dataset: DFHCC///series: DFHCC///id: 252///er: 1///pgr: 1///her2: 0///size: 1.2///node: 1///age: 66///grade: 1///t.dmfs: 1650///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                   1 
   samplename: DFHCC_AR2004021272///dataset: DFHCC///series: DFHCC///id: 206///er: 1///pgr: 1///her2: 0///size: 4.2///node: 1///age: 38///grade: 3///t.dmfs: 60///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                                   1 
   samplename: DFHCC_AR2004021273///dataset: DFHCC///series: DFHCC///id: 230///er: 1///pgr: 1///her2: 0///size: 3///node: 1///age: 46///grade: 3///t.dmfs: 2280///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004021274///dataset: DFHCC///series: DFHCC///id: 142///er: 1///pgr: 1///her2: 0///size: 2.2///node: 1///age: 60///grade: 3///t.dmfs: 2640///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004021275///dataset: DFHCC///series: DFHCC///id: 168///er: 1///pgr: 1///her2: 0///size: 1.3///node: 1///age: 50///grade: 1///t.dmfs: 1410///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004021276///dataset: DFHCC///series: DFHCC///id: 222///er: 1///pgr: 1///her2: 0///size: 1.2///node: 1///age: 37///grade: 1///t.dmfs: 2490///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                   1 
   samplename: DFHCC_AR2004021277///dataset: DFHCC///series: DFHCC///id: 219///er: 1///pgr: 1///her2: 1///size: 2///node: 0///age: 53///grade: 3///t.dmfs: 1620///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                   1 
  samplename: DFHCC_AR2004021278///dataset: DFHCC///series: DFHCC///id: 91///er: 1///pgr: 0///her2: 0///size: 1.1///node: 0///age: 57///grade: 1///t.dmfs: 2370///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                   1 
  samplename: DFHCC_AR2004021279///dataset: DFHCC///series: DFHCC///id: 156///er: 0///pgr: 0///her2: 1///size: 2.9///node: 1///age: 57///grade: 2///t.dmfs: 450///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
  samplename: DFHCC_AR2004021281///dataset: DFHCC///series: DFHCC///id: 227///er: 1///pgr: 1///her2: 0///size: 2.5///node: 1///age: 84///grade: 3///t.dmfs: 540///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004021285///dataset: DFHCC///series: DFHCC///id: 177///er: 1///pgr: 1///her2: 0///size: 1.2///node: 0///age: 46///grade: 2///t.dmfs: 2280///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004021286///dataset: DFHCC///series: DFHCC///id: 231///er: 1///pgr: 1///her2: 0///size: 1.5///node: 1///age: 45///grade: 1///t.dmfs: 2430///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                   1 
   samplename: DFHCC_AR2004021288///dataset: DFHCC///series: DFHCC///id: 150///er: 1///pgr: 1///her2: 0///size: 5///node: 0///age: 62///grade: 1///t.dmfs: 2490///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004021289///dataset: DFHCC///series: DFHCC///id: 221///er: 1///pgr: 1///her2: 0///size: 1.5///node: 1///age: 58///grade: 1///t.dmfs: 1620///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004021290///dataset: DFHCC///series: DFHCC///id: 174///er: 1///pgr: 1///her2: 0///size: 1.4///node: 0///age: 50///grade: 1///t.dmfs: 2370///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004021291///dataset: DFHCC///series: DFHCC///id: 178///er: 1///pgr: 1///her2: 0///size: 1.1///node: 0///age: 74///grade: 1///t.dmfs: 2220///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004021292///dataset: DFHCC///series: DFHCC///id: 192///er: 1///pgr: 1///her2: 0///size: 1.1///node: 1///age: 57///grade: 3///t.dmfs: 2220///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                   1 
  samplename: DFHCC_AR2004021293///dataset: DFHCC///series: DFHCC///id: 87///er: 1///pgr: 0///her2: 0///size: 0.8///node: 0///age: 58///grade: 1///t.dmfs: 2430///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                   1 
     samplename: DFHCC_AR2004021294///dataset: DFHCC///series: DFHCC///id: 7///er: 1///pgr: 1///her2: 0///size: 1///node: 0///age: 59///grade: 1///t.dmfs: 2430///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004021295///dataset: DFHCC///series: DFHCC///id: 229///er: 1///pgr: 1///her2: 0///size: 2.1///node: 0///age: 46///grade: 2///t.dmfs: 2280///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                   1 
   samplename: DFHCC_AR2004021296///dataset: DFHCC///series: DFHCC///id: 9///er: 1///pgr: 1///her2: 1///size: 1.5///node: 0///age: 54///grade: 2///t.dmfs: 1950///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004021297///dataset: DFHCC///series: DFHCC///id: 224///er: 1///pgr: 1///her2: 0///size: 1.3///node: 0///age: 64///grade: 1///t.dmfs: 2490///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004031202///dataset: DFHCC///series: DFHCC///id: 666///er: 1///pgr: 1///her2: 1///size: 2.1///node: 0///age: 64///grade: 3///t.dmfs: 1440///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004031204///dataset: DFHCC///series: DFHCC///id: 643///er: 0///pgr: 0///her2: 1///size: 5.5///node: 1///age: 68///grade: 3///t.dmfs: 1530///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                                   1 
   samplename: DFHCC_AR2004031208///dataset: DFHCC///series: DFHCC///id: 552///er: 1///pgr: 1///her2: 0///size: 2///node: 1///age: 49///grade: 1///t.dmfs: 1680///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                   1 
 samplename: DFHCC_AR2004031209///dataset: DFHCC///series: DFHCC///id: 489///er: 1///pgr: 1///her2: 0///size: 6.5///node: 1///age: 37///grade: 2///t.dmfs: 1770///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 6 
                                                                                                                                                                                                                                                   1 
                                                                                                                                                                                                                                             (Other) 
                                                                                                                                                                                                                                                  16 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19615 }
\keyword{datasets}
