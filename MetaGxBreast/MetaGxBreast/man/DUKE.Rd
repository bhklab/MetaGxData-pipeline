\name{ DUKE_eset }
\alias{ DUKE_eset }
\docType{data}
\title{ DUKE_eset }
\usage{data( DUKE_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/16273092 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse3143 
  PMIDs: 16273092 
  No abstract available.
  notes:
   summary:
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: 1000_at 1001_at ... AFFX-MurIL4_at (12085 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 12085 features, 171 samples
Platform type: 
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

   1 observation deleted due to missingness 
records   n.max n.start  events  median 0.95LCL 0.95UCL 
 170.00  170.00  170.00   43.00    9.01    6.22      NA 

--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
 DUKE_KF-067  DUKE_KF-069  DUKE_KF-070  DUKE_KF-075  DUKE_KF-076  DUKE_KF-077  DUKE_KF-079  DUKE_KF-083  DUKE_KF-090 
           1            1            1            1            1            1            1            1            1 
 DUKE_KF-094  DUKE_KF-095  DUKE_KF-097  DUKE_KF-098  DUKE_KF-101  DUKE_KF-104  DUKE_KF-105  DUKE_KF-106  DUKE_KF-109 
           1            1            1            1            1            1            1            1            1 
 DUKE_KF-110  DUKE_KF-113  DUKE_KF-114  DUKE_KF-115  DUKE_KF-116  DUKE_KF-117  DUKE_KF-118  DUKE_KF-119  DUKE_KF-120 
           1            1            1            1            1            1            1            1            1 
 DUKE_KF-123  DUKE_KF-124  DUKE_KF-125  DUKE_KF-126  DUKE_KF-127  DUKE_KF-128  DUKE_KF-129  DUKE_KF-130  DUKE_KF-131 
           1            1            1            1            1            1            1            1            1 
 DUKE_KF-132  DUKE_KF-133  DUKE_KF-134  DUKE_KF-135  DUKE_KF-136  DUKE_KF-137  DUKE_KF-138  DUKE_KF-139  DUKE_KF-140 
           1            1            1            1            1            1            1            1            1 
 DUKE_KF-141  DUKE_KF-142  DUKE_KF-143  DUKE_KF-144  DUKE_KF-145  DUKE_KF-146  DUKE_KF-147  DUKE_KF-149  DUKE_KF-150 
           1            1            1            1            1            1            1            1            1 
 DUKE_KF-151  DUKE_KF-152  DUKE_KF-153  DUKE_KF-154  DUKE_KF-155  DUKE_KF-157  DUKE_KF-159  DUKE_KF-160  DUKE_KF-176 
           1            1            1            1            1            1            1            1            1 
 DUKE_KF-178  DUKE_KF-182  DUKE_KF-183  DUKE_KF-184  DUKE_KF-187  DUKE_KF-188  DUKE_KF-189  DUKE_KF-193  DUKE_KF-194 
           1            1            1            1            1            1            1            1            1 
 DUKE_KF-195  DUKE_KF-196  DUKE_KF-200  DUKE_KF-201  DUKE_KF-203  DUKE_KF-204  DUKE_KF-207 DUKE_T00-050 DUKE_T00-064 
           1            1            1            1            1            1            1            1            1 
DUKE_T00-070 DUKE_T00-109 DUKE_T00-214 DUKE_T00-454 DUKE_T00-622 DUKE_T01-047 DUKE_T01-052 DUKE_T01-145 DUKE_T01-154 
           1            1            1            1            1            1            1            1            1 
DUKE_T01-156 DUKE_T01-158 DUKE_T01-166 DUKE_T01-194 DUKE_T01-202 DUKE_T01-204 DUKE_T01-255 DUKE_T01-310 DUKE_T01-344 
           1            1            1            1            1            1            1            1            1 
     (Other) 
          72 

alt_sample_name: 
 KF-067  KF-069  KF-070  KF-075  KF-076  KF-077  KF-079  KF-083  KF-090  KF-094  KF-095  KF-097  KF-098  KF-101  KF-104 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
 KF-105  KF-106  KF-109  KF-110  KF-113  KF-114  KF-115  KF-116  KF-117  KF-118  KF-119  KF-120  KF-123  KF-124  KF-125 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
 KF-126  KF-127  KF-128  KF-129  KF-130  KF-131  KF-132  KF-133  KF-134  KF-135  KF-136  KF-137  KF-138  KF-139  KF-140 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
 KF-141  KF-142  KF-143  KF-144  KF-145  KF-146  KF-147  KF-149  KF-150  KF-151  KF-152  KF-153  KF-154  KF-155  KF-157 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
 KF-159  KF-160  KF-176  KF-178  KF-182  KF-183  KF-184  KF-187  KF-188  KF-189  KF-193  KF-194  KF-195  KF-196  KF-200 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
 KF-201  KF-203  KF-204  KF-207 T00-050 T00-064 T00-070 T00-109 T00-214 T00-454 T00-622 T01-047 T01-052 T01-145 T01-154 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
T01-156 T01-158 T01-166 T01-194 T01-202 T01-204 T01-255 T01-310 T01-344 (Other) 
      1       1       1       1       1       1       1       1       1      72 

sample_type: 
tumor 
  171 

er: 
negative positive 
      57      114 

pgr: 
negative positive     NA's 
      23       65       83 

tumor_size: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
   0.20    1.80    2.30    2.74    3.50    8.50      83 

N: 
   0    1 NA's 
  53   36   82 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  171.0   417.0   957.5  1235.0  1852.0  4069.0       1 

vital_status: 
deceased   living     NA's 
      43      127        1 

batch: 
DUKE 
 171 

uncurated_author_metadata: 
   samplename: DUKE_KF-067///dataset: DUKE///series: DUKE///id: KF-067///er: 1///pgr: 1///her2: NA///size: 7.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 629///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
    samplename: DUKE_KF-069///dataset: DUKE///series: DUKE///id: KF-069///er: 1///pgr: 1///her2: NA///size: 4///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1852///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
 samplename: DUKE_KF-070///dataset: DUKE///series: DUKE///id: KF-070///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1915///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
 samplename: DUKE_KF-075///dataset: DUKE///series: DUKE///id: KF-075///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1832///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-076///dataset: DUKE///series: DUKE///id: KF-076///er: 1///pgr: 1///her2: NA///size: 1.7///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1850///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
    samplename: DUKE_KF-077///dataset: DUKE///series: DUKE///id: KF-077///er: 1///pgr: 1///her2: NA///size: 2///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1769///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
    samplename: DUKE_KF-079///dataset: DUKE///series: DUKE///id: KF-079///er: 1///pgr: 1///her2: NA///size: 2///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1600///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
 samplename: DUKE_KF-083///dataset: DUKE///series: DUKE///id: KF-083///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1529///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
 samplename: DUKE_KF-090///dataset: DUKE///series: DUKE///id: KF-090///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1422///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-094///dataset: DUKE///series: DUKE///id: KF-094///er: 1///pgr: 1///her2: NA///size: 1.8///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1462///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-095///dataset: DUKE///series: DUKE///id: KF-095///er: 1///pgr: 1///her2: NA///size: 0.6///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1476///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
   samplename: DUKE_KF-097///dataset: DUKE///series: DUKE///id: KF-097///er: 1///pgr: 1///her2: NA///size: 2.2///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 964///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
 samplename: DUKE_KF-098///dataset: DUKE///series: DUKE///id: KF-098///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1324///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-101///dataset: DUKE///series: DUKE///id: KF-101///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 210///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-104///dataset: DUKE///series: DUKE///id: KF-104///er: 1///pgr: 1///her2: NA///size: 1.2///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1247///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-105///dataset: DUKE///series: DUKE///id: KF-105///er: 1///pgr: 1///her2: NA///size: 1.3///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1112///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                           [1] "DUKE2"
[1] "EMC2"
[1] "EORTC10994"
[1] "EXPO"
[1] "FNCLCC"
[1] "HLP"
[1] "IRB"
[1] "KOO"
[1] "LUND"
[1] "LUND2"
[1] "MAINZ"
[1] "MAQC2"
[1] "MCCC"
[1] "MDA4"
[1] "METABRIC"
[1] "MSK"
[1] "MUG"
[1] "NCCS"
[1] "NCI"
[1] "NKI"
[1] "PNC"
[1] "STK"
[1] "STNO2"
[1] "SUPERTAM_HGU133A"
[1] "SUPERTAM_HGU133PLUS2"
[1] "TCGA"
[1] "TRANSBIG"
[1] "UCSF"
[1] "UNC4"
[1] "UNT"
[1] "UPP"
[1] "VDX"
[1] NA
                                                                                                                                                                                                                  1 
  samplename: DUKE_KF-109///dataset: DUKE///series: DUKE///id: KF-109///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
    samplename: DUKE_KF-110///dataset: DUKE///series: DUKE///id: KF-110///er: 1///pgr: 1///her2: NA///size: 2///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1100///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-113///dataset: DUKE///series: DUKE///id: KF-113///er: 1///pgr: 1///her2: NA///size: 1.5///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1380///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
 samplename: DUKE_KF-114///dataset: DUKE///series: DUKE///id: KF-114///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1240///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
 samplename: DUKE_KF-115///dataset: DUKE///series: DUKE///id: KF-115///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1113///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-116///dataset: DUKE///series: DUKE///id: KF-116///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 655///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-117///dataset: DUKE///series: DUKE///id: KF-117///er: 1///pgr: 1///her2: NA///size: 1.4///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1138///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-118///dataset: DUKE///series: DUKE///id: KF-118///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 358///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-119///dataset: DUKE///series: DUKE///id: KF-119///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 877///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-120///dataset: DUKE///series: DUKE///id: KF-120///er: 1///pgr: 1///her2: NA///size: 1.5///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1331///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-123///dataset: DUKE///series: DUKE///id: KF-123///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 477///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
 samplename: DUKE_KF-124///dataset: DUKE///series: DUKE///id: KF-124///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1083///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
    samplename: DUKE_KF-125///dataset: DUKE///series: DUKE///id: KF-125///er: 1///pgr: 1///her2: NA///size: 4///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1815///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-126///dataset: DUKE///series: DUKE///id: KF-126///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 493///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
 samplename: DUKE_KF-127///dataset: DUKE///series: DUKE///id: KF-127///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 2078///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-128///dataset: DUKE///series: DUKE///id: KF-128///er: 1///pgr: 1///her2: NA///size: 2.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 2905///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
    samplename: DUKE_KF-129///dataset: DUKE///series: DUKE///id: KF-129///er: 1///pgr: 0///her2: NA///size: 2///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 2501///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-130///dataset: DUKE///series: DUKE///id: KF-130///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 900///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
    samplename: DUKE_KF-131///dataset: DUKE///series: DUKE///id: KF-131///er: 1///pgr: 1///her2: NA///size: 3///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1685///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-132///dataset: DUKE///series: DUKE///id: KF-132///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 546///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-133///dataset: DUKE///series: DUKE///id: KF-133///er: 1///pgr: 1///her2: NA///size: 0.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1498///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-134///dataset: DUKE///series: DUKE///id: KF-134///er: 1///pgr: NA///her2: NA///size: NA///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1665///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-135///dataset: DUKE///series: DUKE///id: KF-135///er: 0///pgr: 0///her2: NA///size: 2.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1043///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-136///dataset: DUKE///series: DUKE///id: KF-136///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 670///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-137///dataset: DUKE///series: DUKE///id: KF-137///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 365///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
    samplename: DUKE_KF-138///dataset: DUKE///series: DUKE///id: KF-138///er: 1///pgr: 1///her2: NA///size: 4///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 2271///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-139///dataset: DUKE///series: DUKE///id: KF-139///er: 1///pgr: 1///her2: NA///size: 3.1///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1847///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-140///dataset: DUKE///series: DUKE///id: KF-140///er: 0///pgr: 0///her2: NA///size: 3.7///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 2731///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
 samplename: DUKE_KF-141///dataset: DUKE///series: DUKE///id: KF-141///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1226///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
   samplename: DUKE_KF-142///dataset: DUKE///series: DUKE///id: KF-142///er: 1///pgr: 0///her2: NA///size: 3.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 280///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
 samplename: DUKE_KF-143///dataset: DUKE///series: DUKE///id: KF-143///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1378///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-144///dataset: DUKE///series: DUKE///id: KF-144///er: 1///pgr: 1///her2: NA///size: 0.5///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 2542///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
     samplename: DUKE_KF-145///dataset: DUKE///series: DUKE///id: KF-145///er: 1///pgr: 1///her2: NA///size: 5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 810///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
   samplename: DUKE_KF-146///dataset: DUKE///series: DUKE///id: KF-146///er: 0///pgr: 0///her2: NA///size: 2.2///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 535///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
 samplename: DUKE_KF-147///dataset: DUKE///series: DUKE///id: KF-147///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 2471///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
     samplename: DUKE_KF-149///dataset: DUKE///series: DUKE///id: KF-149///er: 0///pgr: 1///her2: NA///size: 3///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 775///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-150///dataset: DUKE///series: DUKE///id: KF-150///er: 1///pgr: 1///her2: NA///size: 3.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 2124///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-151///dataset: DUKE///series: DUKE///id: KF-151///er: 1///pgr: 1///her2: NA///size: 2.3///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1096///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-152///dataset: DUKE///series: DUKE///id: KF-152///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 842///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-153///dataset: DUKE///series: DUKE///id: KF-153///er: 0///pgr: 0///her2: NA///size: 7.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1248///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
    samplename: DUKE_KF-154///dataset: DUKE///series: DUKE///id: KF-154///er: 1///pgr: 1///her2: NA///size: 3///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 2251///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
   samplename: DUKE_KF-155///dataset: DUKE///series: DUKE///id: KF-155///er: 1///pgr: 1///her2: NA///size: 1.9///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 867///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-157///dataset: DUKE///series: DUKE///id: KF-157///er: 1///pgr: 1///her2: NA///size: 1.8///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1239///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
 samplename: DUKE_KF-159///dataset: DUKE///series: DUKE///id: KF-159///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1405///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
 samplename: DUKE_KF-160///dataset: DUKE///series: DUKE///id: KF-160///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1890///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
     samplename: DUKE_KF-176///dataset: DUKE///series: DUKE///id: KF-176///er: 1///pgr: 1///her2: NA///size: 3///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 603///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-178///dataset: DUKE///series: DUKE///id: KF-178///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 223///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
    samplename: DUKE_KF-182///dataset: DUKE///series: DUKE///id: KF-182///er: 1///pgr: 0///her2: NA///size: 3///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 4043///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-183///dataset: DUKE///series: DUKE///id: KF-183///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 755///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
    samplename: DUKE_KF-184///dataset: DUKE///series: DUKE///id: KF-184///er: 1///pgr: 1///her2: NA///size: 3///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 4069///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
   samplename: DUKE_KF-187///dataset: DUKE///series: DUKE///id: KF-187///er: 1///pgr: 0///her2: NA///size: 4.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 664///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
   samplename: DUKE_KF-188///dataset: DUKE///series: DUKE///id: KF-188///er: 1///pgr: 1///her2: NA///size: 1.8///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 747///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-189///dataset: DUKE///series: DUKE///id: KF-189///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 903///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-193///dataset: DUKE///series: DUKE///id: KF-193///er: 0///pgr: 0///her2: NA///size: 2.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 3500///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-194///dataset: DUKE///series: DUKE///id: KF-194///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 840///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
 samplename: DUKE_KF-195///dataset: DUKE///series: DUKE///id: KF-195///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1732///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-196///dataset: DUKE///series: DUKE///id: KF-196///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 951///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-200///dataset: DUKE///series: DUKE///id: KF-200///er: 0///pgr: 0///her2: NA///size: 1.3///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 2887///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-201///dataset: DUKE///series: DUKE///id: KF-201///er: 1///pgr: 1///her2: NA///size: 8.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1709///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
 samplename: DUKE_KF-203///dataset: DUKE///series: DUKE///id: KF-203///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 3288///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
 samplename: DUKE_KF-204///dataset: DUKE///series: DUKE///id: KF-204///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1751///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
  samplename: DUKE_KF-207///dataset: DUKE///series: DUKE///id: KF-207///er: 1///pgr: 1///her2: NA///size: 0.7///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 1598///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
samplename: DUKE_T00-050///dataset: DUKE///series: DUKE///id: T00-050///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 906///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
 samplename: DUKE_T00-064///dataset: DUKE///series: DUKE///id: T00-064///er: 1///pgr: 1///her2: NA///size: 2.8///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 884///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
 samplename: DUKE_T00-070///dataset: DUKE///series: DUKE///id: T00-070///er: 1///pgr: 1///her2: NA///size: 2.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 871///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
 samplename: DUKE_T00-109///dataset: DUKE///series: DUKE///id: T00-109///er: 1///pgr: 0///her2: NA///size: 5.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 842///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
samplename: DUKE_T00-214///dataset: DUKE///series: DUKE///id: T00-214///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 417///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
samplename: DUKE_T00-454///dataset: DUKE///series: DUKE///id: T00-454///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 306///e.os: 1///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
samplename: DUKE_T00-622///dataset: DUKE///series: DUKE///id: T00-622///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 508///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
 samplename: DUKE_T01-047///dataset: DUKE///series: DUKE///id: T01-047///er: 1///pgr: 1///her2: NA///size: 1.1///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 548///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
samplename: DUKE_T01-052///dataset: DUKE///series: DUKE///id: T01-052///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 456///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
samplename: DUKE_T01-145///dataset: DUKE///series: DUKE///id: T01-145///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 482///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
samplename: DUKE_T01-154///dataset: DUKE///series: DUKE///id: T01-154///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 368///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
samplename: DUKE_T01-156///dataset: DUKE///series: DUKE///id: T01-156///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 501///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
samplename: DUKE_T01-158///dataset: DUKE///series: DUKE///id: T01-158///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 494///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
samplename: DUKE_T01-166///dataset: DUKE///series: DUKE///id: T01-166///er: 0///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 462///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
   samplename: DUKE_T01-194///dataset: DUKE///series: DUKE///id: T01-194///er: 1///pgr: 1///her2: NA///size: 5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 520///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
samplename: DUKE_T01-202///dataset: DUKE///series: DUKE///id: T01-202///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 431///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
samplename: DUKE_T01-204///dataset: DUKE///series: DUKE///id: T01-204///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 487///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
 samplename: DUKE_T01-255///dataset: DUKE///series: DUKE///id: T01-255///er: 1///pgr: 1///her2: NA///size: 2.2///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 438///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
 samplename: DUKE_T01-310///dataset: DUKE///series: DUKE///id: T01-310///er: 1///pgr: 1///her2: NA///size: 1.8///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 417///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
samplename: DUKE_T01-344///dataset: DUKE///series: DUKE///id: T01-344///er: 1///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: 513///e.os: 0///tissue: 1///treatment: 99 
                                                                                                                                                                                                                                                 1 
                                                                                                                                                                                                                                           (Other) 
                                                                                                                                                                                                                                                72 

duplicates: 
   Length     Class      Mode 
      171 character character 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse3143 }
\keyword{datasets}
