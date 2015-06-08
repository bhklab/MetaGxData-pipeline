\name{ KOO_eset }
\alias{ KOO_eset }
\docType{data}
\title{ KOO_eset }
\description{ link does not work, in progress8 }
\usage{data( KOO_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/12747878 
  Title:  
  URL: Unavailable 
  PMIDs: 12747878 

  Abstract: A 6 word abstract is available. Use 'abstract' method.
  notes:
   summary:     
      A new gene signature was used to accurately predict 90% of relapse cases in the study.
   version:     
      2015-04-27 19:13:12
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: 1007_s_at 1053_at ... AFFX-HUMISGF3A/M97935_MB_at (280 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 280 features, 88 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
 KOO_KF_0076  KOO_KF_0077  KOO_KF_0079  KOO_KF_0110   KOO_KF_069   KOO_KF_094   KOO_KF_095   KOO_KF_097   KOO_KF_104 
           1            1            1            1            1            1            1            1            1 
  KOO_KF_105   KOO_KF_113   KOO_KF_117   KOO_KF_120   KOO_KF_125   KOO_KF_128   KOO_KF_129   KOO_KF_131   KOO_KF_133 
           1            1            1            1            1            1            1            1            1 
  KOO_KF_135   KOO_KF_138   KOO_KF_139   KOO_KF_140   KOO_KF_142   KOO_KF_144   KOO_KF_145   KOO_KF_146   KOO_KF_149 
           1            1            1            1            1            1            1            1            1 
  KOO_KF_150   KOO_KF_151   KOO_KF_153   KOO_KF_154   KOO_KF_155   KOO_KF_157   KOO_KF_176   KOO_KF_182   KOO_KF_184 
           1            1            1            1            1            1            1            1            1 
  KOO_KF_187   KOO_KF_188   KOO_KF_193   KOO_KF_200   KOO_KF_201   KOO_KF_207    KOO_KF_67 KOO_T00_0064 KOO_T00_0070 
           1            1            1            1            1            1            1            1            1 
KOO_T00_0109 KOO_T01_0310  KOO_T01_047 KOO_T01_0585  KOO_T01_194  KOO_T01_255  KOO_T01_407  KOO_T01_439  KOO_T01_461 
           1            1            1            1            1            1            1            1            1 
 KOO_T01_472  KOO_T01_502  KOO_T01_514  KOO_T01_618  KOO_T01_838  KOO_T92_037  KOO_T92_058  KOO_T93_002  KOO_T94_071 
           1            1            1            1            1            1            1            1            1 
 KOO_T94_072  KOO_T94_074  KOO_T95_028  KOO_T95_031  KOO_T95_040  KOO_T95_074  KOO_T95_162  KOO_T95_194  KOO_T96_085 
           1            1            1            1            1            1            1            1            1 
 KOO_T96_094  KOO_T96_136  KOO_T96_159  KOO_T96_167  KOO_T96_195  KOO_T96_250 KOO_T97_0075 KOO_T97_0080 KOO_T97_0090 
           1            1            1            1            1            1            1            1            1 
KOO_T97_0094  KOO_T97_017  KOO_T97_019  KOO_T97_043  KOO_T97_071  KOO_T97_191  KOO_T97_281 
           1            1            1            1            1            1            1 

alt_sample_name: 
 KF-0076  KF-0077  KF-0079  KF-0110   KF-069   KF-094   KF-095   KF-097   KF-104   KF-105   KF-113   KF-117   KF-120   KF-125 
       1        1        1        1        1        1        1        1        1        1        1        1        1        1 
  KF-128   KF-129   KF-131   KF-133   KF-135   KF-138   KF-139   KF-140   KF-142   KF-144   KF-145   KF-146   KF-149   KF-150 
       1        1        1        1        1        1        1        1        1        1        1        1        1        1 
  KF-151   KF-153   KF-154   KF-155   KF-157   KF-176   KF-182   KF-184   KF-187   KF-188   KF-193   KF-200   KF-201   KF-207 
       1        1        1        1        1        1        1        1        1        1        1        1        1        1 
   KF-67 T00-0064 T00-0070 T00-0109 T01-0310  T01-047 T01-0585  T01-194  T01-255  T01-407  T01-439  T01-461  T01-472  T01-502 
       1        1        1        1        1        1        1        1        1        1        1        1        1        1 
 T01-514  T01-618  T01-838  T92-037  T92-058  T93-002  T94-071  T94-072  T94-074  T95-028  T95-031  T95-040  T95-074  T95-162 
       1        1        1        1        1        1        1        1        1        1        1        1        1        1 
 T95-194  T96-085  T96-094  T96-136  T96-159  T96-167  T96-195  T96-250 T97-0075 T97-0080 T97-0090 T97-0094  T97-017  T97-019 
       1        1        1        1        1        1        1        1        1        1        1        1        1        1 
 T97-043  T97-071  T97-191  T97-281 
       1        1        1        1 

sample_type: 
tumor 
   88 

er: 
negative positive 
      15       73 

pgr: 
negative positive 
      23       65 

tumor_size: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   0.20    1.80    2.30    2.74    3.50    8.50 

N: 
 0  1 
19 69 

treatment: 
chemotherapy    untreated 
          61           27 

batch: 
KOO 
 88 

uncurated_author_metadata: 
  samplename: KOO_KF_0076///dataset: KOO///series: KOO///id: KF-0076///er: 1///pgr: 1///her2: NA///size: 1.7///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                              1 
    samplename: KOO_KF_0077///dataset: KOO///series: KOO///id: KF-0077///er: 1///pgr: 1///her2: NA///size: 2///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                              1 
    samplename: KOO_KF_0079///dataset: KOO///series: KOO///id: KF-0079///er: 1///pgr: 1///her2: NA///size: 2///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                              1 
    samplename: KOO_KF_0110///dataset: KOO///series: KOO///id: KF-0110///er: 1///pgr: 1///her2: NA///size: 2///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
      samplename: KOO_KF_069///dataset: KOO///series: KOO///id: KF-069///er: 1///pgr: 1///her2: NA///size: 4///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
    samplename: KOO_KF_094///dataset: KOO///series: KOO///id: KF-094///er: 1///pgr: 1///her2: NA///size: 1.8///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
    samplename: KOO_KF_095///dataset: KOO///series: KOO///id: KF-095///er: 1///pgr: 1///her2: NA///size: 0.6///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                              1 
    samplename: KOO_KF_097///dataset: KOO///series: KOO///id: KF-097///er: 1///pgr: 1///her2: NA///size: 2.2///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
    samplename: KOO_KF_104///dataset: KOO///series: KOO///id: KF-104///er: 1///pgr: 1///her2: NA///size: 1.2///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                              1 
    samplename: KOO_KF_105///dataset: KOO///series: KOO///id: KF-105///er: 1///pgr: 1///her2: NA///size: 1.3///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                              1 
    samplename: KOO_KF_113///dataset: KOO///series: KOO///id: KF-113///er: 1///pgr: 1///her2: NA///size: 1.5///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                              1 
    samplename: KOO_KF_117///dataset: KOO///series: KOO///id: KF-117///er: 1///pgr: 1///her2: NA///size: 1.4///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                              1 
    samplename: KOO_KF_120///dataset: KOO///series: KOO///id: KF-120///er: 1///pgr: 1///her2: NA///size: 1.5///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                              1 
      samplename: KOO_KF_125///dataset: KOO///series: KOO///id: KF-125///er: 1///pgr: 1///her2: NA///size: 4///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
    samplename: KOO_KF_128///dataset: KOO///series: KOO///id: KF-128///er: 1///pgr: 1///her2: NA///size: 2.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                              1 
      samplename: KOO_KF_129///dataset: KOO///series: KOO///id: KF-129///er: 1///pgr: 0///her2: NA///size: 2///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
      samplename: KOO_KF_131///dataset: KOO///series: KOO///id: KF-131///er: 1///pgr: 1///her2: NA///size: 3///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
    samplename: KOO_KF_133///dataset: KOO///series: KOO///id: KF-133///er: 1///pgr: 1///her2: NA///size: 0.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
    samplename: KOO_KF_135///dataset: KOO///series: KOO///id: KF-135///er: 0///pgr: 0///her2: NA///size: 2.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
      samplename: KOO_KF_138///dataset: KOO///series: KOO///id: KF-138///er: 1///pgr: 1///her2: NA///size: 4///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
    samplename: KOO_KF_139///dataset: KOO///series: KOO///id: KF-139///er: 1///pgr: 1///her2: NA///size: 3.1///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
    samplename: KOO_KF_140///dataset: KOO///series: KOO///id: KF-140///er: 0///pgr: 0///her2: NA///size: 3.7///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
    samplename: KOO_KF_142///dataset: KOO///series: KOO///id: KF-142///er: 1///pgr: 0///her2: NA///size: 3.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                              1 
    samplename: KOO_KF_144///dataset: KOO///series: KOO///id: KF-144///er: 1///pgr: 1///her2: NA///size: 0.5///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                              1 
      samplename: KOO_KF_145///dataset: KOO///series: KOO///id: KF-145///er: 1///pgr: 1///her2: NA///size: 5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
    samplename: KOO_KF_146///dataset: KOO///series: KOO///id: KF-146///er: 0///pgr: 0///her2: NA///size: 2.2///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
      samplename: KOO_KF_149///dataset: KOO///series: KOO///id: KF-149///er: 0///pgr: 1///her2: NA///size: 3///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
    samplename: KOO_KF_150///dataset: KOO///series: KOO///id: KF-150///er: 1///pgr: 1///her2: NA///size: 3.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
    samplename: KOO_KF_151///dataset: KOO///series: KOO///id: KF-151///er: 1///pgr: 1///her2: NA///size: 2.3///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
    samplename: KOO_KF_153///dataset: KOO///series: KOO///id: KF-153///er: 0///pgr: 0///her2: NA///size: 7.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
      samplename: KOO_KF_154///dataset: KOO///series: KOO///id: KF-154///er: 1///pgr: 1///her2: NA///size: 3///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
    samplename: KOO_KF_155///dataset: KOO///series: KOO///id: KF-155///er: 1///pgr: 1///her2: NA///size: 1.9///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
    samplename: KOO_KF_157///dataset: KOO///series: KOO///id: KF-157///er: 1///pgr: 1///her2: NA///size: 1.8///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
      samplename: KOO_KF_176///dataset: KOO///series: KOO///id: KF-176///er: 1///pgr: 1///her2: NA///size: 3///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                              1 
      samplename: KOO_KF_182///dataset: KOO///series: KOO///id: KF-182///er: 1///pgr: 0///her2: NA///size: 3///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                              1 
      samplename: KOO_KF_184///dataset: KOO///series: KOO///id: KF-184///er: 1///pgr: 1///her2: NA///size: 3///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
    samplename: KOO_KF_187///dataset: KOO///series: KOO///id: KF-187///er: 1///pgr: 0///her2: NA///size: 4.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                              1 
    samplename: KOO_KF_188///dataset: KOO///series: KOO///id: KF-188///er: 1///pgr: 1///her2: NA///size: 1.8///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
    samplename: KOO_KF_193///dataset: KOO///series: KOO///id: KF-193///er: 0///pgr: 0///her2: NA///size: 2.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
    samplename: KOO_KF_200///dataset: KOO///series: KOO///id: KF-200///er: 0///pgr: 0///her2: NA///size: 1.3///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                              1 
    samplename: KOO_KF_201///dataset: KOO///series: KOO///id: KF-201///er: 1///pgr: 1///her2: NA///size: 8.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
    samplename: KOO_KF_207///dataset: KOO///series: KOO///id: KF-207///er: 1///pgr: 1///her2: NA///size: 0.7///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                              1 
      samplename: KOO_KF_67///dataset: KOO///series: KOO///id: KF-67///er: 1///pgr: 1///her2: NA///size: 7.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                              1 
samplename: KOO_T00_0064///dataset: KOO///series: KOO///id: T00-0064///er: 1///pgr: 1///her2: NA///size: 2.8///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
samplename: KOO_T00_0070///dataset: KOO///series: KOO///id: T00-0070///er: 1///pgr: 1///her2: NA///size: 2.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
samplename: KOO_T00_0109///dataset: KOO///series: KOO///id: T00-0109///er: 1///pgr: 0///her2: NA///size: 5.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
samplename: KOO_T01_0310///dataset: KOO///series: KOO///id: T01-0310///er: 1///pgr: 1///her2: NA///size: 1.8///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
  samplename: KOO_T01_047///dataset: KOO///series: KOO///id: T01-047///er: 1///pgr: 1///her2: NA///size: 1.1///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
samplename: KOO_T01_0585///dataset: KOO///series: KOO///id: T01-0585///er: 1///pgr: 1///her2: NA///size: 3.8///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                              1 
    samplename: KOO_T01_194///dataset: KOO///series: KOO///id: T01-194///er: 1///pgr: 1///her2: NA///size: 5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
  samplename: KOO_T01_255///dataset: KOO///series: KOO///id: T01-255///er: 1///pgr: 1///her2: NA///size: 2.2///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
  samplename: KOO_T01_407///dataset: KOO///series: KOO///id: T01-407///er: 1///pgr: 1///her2: NA///size: 1.8///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
  samplename: KOO_T01_439///dataset: KOO///series: KOO///id: T01-439///er: 1///pgr: 1///her2: NA///size: 1.8///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                              1 
  samplename: KOO_T01_461///dataset: KOO///series: KOO///id: T01-461///er: 1///pgr: 1///her2: NA///size: 1.8///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
  samplename: KOO_T01_472///dataset: KOO///series: KOO///id: T01-472///er: 1///pgr: 1///her2: NA///size: 4.2///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
  samplename: KOO_T01_502///dataset: KOO///series: KOO///id: T01-502///er: 1///pgr: 1///her2: NA///size: 6.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
    samplename: KOO_T01_514///dataset: KOO///series: KOO///id: T01-514///er: 1///pgr: 1///her2: NA///size: 2///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
    samplename: KOO_T01_618///dataset: KOO///series: KOO///id: T01-618///er: 1///pgr: 1///her2: NA///size: 2///node: 0///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
  samplename: KOO_T01_838///dataset: KOO///series: KOO///id: T01-838///er: 1///pgr: 1///her2: NA///size: 3.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
    samplename: KOO_T92_037///dataset: KOO///series: KOO///id: T92-037///er: 1///pgr: 0///her2: NA///size: 2///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
  samplename: KOO_T92_058///dataset: KOO///series: KOO///id: T92-058///er: 1///pgr: 0///her2: NA///size: 2.2///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
  samplename: KOO_T93_002///dataset: KOO///series: KOO///id: T93-002///er: 1///pgr: 0///her2: NA///size: 3.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
    samplename: KOO_T94_071///dataset: KOO///series: KOO///id: T94-071///er: 1///pgr: 1///her2: NA///size: 4///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
    samplename: KOO_T94_072///dataset: KOO///series: KOO///id: T94-072///er: 1///pgr: 0///her2: NA///size: 3///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                              1 
    samplename: KOO_T94_074///dataset: KOO///series: KOO///id: T94-074///er: 1///pgr: 1///her2: NA///size: 3///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
  samplename: KOO_T95_028///dataset: KOO///series: KOO///id: T95-028///er: 0///pgr: 0///her2: NA///size: 1.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                              1 
    samplename: KOO_T95_031///dataset: KOO///series: KOO///id: T95-031///er: 1///pgr: 1///her2: NA///size: 2///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                              1 
  samplename: KOO_T95_040///dataset: KOO///series: KOO///id: T95-040///er: 1///pgr: 1///her2: NA///size: 2.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
  samplename: KOO_T95_074///dataset: KOO///series: KOO///id: T95-074///er: 0///pgr: 0///her2: NA///size: 4.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
  samplename: KOO_T95_162///dataset: KOO///series: KOO///id: T95-162///er: 1///pgr: 0///her2: NA///size: 2.3///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                              1 
  samplename: KOO_T95_194///dataset: KOO///series: KOO///id: T95-194///er: 1///pgr: 1///her2: NA///size: 0.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
  samplename: KOO_T96_085///dataset: KOO///series: KOO///id: T96-085///er: 0///pgr: 1///her2: NA///size: 2.2///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
  samplename: KOO_T96_094///dataset: KOO///series: KOO///id: T96-094///er: 1///pgr: 1///her2: NA///size: 1.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                              1 
    samplename: KOO_T96_136///dataset: KOO///series: KOO///id: T96-136///er: 0///pgr: 0///her2: NA///size: 3///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
  samplename: KOO_T96_159///dataset: KOO///series: KOO///id: T96-159///er: 1///pgr: 1///her2: NA///size: 3.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
  samplename: KOO_T96_167///dataset: KOO///series: KOO///id: T96-167///er: 1///pgr: 1///her2: NA///size: 1.2///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
    samplename: KOO_T96_195///dataset: KOO///series: KOO///id: T96-195///er: 1///pgr: 1///her2: NA///size: 3///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
  samplename: KOO_T96_250///dataset: KOO///series: KOO///id: T96-250///er: 0///pgr: 0///her2: NA///size: 2.2///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
samplename: KOO_T97_0075///dataset: KOO///series: KOO///id: T97-0075///er: 1///pgr: 1///her2: NA///size: 1.7///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
samplename: KOO_T97_0080///dataset: KOO///series: KOO///id: T97-0080///er: 1///pgr: 1///her2: NA///size: 2.2///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                              1 
samplename: KOO_T97_0090///dataset: KOO///series: KOO///id: T97-0090///er: 1///pgr: 1///her2: NA///size: 2.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
samplename: KOO_T97_0094///dataset: KOO///series: KOO///id: T97-0094///er: 0///pgr: 0///her2: NA///size: 2.8///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
  samplename: KOO_T97_017///dataset: KOO///series: KOO///id: T97-017///er: 0///pgr: 0///her2: NA///size: 1.8///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
    samplename: KOO_T97_019///dataset: KOO///series: KOO///id: T97-019///er: 1///pgr: 1///her2: NA///size: 2///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
    samplename: KOO_T97_043///dataset: KOO///series: KOO///id: T97-043///er: 1///pgr: 1///her2: NA///size: 4///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                              1 
  samplename: KOO_T97_071///dataset: KOO///series: KOO///id: T97-071///er: 1///pgr: 1///her2: NA///size: 4.5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
  samplename: KOO_T97_191///dataset: KOO///series: KOO///id: T97-191///er: 1///pgr: 1///her2: NA///size: 0.2///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 
    samplename: KOO_T97_281///dataset: KOO///series: KOO///id: T97-281///er: 0///pgr: 0///her2: NA///size: 5///node: 1///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                              1 

duplicates: 
   Length     Class      Mode 
       88 character character 

}}
\source{ Unavailable }
\keyword{datasets}
