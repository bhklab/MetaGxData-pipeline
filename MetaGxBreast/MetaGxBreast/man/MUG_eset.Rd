\name{ MUG_eset }
\alias{ MUG_eset }
\docType{data}
\title{  }
\usage{data( MUG_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=18636107 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse5364 
  PMIDs: 18636107 
  No abstract available.
  notes:
   summary:     
      48 genes were identified that displayed highly restricted levels of expression in tumors compared to normal tissues. This
 was validated in 11 independent cohorts of different cancer types.
   version:     
      2015-04-27 19:21:48
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: H200000001 H200000005 ... opHsV04TC000043 (14288 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 14288 features, 152 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
MUG_IDC_A001 MUG_IDC_A002 MUG_IDC_A003 MUG_IDC_A004 MUG_IDC_A005 MUG_IDC_A006 MUG_IDC_A007 MUG_IDC_A014 MUG_IDC_A015 
           1            1            1            1            1            1            1            1            1 
MUG_IDC_A016 MUG_IDC_A018 MUG_IDC_A019 MUG_IDC_A020 MUG_IDC_A021 MUG_IDC_A022 MUG_IDC_A024 MUG_IDC_A025 MUG_IDC_A026 
           1            1            1            1            1            1            1            1            1 
MUG_IDC_A027 MUG_IDC_A032 MUG_IDC_A033 MUG_IDC_A035 MUG_IDC_A036 MUG_IDC_A037 MUG_IDC_A038 MUG_IDC_A039 MUG_IDC_A040 
           1            1            1            1            1            1            1            1            1 
MUG_IDC_A044 MUG_IDC_A046 MUG_IDC_A047 MUG_IDC_A049 MUG_IDC_A051 MUG_IDC_A052 MUG_IDC_A053 MUG_IDC_A054 MUG_IDC_A056 
           1            1            1            1            1            1            1            1            1 
MUG_IDC_A057 MUG_IDC_A058 MUG_IDC_A059 MUG_IDC_A062 MUG_IDC_A063 MUG_IDC_A068 MUG_IDC_A069 MUG_IDC_A070 MUG_IDC_A074 
           1            1            1            1            1            1            1            1            1 
MUG_IDC_A075 MUG_IDC_A079 MUG_IDC_A080 MUG_IDC_A083 MUG_IDC_A084 MUG_IDC_A085 MUG_IDC_A087 MUG_IDC_A093 MUG_IDC_A096 
           1            1            1            1            1            1            1            1            1 
MUG_IDC_A097 MUG_IDC_A099 MUG_IDC_A102 MUG_IDC_A103 MUG_IDC_A104 MUG_IDC_A105 MUG_IDC_A107 MUG_IDC_A108 MUG_IDC_B013 
           1            1            1            1            1            1            1            1            1 
MUG_IDC_B014 MUG_IDC_B015 MUG_IDC_B017 MUG_IDC_B022 MUG_IDC_B023 MUG_IDC_B024 MUG_IDC_B027 MUG_IDC_B028 MUG_IDC_B031 
           1            1            1            1            1            1            1            1            1 
MUG_IDC_B033 MUG_IDC_B034 MUG_IDC_B036 MUG_IDC_B037 MUG_IDC_B038 MUG_IDC_B040 MUG_IDC_B041 MUG_IDC_B043 MUG_IDC_B044 
           1            1            1            1            1            1            1            1            1 
MUG_IDC_B045 MUG_IDC_B046 MUG_IDC_B047 MUG_IDC_B048 MUG_IDC_B049 MUG_IDC_B050 MUG_IDC_B051 MUG_IDC_B052 MUG_IDC_B054 
           1            1            1            1            1            1            1            1            1 
MUG_IDC_B055 MUG_IDC_B059 MUG_IDC_B060 MUG_IDC_B061 MUG_IDC_B062 MUG_IDC_B063 MUG_IDC_B065 MUG_IDC_B069 MUG_IDC_B071 
           1            1            1            1            1            1            1            1            1 
     (Other) 
          53 

alt_sample_name: 
_IDC_A001 _IDC_A002 _IDC_A003 _IDC_A004 _IDC_A005 _IDC_A006 _IDC_A007 _IDC_A014 _IDC_A015 _IDC_A016 _IDC_A018 _IDC_A019 
        1         1         1         1         1         1         1         1         1         1         1         1 
_IDC_A020 _IDC_A021 _IDC_A022 _IDC_A024 _IDC_A025 _IDC_A026 _IDC_A027 _IDC_A032 _IDC_A033 _IDC_A035 _IDC_A036 _IDC_A037 
        1         1         1         1         1         1         1         1         1         1         1         1 
_IDC_A038 _IDC_A039 _IDC_A040 _IDC_A044 _IDC_A046 _IDC_A047 _IDC_A049 _IDC_A051 _IDC_A052 _IDC_A053 _IDC_A054 _IDC_A056 
        1         1         1         1         1         1         1         1         1         1         1         1 
_IDC_A057 _IDC_A058 _IDC_A059 _IDC_A062 _IDC_A063 _IDC_A068 _IDC_A069 _IDC_A070 _IDC_A074 _IDC_A075 _IDC_A079 _IDC_A080 
        1         1         1         1         1         1         1         1         1         1         1         1 
_IDC_A083 _IDC_A084 _IDC_A085 _IDC_A087 _IDC_A093 _IDC_A096 _IDC_A097 _IDC_A099 _IDC_A102 _IDC_A103 _IDC_A104 _IDC_A105 
        1         1         1         1         1         1         1         1         1         1         1         1 
_IDC_A107 _IDC_A108 _IDC_B013 _IDC_B014 _IDC_B015 _IDC_B017 _IDC_B022 _IDC_B023 _IDC_B024 _IDC_B027 _IDC_B028 _IDC_B031 
        1         1         1         1         1         1         1         1         1         1         1         1 
_IDC_B033 _IDC_B034 _IDC_B036 _IDC_B037 _IDC_B038 _IDC_B040 _IDC_B041 _IDC_B043 _IDC_B044 _IDC_B045 _IDC_B046 _IDC_B047 
        1         1         1         1         1         1         1         1         1         1         1         1 
_IDC_B048 _IDC_B049 _IDC_B050 _IDC_B051 _IDC_B052 _IDC_B054 _IDC_B055 _IDC_B059 _IDC_B060 _IDC_B061 _IDC_B062 _IDC_B063 
        1         1         1         1         1         1         1         1         1         1         1         1 
_IDC_B065 _IDC_B069 _IDC_B071   (Other) 
        1         1         1        53 

sample_type: 
tumor 
  152 

batch: 
MUG 
152 

uncurated_author_metadata: 
samplename: MUG_IDC_A001///dataset: MUG///series: MUG///id: _IDC_A001///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A002///dataset: MUG///series: MUG///id: _IDC_A002///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A003///dataset: MUG///series: MUG///id: _IDC_A003///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A004///dataset: MUG///series: MUG///id: _IDC_A004///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A005///dataset: MUG///series: MUG///id: _IDC_A005///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A006///dataset: MUG///series: MUG///id: _IDC_A006///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A007///dataset: MUG///series: MUG///id: _IDC_A007///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A014///dataset: MUG///series: MUG///id: _IDC_A014///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A015///dataset: MUG///series: MUG///id: _IDC_A015///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A016///dataset: MUG///series: MUG///id: _IDC_A016///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A018///dataset: MUG///series: MUG///id: _IDC_A018///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A019///dataset: MUG///series: MUG///id: _IDC_A019///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A020///dataset: MUG///series: MUG///id: _IDC_A020///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A021///dataset: MUG///series: MUG///id: _IDC_A021///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A022///dataset: MUG///series: MUG///id: _IDC_A022///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A024///dataset: MUG///series: MUG///id: _IDC_A024///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A025///dataset: MUG///series: MUG///id: _IDC_A025///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A026///dataset: MUG///series: MUG///id: _IDC_A026///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A027///dataset: MUG///series: MUG///id: _IDC_A027///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A032///dataset: MUG///series: MUG///id: _IDC_A032///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A033///dataset: MUG///series: MUG///id: _IDC_A033///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A035///dataset: MUG///series: MUG///id: _IDC_A035///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A036///dataset: MUG///series: MUG///id: _IDC_A036///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A037///dataset: MUG///series: MUG///id: _IDC_A037///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A038///dataset: MUG///series: MUG///id: _IDC_A038///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A039///dataset: MUG///series: MUG///id: _IDC_A039///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A040///dataset: MUG///series: MUG///id: _IDC_A040///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A044///dataset: MUG///series: MUG///id: _IDC_A044///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A046///dataset: MUG///series: MUG///id: _IDC_A046///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A047///dataset: MUG///series: MUG///id: _IDC_A047///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A049///dataset: MUG///series: MUG///id: _IDC_A049///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A051///dataset: MUG///series: MUG///id: _IDC_A051///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A052///dataset: MUG///series: MUG///id: _IDC_A052///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A053///dataset: MUG///series: MUG///id: _IDC_A053///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A054///dataset: MUG///series: MUG///id: _IDC_A054///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A056///dataset: MUG///series: MUG///id: _IDC_A056///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A057///dataset: MUG///series: MUG///id: _IDC_A057///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A058///dataset: MUG///series: MUG///id: _IDC_A058///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A059///dataset: MUG///series: MUG///id: _IDC_A059///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A062///dataset: MUG///series: MUG///id: _IDC_A062///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A063///dataset: MUG///series: MUG///id: _IDC_A063///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A068///dataset: MUG///series: MUG///id: _IDC_A068///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A069///dataset: MUG///series: MUG///id: _IDC_A069///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A070///dataset: MUG///series: MUG///id: _IDC_A070///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A074///dataset: MUG///series: MUG///id: _IDC_A074///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A075///dataset: MUG///series: MUG///id: _IDC_A075///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A079///dataset: MUG///series: MUG///id: _IDC_A079///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A080///dataset: MUG///series: MUG///id: _IDC_A080///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A083///dataset: MUG///series: MUG///id: _IDC_A083///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A084///dataset: MUG///series: MUG///id: _IDC_A084///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A085///dataset: MUG///series: MUG///id: _IDC_A085///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A087///dataset: MUG///series: MUG///id: _IDC_A087///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A093///dataset: MUG///series: MUG///id: _IDC_A093///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A096///dataset: MUG///series: MUG///id: _IDC_A096///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A097///dataset: MUG///series: MUG///id: _IDC_A097///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A099///dataset: MUG///series: MUG///id: _IDC_A099///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A102///dataset: MUG///series: MUG///id: _IDC_A102///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A103///dataset: MUG///series: MUG///id: _IDC_A103///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A104///dataset: MUG///series: MUG///id: _IDC_A104///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A105///dataset: MUG///series: MUG///id: _IDC_A105///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A107///dataset: MUG///series: MUG///id: _IDC_A107///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_A108///dataset: MUG///series: MUG///id: _IDC_A108///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B013///dataset: MUG///series: MUG///id: _IDC_B013///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B014///dataset: MUG///series: MUG///id: _IDC_B014///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B015///dataset: MUG///series: MUG///id: _IDC_B015///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B017///dataset: MUG///series: MUG///id: _IDC_B017///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B022///dataset: MUG///series: MUG///id: _IDC_B022///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B023///dataset: MUG///series: MUG///id: _IDC_B023///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B024///dataset: MUG///series: MUG///id: _IDC_B024///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B027///dataset: MUG///series: MUG///id: _IDC_B027///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B028///dataset: MUG///series: MUG///id: _IDC_B028///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B031///dataset: MUG///series: MUG///id: _IDC_B031///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B033///dataset: MUG///series: MUG///id: _IDC_B033///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B034///dataset: MUG///series: MUG///id: _IDC_B034///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B036///dataset: MUG///series: MUG///id: _IDC_B036///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B037///dataset: MUG///series: MUG///id: _IDC_B037///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B038///dataset: MUG///series: MUG///id: _IDC_B038///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B040///dataset: MUG///series: MUG///id: _IDC_B040///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B041///dataset: MUG///series: MUG///id: _IDC_B041///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B043///dataset: MUG///series: MUG///id: _IDC_B043///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B044///dataset: MUG///series: MUG///id: _IDC_B044///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B045///dataset: MUG///series: MUG///id: _IDC_B045///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B046///dataset: MUG///series: MUG///id: _IDC_B046///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B047///dataset: MUG///series: MUG///id: _IDC_B047///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B048///dataset: MUG///series: MUG///id: _IDC_B048///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B049///dataset: MUG///series: MUG///id: _IDC_B049///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B050///dataset: MUG///series: MUG///id: _IDC_B050///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B051///dataset: MUG///series: MUG///id: _IDC_B051///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B052///dataset: MUG///series: MUG///id: _IDC_B052///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B054///dataset: MUG///series: MUG///id: _IDC_B054///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B055///dataset: MUG///series: MUG///id: _IDC_B055///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B059///dataset: MUG///series: MUG///id: _IDC_B059///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B060///dataset: MUG///series: MUG///id: _IDC_B060///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B061///dataset: MUG///series: MUG///id: _IDC_B061///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B062///dataset: MUG///series: MUG///id: _IDC_B062///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B063///dataset: MUG///series: MUG///id: _IDC_B063///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B065///dataset: MUG///series: MUG///id: _IDC_B065///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B069///dataset: MUG///series: MUG///id: _IDC_B069///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
samplename: MUG_IDC_B071///dataset: MUG///series: MUG///id: _IDC_B071///er: NA///pgr: NA///her2: NA///size: NA///node: NA///age: NA///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                                  1 
                                                                                                                                                                                                                                            (Other) 
                                                                                                                                                                                                                                                 53 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse5364 }
\keyword{datasets}
