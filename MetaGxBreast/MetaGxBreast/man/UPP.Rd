\name{ UPP_eset }
\alias{ UPP_eset }
\docType{data}
\title{ UPP_eset }
\usage{data( UPP_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: NA 
  Title:  
  URL: NA 
  PMIDs: NA 

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: 1007_s_at 1053_at ... 244889_at (36178 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 36178 features, 251 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
UPP_100B08 UPP_101B88 UPP_102B06 UPP_103B41 UPP_104B91 UPP_105B13 UPP_106B55  UPP_10B88 UPP_110B34 UPP_111B51 UPP_112B55 
         1          1          1          1          1          1          1          1          1          1          1 
UPP_113B11 UPP_114B68  UPP_11B47 UPP_120B73 UPP_122B81 UPP_124B25 UPP_127B00 UPP_128B48 UPP_130B92 UPP_131B79 UPP_134B33 
         1          1          1          1          1          1          1          1          1          1          1 
UPP_135B40 UPP_136B04 UPP_137B88 UPP_138B34 UPP_139B03  UPP_13B79 UPP_140B91 UPP_142B05 UPP_143B81 UPP_144B49 UPP_145B10 
         1          1          1          1          1          1          1          1          1          1          1 
UPP_146B39 UPP_147B19 UPP_148B49  UPP_14B98 UPP_150B81 UPP_151B84 UPP_152B99 UPP_153B09 UPP_154B42 UPP_155B52 UPP_156B01 
         1          1          1          1          1          1          1          1          1          1          1 
UPP_158B84 UPP_159B47  UPP_15C94 UPP_160B16 UPP_161B31 UPP_162B98 UPP_163B27 UPP_164B81 UPP_165B72 UPP_166B79 UPP_168B51 
         1          1          1          1          1          1          1          1          1          1          1 
UPP_169B79  UPP_16C97 UPP_170B15 UPP_171B77 UPP_172B19 UPP_173B43 UPP_174B41 UPP_175B72 UPP_176B74 UPP_177B67 UPP_178B74 
         1          1          1          1          1          1          1          1          1          1          1 
UPP_179B28  UPP_17C40 UPP_180B38 UPP_181B70 UPP_182B43 UPP_183B75 UPP_184B38 UPP_185B44 UPP_186B22 UPP_187B36 UPP_188B13 
         1          1          1          1          1          1          1          1          1          1          1 
UPP_189B83  UPP_18C56 UPP_191B79 UPP_192B69 UPP_193B72 UPP_194B60 UPP_195B75 UPP_196B81 UPP_197B95 UPP_198B90 UPP_199B55 
         1          1          1          1          1          1          1          1          1          1          1 
 UPP_19C33 UPP_200B47 UPP_201B68 UPP_202B44 UPP_203B49 UPP_204B85 UPP_205B99 UPP_206C05 UPP_207C08 UPP_208C06 UPP_209C10 
         1          1          1          1          1          1          1          1          1          1          1 
   (Other) 
       152 

alt_sample_name: 
 100B08  101B88  102B06  103B41  104B91  105B13  106B55   10B88  110B34  111B51  112B55  113B11  114B68   11B47  120B73 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
 122B81  124B25  127B00  128B48  130B92  131B79  134B33  135B40  136B04  137B88  138B34  139B03   13B79  140B91  142B05 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
 143B81  144B49  145B10  146B39  147B19  148B49   14B98  150B81  151B84  152B99  153B09  154B42  155B52  156B01  158B84 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
 159B47   15C94  160B16  161B31  162B98  163B27  164B81  165B72  166B79  168B51  169B79   16C97  170B15  171B77  172B19 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
 173B43  174B41  175B72  176B74  177B67  178B74  179B28   17C40  180B38  181B70  182B43  183B75  184B38  185B44  186B22 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
 187B36  188B13  189B83   18C56  191B79  192B69  193B72  194B60  195B75  196B81  197B95  198B90  199B55   19C33  200B47 
      1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 
 201B68  202B44  203B49  204B85  205B99  206C05  207C08  208C06  209C10 (Other) 
      1       1       1       1       1       1       1       1       1     152 

sample_type: 
tumor 
  251 

er: 
negative positive     NA's 
      34      213        4 

pgr: 
negative positive 
      61      190 

tumor_size: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  0.200   1.500   2.000   2.243   2.562  13.000 

N: 
   0    1 NA's 
 158   84    9 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  28.00   52.00   64.00   62.11   72.00   93.00 

grade: 
   1    2    3 NA's 
  67  128   54    2 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
     30    1870    3711    3007    3985    4654      17 

recurrence_status: 
deceased_or_recurrence    living_norecurrence                   NA's 
                    55                    181                     15 

treatment: 
hormonotherapy      untreated           NA's 
            80            142             29 

batch: 
UPPT UPPU 
  80  171 

uncurated_author_metadata: 
  samplename: UPP_100B08///dataset: UPP///series: UPPU///id: 100B08///er: 1///pgr: 0///her2: NA///size: 0.9///node: 0///age: 68///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 4319///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_101B88///dataset: UPP///series: UPPU///id: 101B88///er: 0///pgr: 0///her2: NA///size: 1.2///node: 0///age: 40///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 4319///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_102B06///dataset: UPP///series: UPPU///id: 102B06///er: 1///pgr: 1///her2: NA///size: 2.6///node: 0///age: 51///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 4319///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
   samplename: UPP_103B41///dataset: UPP///series: UPPT///id: 103B41///er: 0///pgr: 0///her2: NA///size: 2.2///node: 1///age: 52///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                             1 
 samplename: UPP_104B91///dataset: UPP///series: UPPT///id: 104B91///er: 1///pgr: 1///her2: NA///size: 2.4///node: NA///age: 80///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1308///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                             1 
  samplename: UPP_105B13///dataset: UPP///series: UPPU///id: 105B13///er: 1///pgr: 1///her2: NA///size: 1.3///node: 0///age: 46///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 4289///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
    samplename: UPP_106B55///dataset: UPP///series: UPPU///id: 106B55///er: 1///pgr: 1///her2: NA///size: 6///node: 0///age: 37///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2525///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
      samplename: UPP_10B88///dataset: UPP///series: UPPU///id: 10B88///er: 1///pgr: 1///her2: NA///size: 5///node: 0///age: 70///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 4137///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
    samplename: UPP_110B34///dataset: UPP///series: UPPU///id: 110B34///er: 1///pgr: 1///her2: NA///size: 2///node: 0///age: 74///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 4258///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_111B51///dataset: UPP///series: UPPU///id: 111B51///er: 1///pgr: 1///her2: NA///size: 3.3///node: 0///age: 41///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2616///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
   samplename: UPP_112B55///dataset: UPP///series: UPPT///id: 112B55///er: 1///pgr: 1///her2: NA///size: 2.4///node: 1///age: 61///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 548///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                             1 
  samplename: UPP_113B11///dataset: UPP///series: UPPU///id: 113B11///er: 1///pgr: 1///her2: NA///size: 3.2///node: 0///age: 38///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2008///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_114B68///dataset: UPP///series: UPPT///id: 114B68///er: 1///pgr: 1///her2: NA///size: 1.2///node: 1///age: 67///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 4076///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                             1 
    samplename: UPP_11B47///dataset: UPP///series: UPPU///id: 11B47///er: 0///pgr: 0///her2: NA///size: 1.8///node: 0///age: 69///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2707///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_120B73///dataset: UPP///series: UPPU///id: 120B73///er: 1///pgr: 1///her2: NA///size: 1.6///node: 0///age: 34///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 4228///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_122B81///dataset: UPP///series: UPPU///id: 122B81///er: 1///pgr: 1///her2: NA///size: 1.9///node: 0///age: 79///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 4076///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_124B25///dataset: UPP///series: UPPU///id: 124B25///er: 1///pgr: 1///her2: NA///size: 2.1///node: 0///age: 46///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 3285///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_127B00///dataset: UPP///series: UPPU///id: 127B00///er: 1///pgr: 1///her2: NA///size: 2.2///node: 0///age: 57///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1703///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
samplename: UPP_128B48///dataset: UPP///series: UPPU///id: 128B48///er: NA///pgr: 0///her2: NA///size: 1.2///node: NA///age: 71///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 1673///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_130B92///dataset: UPP///series: UPPT///id: 130B92///er: 1///pgr: 1///her2: NA///size: 1.8///node: 1///age: 73///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2038///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                             1 
  samplename: UPP_131B79///dataset: UPP///series: UPPT///id: 131B79///er: 1///pgr: 1///her2: NA///size: 2.6///node: 1///age: 59///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2129///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                             1 
  samplename: UPP_134B33///dataset: UPP///series: UPPU///id: 134B33///er: 1///pgr: 1///her2: NA///size: 2.8///node: 0///age: 63///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 3011///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_135B40///dataset: UPP///series: UPPT///id: 135B40///er: 1///pgr: 1///her2: NA///size: 2.6///node: 0///age: 64///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 4198///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                             1 
  samplename: UPP_136B04///dataset: UPP///series: UPPU///id: 136B04///er: 1///pgr: 1///her2: NA///size: 1.7///node: 0///age: 54///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 4198///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_137B88///dataset: UPP///series: UPPU///id: 137B88///er: 0///pgr: 0///her2: NA///size: 0.6///node: 0///age: 68///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 4198///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_138B34///dataset: UPP///series: UPPT///id: 138B34///er: 1///pgr: 1///her2: NA///size: 2.3///node: 1///age: 65///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 4198///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                             1 
   samplename: UPP_139B03///dataset: UPP///series: UPPT///id: 139B03///er: 1///pgr: 0///her2: NA///size: 2.8///node: 0///age: 84///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 122///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                             1 
    samplename: UPP_13B79///dataset: UPP///series: UPPU///id: 13B79///er: 1///pgr: 1///her2: NA///size: 2.5///node: 0///age: 67///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 3954///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
   samplename: UPP_140B91///dataset: UPP///series: UPPU///id: 140B91///er: 1///pgr: 1///her2: NA///size: 1.2///node: 0///age: 61///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_142B05///dataset: UPP///series: UPPU///id: 142B05///er: 0///pgr: 0///her2: NA///size: 2.7///node: 0///age: 67///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 4198///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_143B81///dataset: UPP///series: UPPU///id: 143B81///er: 0///pgr: 0///her2: NA///size: 2.9///node: NA///age: 93///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 791///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_144B49///dataset: UPP///series: UPPU///id: 144B49///er: 1///pgr: 1///her2: NA///size: 2.1///node: 0///age: 40///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 4198///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_145B10///dataset: UPP///series: UPPU///id: 145B10///er: 1///pgr: 1///her2: NA///size: 1.8///node: 0///age: 73///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 4167///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
    samplename: UPP_146B39///dataset: UPP///series: UPPU///id: 146B39///er: 0///pgr: 0///her2: NA///size: 5///node: 1///age: 34///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
  samplename: UPP_147B19///dataset: UPP///series: UPPT///id: 147B19///er: 0///pgr: 0///her2: NA///size: 2.4///node: 1///age: 71///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1855///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                             1 
 samplename: UPP_148B49///dataset: UPP///series: UPPT///id: 148B49///er: 1///pgr: 1///her2: NA///size: 2.6///node: 1///age: 71///grade: NA///t.dmfs: NA///e.dmfs: NA///t.rfs: 1065///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                             1 
      samplename: UPP_14B98///dataset: UPP///series: UPPT///id: 14B98///er: 1///pgr: 1///her2: NA///size: 2///node: 1///age: 61///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 3954///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                             1 
  samplename: UPP_150B81///dataset: UPP///series: UPPT///id: 150B81///er: 1///pgr: 1///her2: NA///size: 3.5///node: 1///age: 55///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 4167///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                             1 
  samplename: UPP_151B84///dataset: UPP///series: UPPU///id: 151B84///er: 1///pgr: 1///her2: NA///size: 1.5///node: 0///age: 57///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 4167///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
   samplename: UPP_152B99///dataset: UPP///series: UPPU///id: 152B99///er: 1///pgr: 1///her2: NA///size: 1.8///node: 0///age: 83///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 760///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_153B09///dataset: UPP///series: UPPU///id: 153B09///er: 0///pgr: 0///her2: NA///size: 1.8///node: 0///age: 71///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 4167///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_154B42///dataset: UPP///series: UPPT///id: 154B42///er: 1///pgr: 1///her2: NA///size: 0.9///node: 1///age: 52///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 1673///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                             1 
  samplename: UPP_155B52///dataset: UPP///series: UPPU///id: 155B52///er: 0///pgr: 0///her2: NA///size: 1.3///node: 0///age: 57///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 4167///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_156B01///dataset: UPP///series: UPPU///id: 156B01///er: 1///pgr: 1///her2: NA///size: 1.2///node: 0///age: 70///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 4167///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_158B84///dataset: UPP///series: UPPU///id: 158B84///er: 1///pgr: 1///her2: NA///size: 0.7///node: 0///age: 69///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 3103///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_159B47///dataset: UPP///series: UPPT///id: 159B47///er: 1///pgr: 1///her2: NA///size: 1.4///node: 1///age: 57///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 4167///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                             1 
      samplename: UPP_15C94///dataset: UPP///series: UPPT///id: 15C94///er: 1///pgr: 1///her2: NA///size: 5///node: 0///age: 83///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 1612///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                             1 
  samplename: UPP_160B16///dataset: UPP///series: UPPU///id: 160B16///er: 1///pgr: 1///her2: NA///size: 1.4///node: 0///age: 73///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2008///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_161B31///dataset: UPP///series: UPPU///id: 161B31///er: 1///pgr: 1///her2: NA///size: 1.3///node: 0///age: 69///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 4167///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_162B98///dataset: UPP///series: UPPT///id: 162B98///er: 1///pgr: 0///her2: NA///size: 3.3///node: 1///age: 73///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 3924///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                             1 
  samplename: UPP_163B27///dataset: UPP///series: UPPU///id: 163B27///er: 1///pgr: 1///her2: NA///size: 0.8///node: 0///age: 49///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 4137///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_164B81///dataset: UPP///series: UPPU///id: 164B81///er: 0///pgr: 0///her2: NA///size: 2.3///node: 0///age: 62///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 4137///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
   samplename: UPP_165B72///dataset: UPP///series: UPPT///id: 165B72///er: 1///pgr: 1///her2: NA///size: 2.3///node: 1///age: 74///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 852///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                             1 
  samplename: UPP_166B79///dataset: UPP///series: UPPT///id: 166B79///er: 1///pgr: 1///her2: NA///size: 2.2///node: 1///age: 65///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 4137///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                             1 
  samplename: UPP_168B51///dataset: UPP///series: UPPU///id: 168B51///er: 1///pgr: 1///her2: NA///size: 4.5///node: 0///age: 90///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 1947///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
samplename: UPP_169B79///dataset: UPP///series: UPPT///id: 169B79///er: 1///pgr: 1///her2: NA///size: 1.411///node: 0///age: 57///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 4137///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                             1 
    samplename: UPP_16C97///dataset: UPP///series: UPPU///id: 16C97///er: 1///pgr: 1///her2: NA///size: 1.5///node: 0///age: 80///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 1582///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_170B15///dataset: UPP///series: UPPT///id: 170B15///er: 1///pgr: 1///her2: NA///size: 2.5///node: 1///age: 70///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 1490///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                             1 
  samplename: UPP_171B77///dataset: UPP///series: UPPT///id: 171B77///er: 1///pgr: 1///her2: NA///size: 1.8///node: 0///age: 62///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2555///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                             1 
  samplename: UPP_172B19///dataset: UPP///series: UPPU///id: 172B19///er: 1///pgr: 1///her2: NA///size: 2.3///node: 0///age: 42///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 4106///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_173B43///dataset: UPP///series: UPPT///id: 173B43///er: 1///pgr: 1///her2: NA///size: 1.6///node: 1///age: 72///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 4106///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                             1 
  samplename: UPP_174B41///dataset: UPP///series: UPPU///id: 174B41///er: 1///pgr: 1///her2: NA///size: 1.1///node: 0///age: 74///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2099///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
    samplename: UPP_175B72///dataset: UPP///series: UPPU///id: 175B72///er: 1///pgr: 1///her2: NA///size: 6///node: 1///age: 46///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 152///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
  samplename: UPP_176B74///dataset: UPP///series: UPPU///id: 176B74///er: 1///pgr: 1///her2: NA///size: 2.5///node: 0///age: 77///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2190///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_177B67///dataset: UPP///series: UPPU///id: 177B67///er: 1///pgr: 1///her2: NA///size: 1.8///node: 0///age: 41///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 4076///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_178B74///dataset: UPP///series: UPPU///id: 178B74///er: 1///pgr: 1///her2: NA///size: 1.6///node: 0///age: 82///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2707///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_179B28///dataset: UPP///series: UPPU///id: 179B28///er: 1///pgr: 1///her2: NA///size: 0.8///node: 0///age: 66///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 1460///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
    samplename: UPP_17C40///dataset: UPP///series: UPPT///id: 17C40///er: 1///pgr: 1///her2: NA///size: 2.5///node: NA///age: 83///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 700///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                             1 
  samplename: UPP_180B38///dataset: UPP///series: UPPU///id: 180B38///er: 1///pgr: 1///her2: NA///size: 0.2///node: 0///age: 70///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2646///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
   samplename: UPP_181B70///dataset: UPP///series: UPPU///id: 181B70///er: 1///pgr: 1///her2: NA///size: 2///node: 1///age: 34///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 4045///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
 samplename: UPP_182B43///dataset: UPP///series: UPPU///id: 182B43///er: 0///pgr: 0///her2: NA///size: 2.8///node: 1///age: 51///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 4045///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
  samplename: UPP_183B75///dataset: UPP///series: UPPT///id: 183B75///er: 1///pgr: 1///her2: NA///size: 2.4///node: 1///age: 60///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 3954///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                             1 
    samplename: UPP_184B38///dataset: UPP///series: UPPU///id: 184B38///er: 1///pgr: 1///her2: NA///size: 1///node: 0///age: 63///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 4015///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
   samplename: UPP_185B44///dataset: UPP///series: UPPU///id: 185B44///er: 1///pgr: 1///her2: NA///size: 0.8///node: 0///age: 68///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
   samplename: UPP_186B22///dataset: UPP///series: UPPT///id: 186B22///er: 1///pgr: 1///her2: NA///size: 4.6///node: 1///age: 86///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 517///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                             1 
   samplename: UPP_187B36///dataset: UPP///series: UPPU///id: 187B36///er: 1///pgr: 1///her2: NA///size: 2.3///node: 1///age: 29///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
  samplename: UPP_188B13///dataset: UPP///series: UPPU///id: 188B13///er: 0///pgr: 0///her2: NA///size: 1.4///node: 0///age: 60///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 4015///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
   samplename: UPP_189B83///dataset: UPP///series: UPPT///id: 189B83///er: 0///pgr: 0///her2: NA///size: 3.8///node: 0///age: 58///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                             1 
    samplename: UPP_18C56///dataset: UPP///series: UPPU///id: 18C56///er: 1///pgr: 1///her2: NA///size: 1.9///node: 0///age: 68///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 3924///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_191B79///dataset: UPP///series: UPPT///id: 191B79///er: 1///pgr: 1///her2: NA///size: 2.7///node: 0///age: 59///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 4015///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                             1 
  samplename: UPP_192B69///dataset: UPP///series: UPPU///id: 192B69///er: 1///pgr: 1///her2: NA///size: 0.8///node: 0///age: 73///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2890///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_193B72///dataset: UPP///series: UPPT///id: 193B72///er: 0///pgr: 0///her2: NA///size: 1.6///node: 1///age: 68///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 3985///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                             1 
  samplename: UPP_194B60///dataset: UPP///series: UPPT///id: 194B60///er: 1///pgr: 0///her2: NA///size: 2.2///node: 1///age: 77///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1308///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                             1 
  samplename: UPP_195B75///dataset: UPP///series: UPPU///id: 195B75///er: 1///pgr: 1///her2: NA///size: 1.8///node: 0///age: 87///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 1095///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_196B81///dataset: UPP///series: UPPU///id: 196B81///er: 1///pgr: 1///her2: NA///size: 1.4///node: 0///age: 65///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 3985///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_197B95///dataset: UPP///series: UPPU///id: 197B95///er: 1///pgr: 1///her2: NA///size: 1.6///node: 0///age: 44///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 3985///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
    samplename: UPP_198B90///dataset: UPP///series: UPPU///id: 198B90///er: 1///pgr: 1///her2: NA///size: 2///node: 0///age: 74///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 3985///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_199B55///dataset: UPP///series: UPPU///id: 199B55///er: 1///pgr: 1///her2: NA///size: 2.3///node: 0///age: 54///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 3985///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
   samplename: UPP_19C33///dataset: UPP///series: UPPU///id: 19C33///er: 1///pgr: 1///her2: NA///size: 3.1///node: 1///age: 61///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 1521///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
   samplename: UPP_200B47///dataset: UPP///series: UPPU///id: 200B47///er: 1///pgr: 1///her2: NA///size: 2///node: 1///age: 47///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 3985///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
 samplename: UPP_201B68///dataset: UPP///series: UPPU///id: 201B68///er: 1///pgr: 1///her2: NA///size: 2.4///node: 1///age: 53///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 3985///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: NA 
                                                                                                                                                                                                                                             1 
    samplename: UPP_202B44///dataset: UPP///series: UPPT///id: 202B44///er: 1///pgr: 1///her2: NA///size: 1///node: 0///age: 72///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 3954///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                             1 
    samplename: UPP_203B49///dataset: UPP///series: UPPT///id: 203B49///er: 1///pgr: 1///her2: NA///size: 4///node: 1///age: 51///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 3954///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                             1 
  samplename: UPP_204B85///dataset: UPP///series: UPPU///id: 204B85///er: 1///pgr: 1///her2: NA///size: 1.2///node: 0///age: 70///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 3863///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_205B99///dataset: UPP///series: UPPU///id: 205B99///er: 1///pgr: 1///her2: NA///size: 2.2///node: 0///age: 59///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 1612///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_206C05///dataset: UPP///series: UPPU///id: 206C05///er: 1///pgr: 1///her2: NA///size: 1.8///node: 0///age: 80///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2342///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
    samplename: UPP_207C08///dataset: UPP///series: UPPT///id: 207C08///er: 1///pgr: 1///her2: NA///size: 1///node: 1///age: 59///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 3954///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                             1 
   samplename: UPP_208C06///dataset: UPP///series: UPPU///id: 208C06///er: 1///pgr: 1///her2: NA///size: 2.2///node: NA///age: 67///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 30///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
  samplename: UPP_209C10///dataset: UPP///series: UPPU///id: 209C10///er: 1///pgr: 1///her2: NA///size: 1.4///node: 0///age: 80///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 3954///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                             1 
                                                                                                                                                                                                                                       (Other) 
                                                                                                                                                                                                                                           152 

}}
\keyword{datasets}
