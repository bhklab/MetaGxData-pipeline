\name{ SUPERTAM_HGU133A_eset }
\alias{ SUPERTAM_HGU133A_eset }
\docType{data}
\title{ SUPERTAM_HGU133A_eset }
\usage{data( SUPERTAM_HGU133A_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www-ncbi-nlm-nih-gov.proxy.wexler.hunter.cuny.edu/pubmed/?term=18498629 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6532 
  PMIDs: 18498629 
  No abstract available.
  notes:
   summary:     
      A gene classifier was identified as a predictor of clinical outcome in tamoxifen-treated breast cancer. 
   version:     
      2015-04-27 19:36:00
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
assayData: 21169 features, 856 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
    KIT_104B91     KIT_112B55     KIT_114B68     KIT_130B92     KIT_131B79     KIT_135B40     KIT_138B34     KIT_139B03 
             1              1              1              1              1              1              1              1 
    KIT_148B49      KIT_14B98     KIT_150B81     KIT_154B42     KIT_159B47      KIT_15C94     KIT_162B98     KIT_165B72 
             1              1              1              1              1              1              1              1 
    KIT_166B79     KIT_169B79     KIT_170B15     KIT_171B77     KIT_173B43      KIT_17C40     KIT_183B75     KIT_186B22 
             1              1              1              1              1              1              1              1 
    KIT_191B79     KIT_194B60     KIT_202B44     KIT_203B49     KIT_207C08     KIT_212C21     KIT_216C61     KIT_225C52 
             1              1              1              1              1              1              1              1 
    KIT_230C47     KIT_235C20     KIT_238C87     KIT_239C32     KIT_244C89     KIT_252C64     KIT_254C80     KIT_256C45 
             1              1              1              1              1              1              1              1 
    KIT_271C71      KIT_27C82     KIT_288C57     KIT_289C75     KIT_291C17     KIT_298C47     KIT_301C66     KIT_307C50 
             1              1              1              1              1              1              1              1 
    KIT_311A27     KIT_313A87      KIT_35C29      KIT_39C24      KIT_40C57      KIT_41C65      KIT_48A46      KIT_51A98 
             1              1              1              1              1              1              1              1 
     KIT_53A06      KIT_55A79      KIT_56A94      KIT_63A62      KIT_64A59      KIT_67A43      KIT_69A93       KIT_6B85 
             1              1              1              1              1              1              1              1 
     KIT_72A92      KIT_77A50       KIT_7B96      KIT_82A83      KIT_84A44      KIT_85A03       KIT_8B87      KIT_96A21 
             1              1              1              1              1              1              1              1 
 MDA5_AUST_A01  MDA5_AUST_A02  MDA5_AUST_A03  MDA5_AUST_A04  MDA5_AUST_A05  MDA5_AUST_A06  MDA5_AUST_A07  MDA5_AUST_A08 
             1              1              1              1              1              1              1              1 
 MDA5_AUST_A09  MDA5_AUST_A10 MDA5_AUST_A100 MDA5_AUST_A101  MDA5_AUST_A11  MDA5_AUST_A12  MDA5_AUST_A13  MDA5_AUST_A14 
             1              1              1              1              1              1              1              1 
 MDA5_AUST_A15  MDA5_AUST_A16  MDA5_AUST_A17  MDA5_AUST_A18  MDA5_AUST_A19  MDA5_AUST_A20  MDA5_AUST_A21  MDA5_AUST_A22 
             1              1              1              1              1              1              1              1 
 MDA5_AUST_A23  MDA5_AUST_A24  MDA5_AUST_A25        (Other) 
             1              1              1            757 

sample_type: 
tumor 
  856 

er: 
negative positive 
     140      716 

pgr: 
negative positive     NA's 
       5       65      786 

tumor_size: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  0.200   1.000   2.000   2.146   2.500   9.900     224 

N: 
   0    1 NA's 
 677  165   14 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  26.00   51.00   61.00   59.97   69.00   90.70     194 

grade: 
   1    2    3 NA's 
  84  245  219  308 

dmfs_days: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
      8    1586    2650    2553    3468    5939 

dmfs_status: 
deceased_or_recurrence    living_norecurrence                   NA's 
                   232                    618                      6 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    8.0   646.2  1568.0  1910.0  3398.0  4289.0     778 

recurrence_status: 
deceased_or_recurrence    living_norecurrence                   NA's 
                    34                     44                    778 

treatment: 
     endocrine hormonotherapy      untreated 
           298            214            344 

batch: 
AUST IGRT  KIT OXFT  VDX VDX3 VDXN 
  99   96   72  109  286  136   58 

uncurated_author_metadata: 
    samplename: KIT_104B91///dataset: TAM///series: KIT///id: 104B91///er: 1///pgr: 1///her2: NA///size: 2.4///node: NA///age: 80///grade: 3///t.dmfs: 1308///e.dmfs: 0///t.rfs: 1307///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
       samplename: KIT_112B55///dataset: TAM///series: KIT///id: 112B55///er: 1///pgr: 1///her2: NA///size: 2.4///node: 1///age: 61///grade: 2///t.dmfs: 335///e.dmfs: 1///t.rfs: 336///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
     samplename: KIT_114B68///dataset: TAM///series: KIT///id: 114B68///er: 1///pgr: 1///her2: NA///size: 1.2///node: 1///age: 67///grade: 1///t.dmfs: 3833///e.dmfs: 0///t.rfs: 3832///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
     samplename: KIT_130B92///dataset: TAM///series: KIT///id: 130B92///er: 1///pgr: 1///her2: NA///size: 1.8///node: 1///age: 73///grade: 2///t.dmfs: 1612///e.dmfs: 1///t.rfs: 1613///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
       samplename: KIT_131B79///dataset: TAM///series: KIT///id: 131B79///er: 1///pgr: 1///her2: NA///size: 2.6///node: 1///age: 59///grade: 3///t.dmfs: 639///e.dmfs: 1///t.rfs: 639///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
     samplename: KIT_135B40///dataset: TAM///series: KIT///id: 135B40///er: 1///pgr: 1///her2: NA///size: 2.6///node: 0///age: 64///grade: 1///t.dmfs: 3468///e.dmfs: 1///t.rfs: 3468///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
     samplename: KIT_138B34///dataset: TAM///series: KIT///id: 138B34///er: 1///pgr: 1///her2: NA///size: 2.3///node: 1///age: 65///grade: 1///t.dmfs: 3468///e.dmfs: 0///t.rfs: 3468///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
     samplename: KIT_139B03///dataset: TAM///series: KIT///id: 139B03///er: 1///pgr: 1///her2: NA///size: 2.8///node: 0///age: 84///grade: NA///t.dmfs: 122///e.dmfs: NA///t.rfs: 120///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
    samplename: KIT_148B49///dataset: TAM///series: KIT///id: 148B49///er: 1///pgr: 1///her2: NA///size: 2.6///node: 1///age: 71///grade: NA///t.dmfs: 1034///e.dmfs: 0///t.rfs: 1033///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
         samplename: KIT_14B98///dataset: TAM///series: KIT///id: 14B98///er: 1///pgr: 1///her2: NA///size: 2///node: 1///age: 61///grade: 2///t.dmfs: 3680///e.dmfs: 0///t.rfs: 3679///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
     samplename: KIT_150B81///dataset: TAM///series: KIT///id: 150B81///er: 1///pgr: 1///her2: NA///size: 3.5///node: 1///age: 55///grade: 2///t.dmfs: 3498///e.dmfs: 0///t.rfs: 3497///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
     samplename: KIT_154B42///dataset: TAM///series: KIT///id: 154B42///er: 1///pgr: 1///her2: NA///size: 0.9///node: 1///age: 52///grade: 2///t.dmfs: 1247///e.dmfs: 1///t.rfs: 1248///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
     samplename: KIT_159B47///dataset: TAM///series: KIT///id: 159B47///er: 1///pgr: 1///her2: NA///size: 1.4///node: 1///age: 57///grade: 2///t.dmfs: 2373///e.dmfs: 1///t.rfs: 2372///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
         samplename: KIT_15C94///dataset: TAM///series: KIT///id: 15C94///er: 1///pgr: 1///her2: NA///size: 5///node: 0///age: 83///grade: 2///t.dmfs: 1612///e.dmfs: 0///t.rfs: 1613///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
     samplename: KIT_162B98///dataset: TAM///series: KIT///id: 162B98///er: 1///pgr: 1///her2: NA///size: 3.3///node: 1///age: 73///grade: 3///t.dmfs: 3255///e.dmfs: 1///t.rfs: 3256///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
       samplename: KIT_165B72///dataset: TAM///series: KIT///id: 165B72///er: 1///pgr: 1///her2: NA///size: 2.3///node: 1///age: 74///grade: 2///t.dmfs: 548///e.dmfs: 1///t.rfs: 548///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
     samplename: KIT_166B79///dataset: TAM///series: KIT///id: 166B79///er: 1///pgr: 1///her2: NA///size: 2.2///node: 1///age: 65///grade: 2///t.dmfs: 3589///e.dmfs: 0///t.rfs: 3588///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
   samplename: KIT_169B79///dataset: TAM///series: KIT///id: 169B79///er: 1///pgr: 1///her2: NA///size: 1.411///node: 0///age: 57///grade: 2///t.dmfs: 3589///e.dmfs: 0///t.rfs: 3588///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
     samplename: KIT_170B15///dataset: TAM///series: KIT///id: 170B15///er: 1///pgr: 1///her2: NA///size: 2.5///node: 1///age: 70///grade: 2///t.dmfs: 1490///e.dmfs: 0///t.rfs: 1489///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
       samplename: KIT_171B77///dataset: TAM///series: KIT///id: 171B77///er: 1///pgr: 1///her2: NA///size: 1.8///node: 0///age: 62///grade: 2///t.dmfs: 639///e.dmfs: 0///t.rfs: 639///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
     samplename: KIT_173B43///dataset: TAM///series: KIT///id: 173B43///er: 1///pgr: 1///her2: NA///size: 1.6///node: 1///age: 72///grade: 1///t.dmfs: 3954///e.dmfs: 0///t.rfs: 3953///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
        samplename: KIT_17C40///dataset: TAM///series: KIT///id: 17C40///er: 1///pgr: 1///her2: NA///size: 2.5///node: NA///age: 83///grade: 2///t.dmfs: 700///e.dmfs: 0///t.rfs: 701///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
     samplename: KIT_183B75///dataset: TAM///series: KIT///id: 183B75///er: 1///pgr: 1///her2: NA///size: 2.4///node: 1///age: 60///grade: 2///t.dmfs: 2555///e.dmfs: 1///t.rfs: 2555///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
       samplename: KIT_186B22///dataset: TAM///series: KIT///id: 186B22///er: 1///pgr: 1///her2: NA///size: 4.6///node: 1///age: 86///grade: NA///t.dmfs: 61///e.dmfs: NA///t.rfs: 62///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
     samplename: KIT_191B79///dataset: TAM///series: KIT///id: 191B79///er: 1///pgr: 1///her2: NA///size: 2.7///node: 0///age: 59///grade: 2///t.dmfs: 1612///e.dmfs: 0///t.rfs: 1613///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
     samplename: KIT_194B60///dataset: TAM///series: KIT///id: 194B60///er: 1///pgr: 0///her2: NA///size: 2.2///node: 1///age: 77///grade: 3///t.dmfs: 1247///e.dmfs: 0///t.rfs: 1248///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
       samplename: KIT_202B44///dataset: TAM///series: KIT///id: 202B44///er: 1///pgr: 1///her2: NA///size: 1///node: 0///age: 72///grade: 2///t.dmfs: 3589///e.dmfs: 0///t.rfs: 3588///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
       samplename: KIT_203B49///dataset: TAM///series: KIT///id: 203B49///er: 1///pgr: 1///her2: NA///size: 4///node: 1///age: 51///grade: 2///t.dmfs: 3620///e.dmfs: 0///t.rfs: 3621///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
       samplename: KIT_207C08///dataset: TAM///series: KIT///id: 207C08///er: 1///pgr: 1///her2: NA///size: 1///node: 1///age: 59///grade: 1///t.dmfs: 3407///e.dmfs: 0///t.rfs: 3405///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
    samplename: KIT_212C21///dataset: TAM///series: KIT///id: 212C21///er: 1///pgr: NA///her2: NA///size: 2.4///node: 0///age: 74///grade: 2///t.dmfs: 1369///e.dmfs: 0///t.rfs: 1369///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
     samplename: KIT_216C61///dataset: TAM///series: KIT///id: 216C61///er: 1///pgr: 1///her2: NA///size: 2.2///node: 1///age: 71///grade: 2///t.dmfs: 3285///e.dmfs: 0///t.rfs: 3285///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
     samplename: KIT_225C52///dataset: TAM///series: KIT///id: 225C52///er: 1///pgr: 1///her2: NA///size: 2.7///node: 1///age: 72///grade: 2///t.dmfs: 3741///e.dmfs: 0///t.rfs: 3741///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
       samplename: KIT_230C47///dataset: TAM///series: KIT///id: 230C47///er: 1///pgr: 1///her2: NA///size: 2.2///node: 1///age: 77///grade: 2///t.dmfs: 183///e.dmfs: 1///t.rfs: 182///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
    samplename: KIT_235C20///dataset: TAM///series: KIT///id: 235C20///er: 1///pgr: NA///her2: NA///size: 2.1///node: 0///age: 71///grade: 1///t.dmfs: 3528///e.dmfs: 0///t.rfs: 3530///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
     samplename: KIT_238C87///dataset: TAM///series: KIT///id: 238C87///er: 0///pgr: 1///her2: NA///size: 2.4///node: 0///age: 54///grade: 3///t.dmfs: 3395///e.dmfs: 0///t.rfs: 3376///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
     samplename: KIT_239C32///dataset: TAM///series: KIT///id: 239C32///er: 1///pgr: 1///her2: NA///size: 1.7///node: 1///age: 65///grade: 2///t.dmfs: 3468///e.dmfs: 0///t.rfs: 3468///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
     samplename: KIT_244C89///dataset: TAM///series: KIT///id: 244C89///er: 1///pgr: 1///her2: NA///size: 3.3///node: 1///age: 51///grade: 2///t.dmfs: 2646///e.dmfs: 1///t.rfs: 2646///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
   samplename: KIT_252C64///dataset: TAM///series: KIT///id: 252C64///er: 1///pgr: 1///her2: NA///size: 2.525///node: 0///age: 77///grade: 3///t.dmfs: 1125///e.dmfs: 1///t.rfs: 1124///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
     samplename: KIT_254C80///dataset: TAM///series: KIT///id: 254C80///er: 1///pgr: 1///her2: NA///size: 2.2///node: 0///age: 67///grade: 2///t.dmfs: 3437///e.dmfs: 0///t.rfs: 3438///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
       samplename: KIT_256C45///dataset: TAM///series: KIT///id: 256C45///er: 0///pgr: 1///her2: NA///size: 3.2///node: 1///age: 57///grade: 2///t.dmfs: 475///e.dmfs: 1///t.rfs: 456///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
       samplename: KIT_271C71///dataset: TAM///series: KIT///id: 271C71///er: 0///pgr: 1///her2: NA///size: 2.9///node: 1///age: 66///grade: 3///t.dmfs: 438///e.dmfs: 1///t.rfs: 427///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
       samplename: KIT_27C82///dataset: TAM///series: KIT///id: 27C82///er: 1///pgr: 1///her2: NA///size: 1.3///node: 1///age: 67///grade: 2///t.dmfs: 2160///e.dmfs: 0///t.rfs: 2161///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
     samplename: KIT_288C57///dataset: TAM///series: KIT///id: 288C57///er: 1///pgr: 1///her2: NA///size: 3.8///node: 1///age: 88///grade: 2///t.dmfs: 2312///e.dmfs: 0///t.rfs: 2310///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
       samplename: KIT_289C75///dataset: TAM///series: KIT///id: 289C75///er: 1///pgr: 1///her2: NA///size: 2///node: 1///age: 52///grade: 1///t.dmfs: 3072///e.dmfs: 0///t.rfs: 3073///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
       samplename: KIT_291C17///dataset: TAM///series: KIT///id: 291C17///er: 0///pgr: 1///her2: NA///size: 4.2///node: 0///age: 84///grade: 3///t.dmfs: 329///e.dmfs: 0///t.rfs: 336///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
     samplename: KIT_298C47///dataset: TAM///series: KIT///id: 298C47///er: 1///pgr: 1///her2: NA///size: 2.4///node: 0///age: 62///grade: 2///t.dmfs: 2373///e.dmfs: 1///t.rfs: 2372///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
     samplename: KIT_301C66///dataset: TAM///series: KIT///id: 301C66///er: 1///pgr: 1///her2: NA///size: 3.2///node: 0///age: 74///grade: 2///t.dmfs: 1186///e.dmfs: 0///t.rfs: 1186///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
     samplename: KIT_307C50///dataset: TAM///series: KIT///id: 307C50///er: 1///pgr: 1///her2: NA///size: 2.4///node: 0///age: 66///grade: 2///t.dmfs: 2859///e.dmfs: 0///t.rfs: 2858///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
       samplename: KIT_311A27///dataset: TAM///series: KIT///id: 311A27///er: 1///pgr: 0///her2: NA///size: 3.2///node: 1///age: 70///grade: 3///t.dmfs: 213///e.dmfs: 0///t.rfs: 212///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
     samplename: KIT_313A87///dataset: TAM///series: KIT///id: 313A87///er: 1///pgr: 0///her2: NA///size: 2.8///node: 1///age: 65///grade: 3///t.dmfs: 3407///e.dmfs: 0///t.rfs: 3405///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
         samplename: KIT_35C29///dataset: TAM///series: KIT///id: 35C29///er: 1///pgr: 1///her2: NA///size: 2.3///node: 0///age: 45///grade: 2///t.dmfs: 882///e.dmfs: 0///t.rfs: 883///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
       samplename: KIT_39C24///dataset: TAM///series: KIT///id: 39C24///er: 1///pgr: 1///her2: NA///size: 1.7///node: 1///age: 52///grade: 1///t.dmfs: 3315///e.dmfs: 0///t.rfs: 3314///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
       samplename: KIT_40C57///dataset: TAM///series: KIT///id: 40C57///er: 1///pgr: 1///her2: NA///size: 2.2///node: 1///age: 56///grade: 2///t.dmfs: 3194///e.dmfs: 0///t.rfs: 3194///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
       samplename: KIT_41C65///dataset: TAM///series: KIT///id: 41C65///er: 1///pgr: 1///her2: NA///size: 2.2///node: 1///age: 61///grade: 2///t.dmfs: 3133///e.dmfs: 0///t.rfs: 3132///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
        samplename: KIT_48A46///dataset: TAM///series: KIT///id: 48A46///er: 1///pgr: 1///her2: NA///size: 3.8///node: NA///age: 78///grade: 1///t.dmfs: 669///e.dmfs: 0///t.rfs: 668///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
       samplename: KIT_51A98///dataset: TAM///series: KIT///id: 51A98///er: 1///pgr: 1///her2: NA///size: 3.2///node: 1///age: 59///grade: 2///t.dmfs: 3680///e.dmfs: 0///t.rfs: 3679///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
         samplename: KIT_53A06///dataset: TAM///series: KIT///id: 53A06///er: 1///pgr: 1///her2: NA///size: 2.1///node: 1///age: 82///grade: 2///t.dmfs: 943///e.dmfs: 0///t.rfs: 942///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
       samplename: KIT_55A79///dataset: TAM///series: KIT///id: 55A79///er: 1///pgr: 1///her2: NA///size: 4.1///node: 1///age: 62///grade: 3///t.dmfs: 1186///e.dmfs: 0///t.rfs: 1186///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
         samplename: KIT_56A94///dataset: TAM///series: KIT///id: 56A94///er: 1///pgr: 1///her2: NA///size: 4.5///node: 1///age: 52///grade: 2///t.dmfs: 395///e.dmfs: 1///t.rfs: 394///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
          samplename: KIT_63A62///dataset: TAM///series: KIT///id: 63A62///er: 1///pgr: 1///her2: NA///size: 2.4///node: 1///age: 55///grade: 2///t.dmfs: 61///e.dmfs: NA///t.rfs: 62///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
       samplename: KIT_64A59///dataset: TAM///series: KIT///id: 64A59///er: 1///pgr: 0///her2: NA///size: 2.8///node: 0///age: 61///grade: 2///t.dmfs: 3650///e.dmfs: 0///t.rfs: 3650///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
         samplename: KIT_67A43///dataset: TAM///series: KIT///id: 67A43///er: 0///pgr: 1///her2: NA///size: 2.8///node: 1///age: 79///grade: 3///t.dmfs: 329///e.dmfs: 1///t.rfs: 336///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
         samplename: KIT_69A93///dataset: TAM///series: KIT///id: 69A93///er: 1///pgr: 1///her2: NA///size: 5///node: 1///age: 62///grade: 3///t.dmfs: 4289///e.dmfs: 0///t.rfs: 4289///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
           samplename: KIT_6B85///dataset: TAM///series: KIT///id: 6B85///er: 1///pgr: 1///her2: NA///size: 2.6///node: 1///age: 71///grade: 1///t.dmfs: 213///e.dmfs: 1///t.rfs: 212///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
       samplename: KIT_72A92///dataset: TAM///series: KIT///id: 72A92///er: 1///pgr: 1///her2: NA///size: 1.6///node: 0///age: 55///grade: 1///t.dmfs: 3072///e.dmfs: 0///t.rfs: 3073///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
           samplename: KIT_77A50///dataset: TAM///series: KIT///id: 77A50///er: 1///pgr: 0///her2: NA///size: 2///node: 1///age: 75///grade: 2///t.dmfs: 395///e.dmfs: 1///t.rfs: 394///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
           samplename: KIT_7B96///dataset: TAM///series: KIT///id: 7B96///er: 1///pgr: 1///her2: NA///size: 1.9///node: 1///age: 52///grade: 2///t.dmfs: 882///e.dmfs: 1///t.rfs: 883///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
         samplename: KIT_82A83///dataset: TAM///series: KIT///id: 82A83///er: 1///pgr: 1///her2: NA///size: 2.7///node: 1///age: 72///grade: 3///t.dmfs: 943///e.dmfs: 0///t.rfs: 942///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
       samplename: KIT_84A44///dataset: TAM///series: KIT///id: 84A44///er: 1///pgr: 1///her2: NA///size: 2.8///node: 0///age: 84///grade: 2///t.dmfs: 2281///e.dmfs: 0///t.rfs: 2281///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
         samplename: KIT_85A03///dataset: TAM///series: KIT///id: 85A03///er: 1///pgr: 1///her2: NA///size: 1.5///node: 1///age: 86///grade: 2///t.dmfs: 760///e.dmfs: 0///t.rfs: 759///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
         samplename: KIT_8B87///dataset: TAM///series: KIT///id: 8B87///er: 1///pgr: 1///her2: NA///size: 1.7///node: 0///age: 58///grade: 1///t.dmfs: 3620///e.dmfs: 0///t.rfs: 3621///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
           samplename: KIT_96A21///dataset: TAM///series: KIT///id: 96A21///er: 1///pgr: 1///her2: NA///size: 3.8///node: 1///age: 63///grade: 2///t.dmfs: 91///e.dmfs: 1///t.rfs: 29///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                                1 
 samplename: MDA5_AUST_A01///dataset: MDA5///series: AUST///id: A01///er: 1///pgr: NA///her2: NA///size: 2.5///node: 1///age: 41.3///grade: 3///t.dmfs: 3876///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 5 
                                                                                                                                                                                                                                                1 
 samplename: MDA5_AUST_A02///dataset: MDA5///series: AUST///id: A02///er: 1///pgr: NA///her2: NA///size: 2.8///node: 0///age: 69.1///grade: 1///t.dmfs: 1278///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 5 
                                                                                                                                                                                                                                                1 
 samplename: MDA5_AUST_A03///dataset: MDA5///series: AUST///id: A03///er: 1///pgr: NA///her2: NA///size: 9.9///node: 1///age: 68.1///grade: 1///t.dmfs: 1701///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 5 
                                                                                                                                                                                                                                                1 
 samplename: MDA5_AUST_A04///dataset: MDA5///series: AUST///id: A04///er: 1///pgr: NA///her2: NA///size: 2.5///node: 1///age: 58.3///grade: 1///t.dmfs: 5555///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 5 
                                                                                                                                                                                                                                                1 
 samplename: MDA5_AUST_A05///dataset: MDA5///series: AUST///id: A05///er: 1///pgr: NA///her2: NA///size: 1.2///node: 0///age: 62.1///grade: 2///t.dmfs: 5267///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 5 
                                                                                                                                                                                                                                                1 
  samplename: MDA5_AUST_A06///dataset: MDA5///series: AUST///id: A06///er: 1///pgr: NA///her2: NA///size: 2.5///node: 1///age: 68///grade: NA///t.dmfs: 5084///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 5 
                                                                                                                                                                                                                                                1 
   samplename: MDA5_AUST_A07///dataset: MDA5///series: AUST///id: A07///er: 1///pgr: NA///her2: NA///size: 1.6///node: 1///age: 58///grade: 1///t.dmfs: 2402///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 5 
                                                                                                                                                                                                                                                1 
   samplename: MDA5_AUST_A08///dataset: MDA5///series: AUST///id: A08///er: 1///pgr: NA///her2: NA///size: 3.5///node: 1///age: 53///grade: 2///t.dmfs: 1909///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 5 
                                                                                                                                                                                                                                                1 
  samplename: MDA5_AUST_A09///dataset: MDA5///series: AUST///id: A09///er: 1///pgr: NA///her2: NA///size: 5///node: NA///age: 67.4///grade: 2///t.dmfs: 4603///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 5 
                                                                                                                                                                                                                                                1 
samplename: MDA5_AUST_A10///dataset: MDA5///series: AUST///id: A10///er: 1///pgr: NA///her2: NA///size: 0.8///node: NA///age: 52.8///grade: 1///t.dmfs: 4281///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 5 
                                                                                                                                                                                                                                                1 
  samplename: MDA5_AUST_A100///dataset: MDA5///series: AUST///id: A100///er: 1///pgr: NA///her2: NA///size: NA///node: 0///age: 67///grade: 3///t.dmfs: 5267///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 5 
                                                                                                                                                                                                                                                1 
  samplename: MDA5_AUST_A101///dataset: MDA5///series: AUST///id: A101///er: 1///pgr: NA///her2: NA///size: NA///node: 1///age: 56///grade: 2///t.dmfs: 5939///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 5 
                                                                                                                                                                                                                                                1 
 samplename: MDA5_AUST_A11///dataset: MDA5///series: AUST///id: A11///er: 1///pgr: NA///her2: NA///size: 9.9///node: 1///age: 45.8///grade: 3///t.dmfs: 4395///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 5 
                                                                                                                                                                                                                                                1 
   samplename: MDA5_AUST_A12///dataset: MDA5///series: AUST///id: A12///er: 1///pgr: NA///her2: NA///size: 3///node: 0///age: 70.4///grade: 1///t.dmfs: 2949///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 5 
                                                                                                                                                                                                                                                1 
 samplename: MDA5_AUST_A13///dataset: MDA5///series: AUST///id: A13///er: 1///pgr: NA///her2: NA///size: 2.1///node: 0///age: 76.3///grade: 2///t.dmfs: 4249///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 5 
                                                                                                                                                                                                                                                1 
 samplename: MDA5_AUST_A14///dataset: MDA5///series: AUST///id: A14///er: 1///pgr: NA///her2: NA///size: 2.4///node: 0///age: 64.6///grade: 2///t.dmfs: 4457///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 5 
                                                                                                                                                                                                                                                1 
 samplename: MDA5_AUST_A15///dataset: MDA5///series: AUST///id: A15///er: 1///pgr: NA///her2: NA///size: 2.5///node: 0///age: 42.3///grade: 1///t.dmfs: 4358///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 5 
                                                                                                                                                                                                                                                1 
 samplename: MDA5_AUST_A16///dataset: MDA5///series: AUST///id: A16///er: 1///pgr: NA///her2: NA///size: 2.2///node: 1///age: 79.9///grade: 1///t.dmfs: 2727///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 5 
                                                                                                                                                                                                                                                1 
   samplename: MDA5_AUST_A17///dataset: MDA5///series: AUST///id: A17///er: 1///pgr: NA///her2: NA///size: 2///node: 0///age: 64.6///grade: 2///t.dmfs: 2730///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 5 
                                                                                                                                                                                                                                                1 
   samplename: MDA5_AUST_A18///dataset: MDA5///series: AUST///id: A18///er: 1///pgr: NA///her2: NA///size: 6///node: 0///age: 77.2///grade: 1///t.dmfs: 3916///e.dmfs: 1///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 5 
                                                                                                                                                                                                                                                1 
   samplename: MDA5_AUST_A19///dataset: MDA5///series: AUST///id: A19///er: 1///pgr: NA///her2: NA///size: 9.9///node: 0///age: 70///grade: 2///t.dmfs: 3989///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 5 
                                                                                                                                                                                                                                                1 
 samplename: MDA5_AUST_A20///dataset: MDA5///series: AUST///id: A20///er: 1///pgr: NA///her2: NA///size: 1.6///node: 1///age: 43.3///grade: 1///t.dmfs: 3792///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 5 
                                                                                                                                                                                                                                                1 
 samplename: MDA5_AUST_A21///dataset: MDA5///series: AUST///id: A21///er: 1///pgr: NA///her2: NA///size: 3.5///node: 0///age: 52.5///grade: 2///t.dmfs: 4117///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 5 
                                                                                                                                                                                                                                                1 
   samplename: MDA5_AUST_A22///dataset: MDA5///series: AUST///id: A22///er: 1///pgr: NA///her2: NA///size: 2///node: 0///age: 50.1///grade: 1///t.dmfs: 3460///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 5 
                                                                                                                                                                                                                                                1 
 samplename: MDA5_AUST_A23///dataset: MDA5///series: AUST///id: A23///er: 1///pgr: NA///her2: NA///size: 0.7///node: 0///age: 62.7///grade: 2///t.dmfs: 3687///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 5 
                                                                                                                                                                                                                                                1 
 samplename: MDA5_AUST_A24///dataset: MDA5///series: AUST///id: A24///er: 1///pgr: NA///her2: NA///size: 9.9///node: 0///age: 54.8///grade: 2///t.dmfs: 3405///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 5 
                                                                                                                                                                                                                                                1 
   samplename: MDA5_AUST_A25///dataset: MDA5///series: AUST///id: A25///er: 1///pgr: NA///her2: NA///size: 2///node: 1///age: 36.9///grade: 2///t.dmfs: 3314///e.dmfs: 0///t.rfs: NA///e.rfs: NA///t.os: NA///e.os: NA///tissue: 1///treatment: 5 
                                                                                                                                                                                                                                                1 
                                                                                                                                                                                                                                          (Other) 
                                                                                                                                                                                                                                              757 

duplicates: 
SUPERTAM_HGU133A.SUPERTAM_HGU133A_MDA5_AUST_A53 SUPERTAM_HGU133A.SUPERTAM_HGU133A_MDA5_AUST_A92 
                                              1                                               1 
                                           NA's 
                                            854 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6532 }
\keyword{datasets}
