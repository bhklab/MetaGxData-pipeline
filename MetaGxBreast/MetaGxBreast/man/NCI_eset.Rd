\name{ NCI_eset }
\alias{ NCI_eset }
\docType{data}
\title{  }
\usage{data( NCI_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=12490681; http://www.ncbi.nlm.nih.gov/pubmed/?term=11823860 
  Title:  
  URL: Not available 
  PMIDs: 12490681, 11823860  
  No abstract available.
  notes:
   summary:     
      It was found that the gene expression profile that was studied was more powerful in predicting outcome of disease in youn
ger patients than using standard clinical and pathological criteria.
   version:     
      2015-04-27 19:23:18
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: AF106966 AF217974 ... Y12473 (5154 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 5154 features, 99 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
NCI_21585 NCI_21586 NCI_21588 NCI_21589 NCI_21591 NCI_21592 NCI_21593 NCI_21594 NCI_21595 NCI_21596 NCI_21597 NCI_21598 
        1         1         1         1         1         1         1         1         1         1         1         1 
NCI_21599 NCI_21600 NCI_21603 NCI_21604 NCI_21605 NCI_21606 NCI_21607 NCI_21608 NCI_21609 NCI_21610 NCI_21611 NCI_21612 
        1         1         1         1         1         1         1         1         1         1         1         1 
NCI_21613 NCI_21614 NCI_21615 NCI_21616 NCI_21617 NCI_21618 NCI_21619 NCI_21620 NCI_21621 NCI_21622 NCI_21623 NCI_21624 
        1         1         1         1         1         1         1         1         1         1         1         1 
NCI_21625 NCI_21626 NCI_21627 NCI_21628 NCI_21629 NCI_21630 NCI_21631 NCI_21632 NCI_21633 NCI_21634 NCI_21635 NCI_21636 
        1         1         1         1         1         1         1         1         1         1         1         1 
NCI_21637 NCI_21638 NCI_21639 NCI_21640 NCI_21641 NCI_21642 NCI_21643 NCI_21644 NCI_21645 NCI_21647 NCI_21648 NCI_21649 
        1         1         1         1         1         1         1         1         1         1         1         1 
NCI_21650 NCI_21652 NCI_21653 NCI_21655 NCI_21656 NCI_21657 NCI_21658 NCI_21659 NCI_21660 NCI_21661 NCI_21662 NCI_21663 
        1         1         1         1         1         1         1         1         1         1         1         1 
NCI_21664 NCI_21665 NCI_21666 NCI_21667 NCI_21668 NCI_21669 NCI_21670 NCI_21671 NCI_21672 NCI_21674 NCI_21675 NCI_21676 
        1         1         1         1         1         1         1         1         1         1         1         1 
NCI_21677 NCI_21678 NCI_21679 NCI_21680 NCI_21681 NCI_21682 NCI_21823 NCI_21824 NCI_21825 NCI_21826 NCI_21827 NCI_21828 
        1         1         1         1         1         1         1         1         1         1         1         1 
NCI_21829 NCI_21830 NCI_21831 
        1         1         1 

alt_sample_name: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  21580   21610   21640   21650   21670   21830 

sample_type: 
tumor 
   99 

er: 
negative positive 
      34       65 

tumor_size: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   0.80    1.80    2.50    2.82    3.00    8.00 

N: 
 0  1 
46 53 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  33.00   49.00   57.00   57.47   64.50   90.00 

grade: 
 1  2  3 
16 38 45 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
      8     967    2057    1969    2930    4067 

recurrence_status: 
deceased_or_recurrence    living_norecurrence 
                    45                     54 

treatment: 
  chemotherapy hormonotherapy      untreated 
            10             78             11 

batch: 
NCI 
 99 

uncurated_author_metadata: 
samplename: NCI_21585///dataset: NCI///series: NCI///id: 21585///er: 0///pgr: NA///her2: NA///size: 2.7///node: 1///age: 65///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 4067///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                         1 
samplename: NCI_21586///dataset: NCI///series: NCI///id: 21586///er: 1///pgr: NA///her2: NA///size: 2.1///node: 0///age: 64///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1875///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                         1 
   samplename: NCI_21588///dataset: NCI///series: NCI///id: 21588///er: 1///pgr: NA///her2: NA///size: 7///node: 1///age: 49///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 184///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                         1 
 samplename: NCI_21589///dataset: NCI///series: NCI///id: 21589///er: 0///pgr: NA///her2: NA///size: 4.5///node: 0///age: 73///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 252///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                         1 
 samplename: NCI_21591///dataset: NCI///series: NCI///id: 21591///er: 0///pgr: NA///her2: NA///size: 1.8///node: 1///age: 62///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 574///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21592///dataset: NCI///series: NCI///id: 21592///er: 0///pgr: NA///her2: NA///size: 1.7///node: 1///age: 33///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 3608///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                         1 
  samplename: NCI_21593///dataset: NCI///series: NCI///id: 21593///er: 1///pgr: NA///her2: NA///size: 4///node: 1///age: 51///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2717///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21594///dataset: NCI///series: NCI///id: 21594///er: 1///pgr: NA///her2: NA///size: 2.2///node: 1///age: 61///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2800///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
   samplename: NCI_21595///dataset: NCI///series: NCI///id: 21595///er: 1///pgr: NA///her2: NA///size: 6///node: 0///age: 53///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 455///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
  samplename: NCI_21596///dataset: NCI///series: NCI///id: 21596///er: 1///pgr: NA///her2: NA///size: 2///node: 1///age: 33///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 3639///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
  samplename: NCI_21597///dataset: NCI///series: NCI///id: 21597///er: 1///pgr: NA///her2: NA///size: 3///node: 0///age: 64///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 3491///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
   samplename: NCI_21598///dataset: NCI///series: NCI///id: 21598///er: 1///pgr: NA///her2: NA///size: 5///node: 1///age: 43///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 770///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
  samplename: NCI_21599///dataset: NCI///series: NCI///id: 21599///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 67///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 1382///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21600///dataset: NCI///series: NCI///id: 21600///er: 1///pgr: NA///her2: NA///size: 1.5///node: 0///age: 57///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2550///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                         1 
  samplename: NCI_21603///dataset: NCI///series: NCI///id: 21603///er: 1///pgr: NA///her2: NA///size: 4///node: 1///age: 35///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 3222///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                         1 
samplename: NCI_21604///dataset: NCI///series: NCI///id: 21604///er: 0///pgr: NA///her2: NA///size: 1.5///node: 1///age: 64///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1038///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                         1 
   samplename: NCI_21605///dataset: NCI///series: NCI///id: 21605///er: 0///pgr: NA///her2: NA///size: 2///node: 1///age: 55///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 317///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                         1 
 samplename: NCI_21606///dataset: NCI///series: NCI///id: 21606///er: 1///pgr: NA///her2: NA///size: 1.8///node: 0///age: 54///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 996///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                         1 
samplename: NCI_21607///dataset: NCI///series: NCI///id: 21607///er: 0///pgr: NA///her2: NA///size: 2.6///node: 0///age: 53///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 3449///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                         1 
samplename: NCI_21608///dataset: NCI///series: NCI///id: 21608///er: 1///pgr: NA///her2: NA///size: 1.3///node: 0///age: 65///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 3334///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                         1 
  samplename: NCI_21609///dataset: NCI///series: NCI///id: 21609///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 59///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 1476///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
  samplename: NCI_21610///dataset: NCI///series: NCI///id: 21610///er: 0///pgr: NA///her2: NA///size: 4///node: 0///age: 45///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 3489///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                         1 
   samplename: NCI_21611///dataset: NCI///series: NCI///id: 21611///er: 1///pgr: NA///her2: NA///size: 4///node: 1///age: 60///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 213///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
  samplename: NCI_21612///dataset: NCI///series: NCI///id: 21612///er: 1///pgr: NA///her2: NA///size: 6///node: 1///age: 53///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 3556///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
  samplename: NCI_21613///dataset: NCI///series: NCI///id: 21613///er: 0///pgr: NA///her2: NA///size: 3///node: 0///age: 43///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 3616///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                         1 
samplename: NCI_21614///dataset: NCI///series: NCI///id: 21614///er: 0///pgr: NA///her2: NA///size: 1.8///node: 0///age: 44///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 3004///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                         1 
  samplename: NCI_21615///dataset: NCI///series: NCI///id: 21615///er: 1///pgr: NA///her2: NA///size: 2///node: 0///age: 63///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 3267///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                         1 
  samplename: NCI_21616///dataset: NCI///series: NCI///id: 21616///er: 1///pgr: NA///her2: NA///size: 3///node: 0///age: 41///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 3527///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21617///dataset: NCI///series: NCI///id: 21617///er: 1///pgr: NA///her2: NA///size: 2.2///node: 0///age: 48///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 3391///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21618///dataset: NCI///series: NCI///id: 21618///er: 1///pgr: NA///her2: NA///size: 1.5///node: 1///age: 63///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 3220///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
 samplename: NCI_21619///dataset: NCI///series: NCI///id: 21619///er: 1///pgr: NA///her2: NA///size: 1.9///node: 0///age: 49///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 938///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
  samplename: NCI_21620///dataset: NCI///series: NCI///id: 21620///er: 1///pgr: NA///her2: NA///size: 2///node: 0///age: 61///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 3131///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21621///dataset: NCI///series: NCI///id: 21621///er: 0///pgr: NA///her2: NA///size: 1.5///node: 0///age: 65///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 3014///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
   samplename: NCI_21622///dataset: NCI///series: NCI///id: 21622///er: 0///pgr: NA///her2: NA///size: 3///node: 1///age: 52///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 405///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21623///dataset: NCI///series: NCI///id: 21623///er: 0///pgr: NA///her2: NA///size: 4.5///node: 1///age: 57///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2057///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21624///dataset: NCI///series: NCI///id: 21624///er: 0///pgr: NA///her2: NA///size: 2.5///node: 0///age: 54///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 3374///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21625///dataset: NCI///series: NCI///id: 21625///er: 1///pgr: NA///her2: NA///size: 1.3///node: 0///age: 48///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 3370///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21626///dataset: NCI///series: NCI///id: 21626///er: 1///pgr: NA///her2: NA///size: 0.8///node: 1///age: 83///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2690///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21627///dataset: NCI///series: NCI///id: 21627///er: 1///pgr: NA///her2: NA///size: 4.5///node: 1///age: 61///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 1438///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
 samplename: NCI_21628///dataset: NCI///series: NCI///id: 21628///er: 0///pgr: NA///her2: NA///size: 2.5///node: 1///age: 70///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 853///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21629///dataset: NCI///series: NCI///id: 21629///er: 0///pgr: NA///her2: NA///size: 2.5///node: 0///age: 54///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2106///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
  samplename: NCI_21630///dataset: NCI///series: NCI///id: 21630///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 70///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2809///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21631///dataset: NCI///series: NCI///id: 21631///er: 1///pgr: NA///her2: NA///size: 1.9///node: 0///age: 61///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 3125///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21632///dataset: NCI///series: NCI///id: 21632///er: 1///pgr: NA///her2: NA///size: 2.6///node: 1///age: 56///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2728///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
 samplename: NCI_21633///dataset: NCI///series: NCI///id: 21633///er: 0///pgr: NA///her2: NA///size: 3.6///node: 1///age: 54///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 879///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21634///dataset: NCI///series: NCI///id: 21634///er: 0///pgr: NA///her2: NA///size: 2.5///node: 0///age: 64///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 1805///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21635///dataset: NCI///series: NCI///id: 21635///er: 1///pgr: NA///her2: NA///size: 1.5///node: 1///age: 34///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 3098///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21636///dataset: NCI///series: NCI///id: 21636///er: 1///pgr: NA///her2: NA///size: 1.7///node: 0///age: 57///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2491///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21637///dataset: NCI///series: NCI///id: 21637///er: 1///pgr: NA///her2: NA///size: 2.8///node: 1///age: 73///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2205///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21638///dataset: NCI///series: NCI///id: 21638///er: 1///pgr: NA///her2: NA///size: 0.8///node: 0///age: 56///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2734///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21639///dataset: NCI///series: NCI///id: 21639///er: 0///pgr: NA///her2: NA///size: 2.4///node: 1///age: 40///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2230///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21640///dataset: NCI///series: NCI///id: 21640///er: 1///pgr: NA///her2: NA///size: 2.3///node: 1///age: 68///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1548///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
 samplename: NCI_21641///dataset: NCI///series: NCI///id: 21641///er: 1///pgr: NA///her2: NA///size: 2.5///node: 1///age: 72///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 274///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21642///dataset: NCI///series: NCI///id: 21642///er: 1///pgr: NA///her2: NA///size: 1.5///node: 0///age: 48///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2981///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
  samplename: NCI_21643///dataset: NCI///series: NCI///id: 21643///er: 0///pgr: NA///her2: NA///size: 1///node: 0///age: 48///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2107///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
   samplename: NCI_21644///dataset: NCI///series: NCI///id: 21644///er: 0///pgr: NA///her2: NA///size: 8///node: 1///age: 59///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 657///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21645///dataset: NCI///series: NCI///id: 21645///er: 1///pgr: NA///her2: NA///size: 1.8///node: 0///age: 77///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2560///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
 samplename: NCI_21647///dataset: NCI///series: NCI///id: 21647///er: 0///pgr: NA///her2: NA///size: 2.5///node: 1///age: 74///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 384///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
 samplename: NCI_21648///dataset: NCI///series: NCI///id: 21648///er: 1///pgr: NA///her2: NA///size: 2.8///node: 1///age: 60///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 132///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21649///dataset: NCI///series: NCI///id: 21649///er: 1///pgr: NA///her2: NA///size: 2.2///node: 0///age: 54///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2153///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21650///dataset: NCI///series: NCI///id: 21650///er: 1///pgr: NA///her2: NA///size: 2.9///node: 0///age: 59///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2880///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21652///dataset: NCI///series: NCI///id: 21652///er: 1///pgr: NA///her2: NA///size: 0.9///node: 0///age: 65///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 3038///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
 samplename: NCI_21653///dataset: NCI///series: NCI///id: 21653///er: 0///pgr: NA///her2: NA///size: 2.6///node: 0///age: 51///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 567///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
  samplename: NCI_21655///dataset: NCI///series: NCI///id: 21655///er: 0///pgr: NA///her2: NA///size: 2///node: 0///age: 47///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2142///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
  samplename: NCI_21656///dataset: NCI///series: NCI///id: 21656///er: 1///pgr: NA///her2: NA///size: 1///node: 0///age: 61///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2861///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
     samplename: NCI_21657///dataset: NCI///series: NCI///id: 21657///er: 1///pgr: NA///her2: NA///size: 6///node: 1///age: 73///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 8///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 0 
                                                                                                                                                                                                                                         1 
samplename: NCI_21658///dataset: NCI///series: NCI///id: 21658///er: 1///pgr: NA///her2: NA///size: 2.8///node: 1///age: 45///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 1926///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21659///dataset: NCI///series: NCI///id: 21659///er: 0///pgr: NA///her2: NA///size: 3.2///node: 0///age: 64///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 1367///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
 samplename: NCI_21660///dataset: NCI///series: NCI///id: 21660///er: 0///pgr: NA///her2: NA///size: 6.5///node: 0///age: 49///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 233///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                         1 
samplename: NCI_21661///dataset: NCI///series: NCI///id: 21661///er: 1///pgr: NA///her2: NA///size: 2.6///node: 1///age: 44///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1918///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21662///dataset: NCI///series: NCI///id: 21662///er: 1///pgr: NA///her2: NA///size: 1.6///node: 1///age: 61///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2072///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
   samplename: NCI_21663///dataset: NCI///series: NCI///id: 21663///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 50///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 861///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
   samplename: NCI_21664///dataset: NCI///series: NCI///id: 21664///er: 1///pgr: NA///her2: NA///size: 5///node: 1///age: 53///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 349///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
  samplename: NCI_21665///dataset: NCI///series: NCI///id: 21665///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 46///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2071///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
  samplename: NCI_21666///dataset: NCI///series: NCI///id: 21666///er: 1///pgr: NA///her2: NA///size: 3///node: 1///age: 56///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 1841///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21667///dataset: NCI///series: NCI///id: 21667///er: 1///pgr: NA///her2: NA///size: 2.5///node: 0///age: 67///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 1528///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
  samplename: NCI_21668///dataset: NCI///series: NCI///id: 21668///er: 1///pgr: NA///her2: NA///size: 5///node: 1///age: 48///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1469///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
  samplename: NCI_21669///dataset: NCI///series: NCI///id: 21669///er: 0///pgr: NA///her2: NA///size: 3///node: 1///age: 61///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1141///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
  samplename: NCI_21670///dataset: NCI///series: NCI///id: 21670///er: 1///pgr: NA///her2: NA///size: 3///node: 0///age: 71///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 1909///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21671///dataset: NCI///series: NCI///id: 21671///er: 0///pgr: NA///her2: NA///size: 2.6///node: 1///age: 50///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1979///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
  samplename: NCI_21672///dataset: NCI///series: NCI///id: 21672///er: 1///pgr: NA///her2: NA///size: 2///node: 0///age: 56///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 1209///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
 samplename: NCI_21674///dataset: NCI///series: NCI///id: 21674///er: 1///pgr: NA///her2: NA///size: 5.5///node: 0///age: 75///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 903///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
  samplename: NCI_21675///dataset: NCI///series: NCI///id: 21675///er: 0///pgr: NA///her2: NA///size: 3///node: 0///age: 73///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1990///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
  samplename: NCI_21676///dataset: NCI///series: NCI///id: 21676///er: 1///pgr: NA///her2: NA///size: 3///node: 0///age: 41///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2242///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21677///dataset: NCI///series: NCI///id: 21677///er: 1///pgr: NA///her2: NA///size: 1.8///node: 1///age: 57///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2262///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
  samplename: NCI_21678///dataset: NCI///series: NCI///id: 21678///er: 0///pgr: NA///her2: NA///size: 4///node: 1///age: 54///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1462///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
 samplename: NCI_21679///dataset: NCI///series: NCI///id: 21679///er: 0///pgr: NA///her2: NA///size: 3.8///node: 1///age: 90///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 540///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
  samplename: NCI_21680///dataset: NCI///series: NCI///id: 21680///er: 0///pgr: NA///her2: NA///size: 1///node: 0///age: 33///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 1197///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                         1 
  samplename: NCI_21681///dataset: NCI///series: NCI///id: 21681///er: 0///pgr: NA///her2: NA///size: 2///node: 1///age: 50///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1136///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                         1 
  samplename: NCI_21682///dataset: NCI///series: NCI///id: 21682///er: 1///pgr: NA///her2: NA///size: 2///node: 1///age: 70///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 2226///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
   samplename: NCI_21823///dataset: NCI///series: NCI///id: 21823///er: 1///pgr: NA///her2: NA///size: 4///node: 1///age: 82///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 905///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21824///dataset: NCI///series: NCI///id: 21824///er: 1///pgr: NA///her2: NA///size: 1.5///node: 1///age: 59///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 2431///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
  samplename: NCI_21825///dataset: NCI///series: NCI///id: 21825///er: 1///pgr: NA///her2: NA///size: 3///node: 0///age: 63///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 3515///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21826///dataset: NCI///series: NCI///id: 21826///er: 1///pgr: NA///her2: NA///size: 1.4///node: 0///age: 49///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 3226///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21827///dataset: NCI///series: NCI///id: 21827///er: 1///pgr: NA///her2: NA///size: 1.5///node: 0///age: 71///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 2768///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
 samplename: NCI_21828///dataset: NCI///series: NCI///id: 21828///er: 1///pgr: NA///her2: NA///size: 2.4///node: 0///age: 75///grade: 1///t.dmfs: NA///e.dmfs: NA///t.rfs: 623///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
samplename: NCI_21829///dataset: NCI///series: NCI///id: 21829///er: 1///pgr: NA///her2: NA///size: 1.5///node: 1///age: 38///grade: 2///t.dmfs: NA///e.dmfs: NA///t.rfs: 1918///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 1 
                                                                                                                                                                                                                                         1 
   samplename: NCI_21830///dataset: NCI///series: NCI///id: 21830///er: 1///pgr: NA///her2: NA///size: 4///node: 1///age: 75///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 597///e.rfs: 1///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 
  samplename: NCI_21831///dataset: NCI///series: NCI///id: 21831///er: 0///pgr: NA///her2: NA///size: 5///node: 1///age: 59///grade: 3///t.dmfs: NA///e.dmfs: NA///t.rfs: 1858///e.rfs: 0///t.os: NA///e.os: NA///tissue: 1///treatment: 2 
                                                                                                                                                                                                                                         1 

}}
\source{ Not available }
\keyword{datasets}
