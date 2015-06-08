# README #

The status of curation is summarized in a Google Spreadsheet [here](https://docs.google.com/spreadsheets/d/1R-n2huJ8xU1eM1s5pi5okZEBkJALsFRh5ji3TYDZZiQ/edit?usp=sharing).

Outstanding TODO items and known issues are recorded in the "Issues" section.

Please contact Levi Waldron (levi.waldron@hunter.cuny.edu) or Benjamin Haibe-Kains (benjamin.haibe.kains@utoronto.ca) with questions or comments.

## Downloading datasets from GEO:

Download the R script gse_PROCESSED.R from the curatedOvarianData pipeline [here](https://bitbucket.org/lwaldron/curatedovariandata/src/tip/src/).  Then you can run this R script from a terminal command line as follows, e.g. to get expression and clinical data for GSE12418:

```
#!r

R CMD BATCH --vanilla "--args GSE12418 ./tmp" gse_PROCESSED.R GSE12418_downloadPROCESSED.log
```

If the series has multiple platforms (GPL) and you only want one platform, invoke gse_PROCESSED.R like this:

```
#!r

R CMD BATCH --vanilla "--args GSE32062 ./tmp GPL6480" gse_PROCESSED.R GSE32062-GPL6480_downloadPROCESSED.log
```

Note that ./tmp will put results inside the tmp directory underneath the present working directory where this command is run.  This command assumes that gse_PROCESSED.R is in the directory where you run the above command.  Any of the filenames (output directory, location of gse_PROCESSED.R, and name+location of the logfile can be changed.


