### module load R/3.3.0
.libPaths("/mnt/work1/users/bhklab/Rlib")

### install CRAN and Bioconductor packages
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("devtools", "gdata", "knitr", "HiDimDA", "survival", "reshape2", "genefu", "annotate", "hgu133plus2.db", "stringr", "survMisc", "xtable", "gridExtra", "Biobase", "GSVA", "sparsediscrim", "survcomp", "NMF", "ggplot2", "e1071", "randomForest", "clue", "GEOquery", "logging", "metafor"))

### build and install MetaGx
# library(devtools)
# devtools::install_github(repo="bhklab/MetaGx")
library(MetaGx)

### Load compendium of ovarian cancer datasets
# install.packages("/mnt/work1/users/bhklab/Data/MetaGxData-private/MetaGxOvarian_0.99.0.tar.gz", repos=NULL)
library(MetaGxBreast)

library(clue)
library(genefu)
library(survcomp)
library(Hmisc)

########################
### Functions
########################

### from Ravi Varadhan rvaradhan at jhmi.edu: https://stat.ethz.ch/pipermail/r-help/2010-April/236664.html
pMatrix.min <- function(A, B) {
# finds the permutation P of A such that ||PA - B|| is minimum in Frobenius norm 
# Uses the linear-sum assignment problem (LSAP) solver in the "clue" package
# Returns P%*%A and the permutation vector `pvec' such that 
# A[pvec, ] is the permutation of A closest to B
	n <- nrow(A)
	D <- matrix(NA, n, n)
	for (i in 1:n) {
	  for (j in 1:n) {
	    D[j, i] <- (sum((B[j, ] - A[i, ])^2))
	  }
  }
  vec <- c(clue::solve_LSAP(D))
  return (list(A=A[vec,], pvec=vec))
}

### add Surv objects as phenoData label "y" to the esets
demo.survdata <- function(df, survType=c("os", "rfs", "dmfs"), time.cens) {
 # survType <- match.arg(survType)
 st <- se <- NULL
 for (i in 1:length(survType)) {
   switch (survType[i],
     "os" = {
       if (is.null(st) || is.null(se)) {
         st <- as.numeric(df$t.os)
         se <- df$e.os == "1"
       } else {
         ccix <- complete.cases(st, se)
         st[!ccix] <- as.numeric(df$t.os[!ccix])
         se[!ccix] <- df$e.os[!ccix] == "1"
       }
      
     },
     "rfs" = {
       if (is.null(st) || is.null(se)) {
         st <- as.numeric(df$t.rfs)
         se <- df$e.rfs == "1"
       } else {
         ccix <- complete.cases(st, se)
         st[!ccix] <- as.numeric(df$t.rfs[!ccix])
         se[!ccix] <- df$e.rfs[!ccix] == "1"
       }
     },
     "dmfs" = {
       if (is.null(st) || is.null(se)) {
         st <- as.numeric(df$t.dmfs)
         se <- df$e.dmfs == "1"
       } else {
         ccix <- complete.cases(st, se)
         st[!ccix] <- as.numeric(df$t.dmfs[!ccix])
         se[!ccix] <- df$e.dmfs[!ccix] == "1"
       }
   })
 }
 ### time censoring
 cc.ix <- complete.cases(st, se)
 if(time.cens != 0) { 
 	st[cc.ix][st[cc.ix] > time.cens] <- time.cens
 	se[cc.ix][st[cc.ix] > time.cens] <- 0
 }
  names(st) <- names(se) <- rownames(df)
  return (Surv(st, se))
}
 
### add Surv objects as phenoData label "y" to the esets
eset.survdata <- function(eset, survType=c("os", "rfs", "dmfs"), time.cens=3650) {
  # survType <- match.arg(survType)
  st <- se <- NULL
  for (i in 1:length(survType)) {
    switch (survType[i],
      "os" = {
        if (is.null(st) || is.null(se)) {
          st <- as.numeric(pData(eset)$days_to_death)
          se <- pData(eset)$vital_status == "deceased"
        } else {
          ccix <- complete.cases(st, se)
          st[!ccix] <- as.numeric(pData(eset)$days_to_death[!ccix])
          se[!ccix] <- pData(eset)$vital_status[!ccix] == "deceased"
        }
        
      },
      "rfs" = {
        if (is.null(st) || is.null(se)) {
          st <- as.numeric(pData(eset)$days_to_tumor_recurrence)
          se <- pData(eset)$recurrence_status == "recurrence"
        } else {
          ccix <- complete.cases(st, se)
          st[!ccix] <- as.numeric(pData(eset)$days_to_tumor_recurrence[!ccix])
          se[!ccix] <- pData(eset)$recurrence_status[!ccix] == "recurrence"
        }
      },
      "dmfs" = {
        if (is.null(st) || is.null(se)) {
          st <- as.numeric(pData(eset)$dmfs_days)
          se <- pData(eset)$dmfs_status == "recurrence"
        } else {
          ccix <- complete.cases(st, se)
          st[!ccix] <- as.numeric(pData(eset)$dmfs_days[!ccix])
          se[!ccix] <- pData(eset)$dmfs_status[!ccix] == "recurrence"
        }
    })
  }
  ### time censoring
  cc.ix <- complete.cases(st, se)
  if(time.cens != 0) { 
  	st[cc.ix][st[cc.ix] > time.cens] <- time.cens
  	se[cc.ix][st[cc.ix] > time.cens] <- 0
  }
  names(st) <- names(se) <- rownames(pData(eset))
  return (Surv(st, se))
}

########################

### original microarray data
old.data.path <- "/mnt/work1/users/bhklab/Data/microarray_data/BREAST"

### check clinical data
checkClinical <- function (data.list, data.name) {
  data <- mydatasets$Original$data
  annot <- mydatasets$Original$annot
  demo <- mydatasets$Original$demo
  eset <- mydatasets$MetaGx
  
  message("Check Clinical data")
  
  ii <- intersect(rownames(demo), rownames(pData(eset)))
  if((length(ii) != ncol(eset)) || (length(ii) != nrow(demo))) {
    warning(sprintf("different sample names or number of samples in %s", data.name))
  }
  ### reorder samples
  demo <- demo[ii, , drop=FALSE]
  data <- data[ii, , drop=FALSE]
  eset <- eset[ , ii]

  ########################
  ### check er, pgr, her2, tumour_size, age, rfs, dmfs and os
  ########################
  clinvar <- cbind(c("er", "pgr", "her2", "age", "size", "grade", "t.os", "e.os", "t.rfs", "e.rfs", "t.dmfs", "e.dmfs"),
    c("er", "pgr", "her2", "age_at_initial_pathologic_diagnosis", "tumor_size", "grade", "days_to_death", "vital_status", "days_to_tumor_recurrence", "recurrence_status", "dmfs_days", "dmfs_status")
    )
  colnames(clinvar) <- c("original", "MetaGx")
  for (i in 1:nrow(clinvar)) {
    cc1 <- demo[ , clinvar[i, "original"]]
    cc1[is.na(cc1)] <- "NA"
    cc2 <- pData(eset)[ , clinvar[i, "MetaGx"]]
    cc2[is.na(cc2)] <- "NA"
    tt <- table(cc1, cc2) 
    if (sum(diag(tt)) != length(ii)) {
      warning(sprintf("metadata mismatch for %s/%s in %s", clinvar[i, "original"], clinvar[i, "MetaGx"], data.name))
    }
  }
}

########################
### check subtypes
########################
checkSubtyping <- function (data.list, data.name) {
  data <- mydatasets$Original$data
  annot <- mydatasets$Original$annot
  demo <- mydatasets$Original$demo
  eset <- mydatasets$MetaGx
  
  ii <- intersect(rownames(demo), rownames(pData(eset)))
  demo <- demo[ii, , drop=FALSE]
  data <- data[ii, , drop=FALSE]
  eset <- eset[ , ii]
  
  message("Check subtyping")
  ### pam50
  message("\tPAM50")
  sbt <- genefu::molecular.subtyping(sbt.model="pam50", data=data, annot=annot, do.mapping=TRUE)$subtype
  eset2 <- MetaGx::subtypeClassification(eset=eset, model="pam50")
  sbt2 <- MetaGx::getSubtype(eset2, method="class")
  tt <- table(sbt, sbt2)
  tt <- pMatrix.min(tt, diag(1, nrow(tt)))$A
  if (sum(tt[upper.tri(tt)]) + sum(sum(tt[lower.tri(tt)])) != 0) {
    warning(sprintf("Different PAM50 subtypes in %s", data.name))
  }
  ### scmgene
  message("\tSCMGENE")
  sbt <- genefu::molecular.subtyping(sbt.model="scmgene", data=data, annot=annot, do.mapping=TRUE)$subtype
  eset2 <- MetaGx::subtypeClassification(eset=eset, model="scmgene")
  sbt2 <- MetaGx::getSubtype(eset2, method="class")
  tt <- table(sbt, sbt2)
  tt <- pMatrix.min(tt, diag(1, nrow(tt)))$A
  if (sum(tt[upper.tri(tt)]) + sum(sum(tt[lower.tri(tt)])) != 0) {
    warning(sprintf("Different SCMGENE subtypes in %s", data.name))
  }
  ### scmod2
  message("\tSMOD2")
  sbt <- genefu::molecular.subtyping(sbt.model="scmod2", data=data, annot=annot, do.mapping=TRUE)$subtype
  eset2 <- MetaGx::subtypeClassification(eset=eset, model="scmod2")
  sbt2 <- MetaGx::getSubtype(eset2, method="class")
  tt <- table(sbt, sbt2)
  tt <- pMatrix.min(tt, diag(1, nrow(tt)))$A
  if (sum(tt[upper.tri(tt)]) + sum(sum(tt[lower.tri(tt)])) != 0) {
    warning(sprintf("Different SCMOD2 subtypes in %s", data.name))
  }
}

########################
### check prognostic value 
########################
checkPrognosis <- function (data.list, data.name, survType=c("rfs", "dmfs", "os"), time.cens=3650) {
  data <- mydatasets$Original$data
  annot <- mydatasets$Original$annot
  demo <- mydatasets$Original$demo
  eset <- mydatasets$MetaGx
  
  ii <- intersect(rownames(demo), rownames(pData(eset)))
  demo <- demo[ii, , drop=FALSE]
  data <- data[ii, , drop=FALSE]
  eset <- eset[ , ii]
  
  message("Check prognostic value")
  
  ### survival data
  ### original dataset
  ss <- demo.survdata(df=demo, survType=survType, time.cens=10 * 365)
  ss2 <- eset.survdata(eset=eset, survType=survType, time.cens=10 * 365)
  
  ### list of genes/signatures of interest
  # goi <- genefu::mod2$AURKA
  goi <- list("AURKA"=cbind("probe"="AURKA", "EntrezGene.ID"="6790", "coefficient"="1"),
    "STAT1"=cbind("probe"="STAT1", "EntrezGene.ID"="6772", "coefficient"="1"),
    "PLAU"=cbind("probe"="PLAU", "EntrezGene.ID"="5328", "coefficient"="1")
    )
  ### AURKA should be highly prognostic in the global population and ER+/HER2-
  ### STAT1 should be prognostic in the ER-/HER2-
  ### PLAU should be prognostic in the HER2+

  sbt <- genefu::molecular.subtyping(sbt.model="scmod2", data=data, annot=annot, do.mapping=TRUE)$subtype
  sbtl <- list("Global"=c("ER-/HER2-", "ER+/HER2- High Prolif", "ER+/HER2- Low Prolif", "HER2+"),
    "Luminals"=c("ER+/HER2- High Prolif", "ER+/HER2- Low Prolif"),
    "LuminalB"= "ER+/HER2- High Prolif",
    "LuminalA"= "ER+/HER2- Low Prolif",
    "Basal-like"="ER-/HER2-",
    "Her2"="HER2+"
    )

  eset2 <- MetaGx::subtypeClassification(eset=eset, model="scmod2")
  sbt2 <- MetaGx::getSubtype(eset2, method="class")
  sbtl2 <-  list("Global"=c("Basal", "LumB", "LumA", "Her2"),
    "Luminals"=c("LumB", "LumA"),
    "LuminalB"= "LumB",
    "LuminalA"= "LumA",
    "Basal-like"="Basal",
    "Her2"="Her2"
    )

  for (i in 1:length(goi)) {
    message(sprintf("\t%s", names(goi)[i]))
    ### score
    xx <- genefu::sig.score(x=goi[[i]], data=data, annot=annot, do.mapping=TRUE)
    xx2 <- genefu::sig.score(x=goi[[i]], data=t(exprs(eset)), annot=fData(eset), do.mapping=TRUE)
    pdf(sprintf("km_surv_curve_%s_%s.pdf", data.name, names(goi)[i]), width=6, height=21)
    par(mfrow=c(7,2))
    ### scatterplot of scores
    plot(x=xx$score, y=xx2$score, xlab="Original", ylab="MetaGx")
    cort <- cor.test(x=xx$score, y=xx2$score, method="pearson", use="complete.obs")
    legend("topleft", legend=sprintf("Pearson cor = %.3g, p = %.1E", cort$estimate, cort$p.value), bty="n")
    plot.new()
    ### survival curves
    for (j in 1:length(sbtl)) {
      sbtix <- is.element(sbt, sbtl[[j]])
      sbtix2 <- is.element(sbt2, sbtl2[[j]])
      ### D index
      dindex <- survcomp::D.index(x=xx$score[sbtix], surv.time=ss[sbtix, 1], surv.event=ss[sbtix, 2], na.rm=TRUE)
      dindex2 <- survcomp::D.index(x=xx$score[sbtix2], surv.time=ss2[sbtix2, 1], surv.event=ss2[sbtix2, 2], na.rm=TRUE)
      ### survival curves
      xxbin <- Hmisc::cut2(xx$score, g=2)
      levels(xxbin) <- c("Low", "High")
      km.coxph.plot(Surv(st, se) ~ xx, data.s=data.frame("xx"=xxbin, "st"=ss[ , 1] / 365, "se"=ss[ , 2])[sbtix, , drop=FALSE], x.label="Time (years)", y.label="Disease-free survival", main.title=sprintf("%s in %s", names(goi)[i], names(sbtl)[j]))
      legend("topright", legend=sprintf("D index = %.3g, p = %.1E", dindex[[1]], dindex[[6]]), bty="n")
      xxbin2 <- Hmisc::cut2(xx2$score, g=2)
      levels(xxbin2) <- c("Low", "High")
      km.coxph.plot(Surv(st, se) ~ xx, data.s=data.frame("xx"=xxbin2, "st"=ss2[ , 1] / 365, "se"=ss[ , 2])[sbtix2, , drop=FALSE], x.label="Time (years)", y.label="Disease-free survival", main.title=sprintf("%s in %s", names(goi)[i], names(sbtl2)[j]))
      legend("topright", legend=sprintf("D index = %.3g, p = %.1E", dindex2[[1]], dindex2[[6]]), bty="n")
    }
    dev.off()
  }
}


########################
### VDX
########################
message("----------\nVDX\n----------")
### original dataset
dd <- "minn2007"
load(file.path(old.data.path, dd, paste(dd, "RData", sep=".")))

### comparison with MetaGxBreast
dd2 <- "VDX"
eset <- get(dd2)

mydatasets <- list("Original"=list("data"=data, "annot"=annot, "demo"=demo), "MetaGx"=eset)

checkClinical (mydatasets, "VDX")
checkSubtyping(mydatasets, "VDX")
checkPrognosis(mydatasets, "VDX", survType=c("rfs", "dmfs"), time.cens=10 * 365)

rm(list=c("eset", "data", "annot", "demo"))
gc(T)

########################
### NKI
########################
message("----------\nNKI\n----------")
### original dataset
dd <- "nki2002"
load(file.path(old.data.path, dd, paste(dd, "RData", sep=".")))

### comparison with MetaGxBreast
dd2 <- "NKI"
eset <- get(dd2)

mydatasets <- list("Original"=list("data"=data, "annot"=annot, "demo"=demo), "MetaGx"=eset)

checkClinical (mydatasets, "NKI")
checkSubtyping(mydatasets, "NKI")
checkPrognosis(mydatasets, "NKI", survType=c("rfs", "dmfs"), time.cens=10 * 365)

rm(list=c("eset", "data", "annot", "demo"))
gc(T)

########################
### MAINZ
########################
message("----------\nMAINZ\n----------")
### original dataset
dd <- "schmidt2008"
load(file.path(old.data.path, dd, paste(dd, "RData", sep=".")))

### comparison with MetaGxBreast
dd2 <- "MAINZ"
eset <- get(dd2)

mydatasets <- list("Original"=list("data"=data, "annot"=annot, "demo"=demo), "MetaGx"=eset)

checkClinical (mydatasets, "MAINZ")
checkSubtyping(mydatasets, "MAINZ")
checkPrognosis(mydatasets, "MAINZ", survType=c("rfs", "dmfs"), time.cens=10 * 365)

rm(list=c("eset", "data", "annot", "demo"))
gc(T)

########################
### METABRIC
########################
message("----------\nMETABRIC\n----------")
### original dataset
dd <- "METABRIC"
load(file.path(old.data.path, dd, "METABRIC_GE_ENTREZ.RData"))
data <- tumor.ge
annot <- annot.ge
annot <- cbind(annot, "probe"=rownames(annot), "EntrezGene.ID"=annot[ , "EntrezID"])
demo <- tumor.info

### comparison with MetaGxBreast
dd2 <- "METABRIC"
eset <- get(dd2)

mydatasets <- list("Original"=list("data"=data, "annot"=annot, "demo"=demo), "MetaGx"=eset)

checkClinical (mydatasets, "METABRIC")
checkSubtyping(mydatasets, "METABRIC")
checkPrognosis(mydatasets, "METABRIC", survType=c("os"), time.cens=10 * 365)

rm(list=c("eset", "data", "annot", "demo"))
gc(T)








### end

