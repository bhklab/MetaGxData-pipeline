## curation of METABRIC
## Benjamin Haibe-Kains

# EGAD00010000210   Illumina HT 12     997   Normalized expression data; discovery set
# EGAD00010000211   Illumina HT 12     995   Normalized expression data; validation set
# EGAD00010000212   Illumina HT 12     144   Normalized expression data; normals

.libPaths(c("/mnt/work1/users/bhklab/Rlib", .libPaths()))

## remove all existing objects from your workspace
rm(list=ls())

original.dir <- getwd()
setwd("/mnt/work1/users/bhklab/Data/METABRIC")

## create a directory to save all the results we generate
saveres <- "datasets"
if(!file.exists(saveres)) { dir.create(file.path(saveres), showWarnings=FALSE) }
  
## load required libraries
library(lumi)
library(genefu)

badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"


## read gene expression data
myfn <- file.path(saveres, "dd1.RData")
if(!file.exists(myfn)) {
  dd1 <- t(read.csv(file=file.path("EGAD00010000210", "discovery_ExpressionMatrix.txt"), sep="\t"))
  rownames(dd1) <- toupper(gsub(pattern=badchars, replacement="_", x=rownames(dd1)))
  save(list=c("dd1"), compress=TRUE, file=myfn)
} else { load(myfn) }

## weird separation by space " "
myfn <- file.path(saveres, "dd2.RData")
if(!file.exists(myfn)) {
  dd2 <- t(read.csv(file=file.path("EGAD00010000211", "validation_ExpressionMatrix.txt"), sep=" "))
  rownames(dd2) <- toupper(gsub(pattern=badchars, replacement="_", x=rownames(dd2)))
  save(list=c("dd2"), compress=TRUE, file=myfn)
} else { load(myfn) }

myfn <- file.path(saveres, "dd3.RData")
if(!file.exists(myfn)) {
  dd3 <- t(read.csv(file=file.path("EGAD00010000212", "normals_ExpressionMatrix.txt"), sep="\t"))
  metabricid.dd3 <- rownames(dd3)
  rownames(dd3) <- toupper(gsub(pattern=badchars, replacement="_", x=rownames(dd3)))
  save(list=c("dd3", "metabricid.dd3"), compress=TRUE, file=myfn)
} else { load(myfn) }

## merge gene expression data
cc <- intersect(intersect(colnames(dd1), colnames(dd2)), colnames(dd3))
tumor.ge <- rbind(dd1[ , cc, drop=FALSE], dd2[ , cc, drop=FALSE])
normal.ge <- dd3[ , cc, drop=FALSE]

## annotations
annot.ge <- lumi::IlluminaID2nuID(IlluminaID=cc, lib.mapping="lumiHumanIDMapping")
annot.ge <- cbind(annot.ge, "EntrezGene.ID"=lumi::nuID2EntrezID(nuID=as.character(annot.ge[ , "nuID"]), lib.mapping="lumiHumanIDMapping"), "EntrezID"=lumi::nuID2EntrezID(nuID=as.character(annot.ge[ , "nuID"]), lib.mapping="lumiHumanIDMapping"))
annot.ge[annot.ge == ""] <- NA

## read clinical information
ss1 <- read.csv(file.path("clinical_information", "table_S2_revised.txt"), stringsAsFactors=FALSE, sep="\t")
rownames(ss1) <- toupper(gsub(pattern=badchars, replacement="_", x=ss1[ ,1]))
ss2 <- read.csv(file.path("clinical_information", "table_S3_revised.txt"), stringsAsFactors=FALSE, sep="\t")
rownames(ss2) <- toupper(gsub(pattern=badchars, replacement="_", x=ss2[ ,1]))

## merge clinical information
cc <- intersect(colnames(ss1), colnames(ss2))
tumor.info <- rbind(ss1[ , cc, drop=FALSE], ss2[ , cc, drop=FALSE])
tumor.info <- cbind(tumor.info, "dataset2"=c(rep("discovery", nrow(ss1)), rep("validation", nrow(ss2))), "tissue2"="tumor")
ss3 <- matrix(NA, nrow=nrow(dd3), ncol=ncol(tumor.info), dimnames=list(rownames(dd3), colnames(tumor.info)))
ss3[ ,"tissue2"] <- "normal"
ss3[ , "METABRIC_ID"] <- metabricid.dd3
normal.info <- ss3

intclust.temp <- read.csv(file.path("clinical_information", "table_S44.csv"), stringsAsFactors=FALSE, row.names=1)
rownames(intclust.temp) <- toupper(gsub(pattern=badchars, replacement="_", x=rownames(intclust.temp)))
iix <- intersect(rownames(tumor.info), rownames(intclust.temp))
tumor.info[iix, "IntClustMemb"] <- intclust.temp[iix, "IntClustMemb"]

## add column with standard names
# cinfo <- c("samplename", "dataset", "series", "id", "age", "size", "node", "er", "pgr", "her2", "grade", "t.rfs", "e.rfs", "t.dmfs", "e.dmfs", "t.os", "e.os", "treatment", "tissue")
size <- tumor.info[ , "size"]
size[size == "null"] <- NA
size <- as.numeric(size) / 10
ers <- factor(tumor.info[ , "ER_IHC_status"], levels=c("neg", "pos"))
levels(ers) <- c(0, 1)
her2s <- factor(tumor.info[ , "HER2_IHC_status"], levels=c(0, 1, 2, 3))
levels(her2s) <- c(0, 0, 1, 1)
grade <- factor(tumor.info[ , "grade"], levels=c(1, 2, 3))
tt <- factor(tumor.info[ , "Treatment"], levels=c("NONE", "CT", "CT/HT", "CT/HT/RT", "CT/RT", "HT", "HT/RT", "RT"))
levels(tt) <- c(0, 1, 6, 6, 1, 2, 2, 0)
## tumor
tumor.info <- data.frame("samplename"=rownames(tumor.info), "id"=rownames(tumor.info), "dataset"="METABRIC", "series"=as.character(tumor.info[ , "Site"]), "node"=as.numeric(tumor.info[ , "lymph_nodes_positive"] > 0), "age"=as.numeric(tumor.info[ , "age_at_diagnosis"]), "size"=size, "er"=as.numeric(ers)-1, "her2"=as.numeric(her2s)-1, "pgr"=NA, "grade"=as.numeric(grade), "treatment"=tt, "tissue"=1, "brca.mutation"=NA, tumor.info)
## normal
normal.info <- data.frame("samplename"=rownames(normal.info), "id"=rownames(normal.info), "dataset"="METABRIC", "series"=as.character(normal.info[ , "Site"]), "node"=NA, "age"=NA, "size"=NA, "er"=NA, "her2"=NA, "pgr"=NA, "grade"=NA, "treatment"=0, "tissue"=0, "brca.mutation"=NA, normal.info)


## suvival data
## a, d, d-d.s, d-o.c  alive, dead, dead (disease-Specific) and dead (other causes)

## overall survival
t.os <- tumor.info[ , "T"]
e.os <- rep(NA, nrow(tumor.info))
names(e.os) <- names(t.os) <- rownames(tumor.info)
e.os[!is.na(tumor.info[ , "last_follow_up_status"]) & is.element(tumor.info[ , "last_follow_up_status"], c("d", "d-d.s.", "d-o.c."))] <- 1
e.os[!is.na(tumor.info[ , "last_follow_up_status"]) & is.element(tumor.info[ , "last_follow_up_status"], c("a"))] <- 0
## plot survival curves
mysurv <- survcomp::censor.time(surv.time=t.os / 30, surv.event=e.os, time.cens=15 * 12)
myx <- tumor.info[ , "dataset2"] == "discovery"
## discovery set
message("Discovery set (OS)")
sf <- survfit(Surv(mysurv[[1]], mysurv[[2]])[myx] ~ 1)
# plot(sf)
## median survival time
iix <- max(which(sf$surv > 0.5))
message(sprintf("median survival time (months): %.3g", sf$time[iix]))
## Integrated Clustering members
intclust <- factor(tumor.info[myx, "IntClustMemb"], levels=1:length(table(tumor.info[myx, "IntClustMemb"])), ordered=FALSE)
names(intclust) <- rownames(tumor.info)[myx]
intclust.n <- table(tumor.info[names(intclust), "IntClustMemb"])
intclust.nd <- sapply(levels(intclust), function(x, y, z) { return(sum(z[!is.na(y) & y == x], na.rm=TRUE)) }, y=intclust, z=mysurv[[2]][names(intclust)])
sf <- survfit(Surv(mysurv[[1]], mysurv[[2]])[myx] ~ intclust)
mycol <- rainbow(length(levels(intclust)))
lrp <- survdiff(Surv(mysurv[[1]], mysurv[[2]])[myx] ~ intclust)
lrp <- 1 - pchisq(lrp$chisq, df=length(lrp$n) - 1)
pdf(file.path("datasets", "intclust_os_discoveryset.pdf"))
plot(sf, col=mycol) 
legend("bottomleft", legend=c(sprintf("IntClust%s %i(%i)", levels(intclust), intclust.n, intclust.nd), sprintf("Logrank p-value = %.2E", lrp)), col=c(mycol, "white"), lty=1, bty="n")
dev.off()

## disease-specific survival
t.dfs <- tumor.info[ , "T"]
e.dfs <- rep(NA, nrow(tumor.info))
names(e.dfs) <- names(t.dfs) <- rownames(tumor.info)
e.dfs[!is.na(tumor.info[ , "last_follow_up_status"]) & is.element(tumor.info[ , "last_follow_up_status"], c("d-d.s."))] <- 1
e.dfs[!is.na(tumor.info[ , "last_follow_up_status"]) & is.element(tumor.info[ , "last_follow_up_status"], c("a", "d", "d-o.c."))] <- 0
## plot survival curves
mysurv <- survcomp::censor.time(surv.time=t.dfs / 30, surv.event=e.dfs, time.cens=300)
myx <- tumor.info[ , "dataset2"] == "discovery"
## discovery set
message("Discovery set (DFS)")
sf <- survfit(Surv(mysurv[[1]], mysurv[[2]])[myx] ~ 1)
plot(sf)
## median survival time
iix <- max(which(sf$surv > 0.5))
message(sprintf("median survival time (months): %.3g", sf$time[iix]))
## Integrated Clustering members
intclust <- factor(tumor.info[myx, "IntClustMemb"], levels=1:length(table(tumor.info[myx, "IntClustMemb"])), ordered=FALSE)
names(intclust) <- rownames(tumor.info)[myx]
intclust.n <- table(tumor.info[names(intclust), "IntClustMemb"])
intclust.nd <- sapply(levels(intclust), function(x, y, z) { return(sum(z[!is.na(y) & y == x], na.rm=TRUE)) }, y=intclust, z=e.dfs[names(intclust)])
sf <- survfit(Surv(mysurv[[1]], mysurv[[2]])[myx] ~ intclust)
lrp <- survdiff(Surv(mysurv[[1]], mysurv[[2]])[myx] ~ intclust)
lrp <- 1 - pchisq(lrp$chisq, df=length(lrp$n) - 1)
mycol <- rainbow(length(levels(intclust)))
pdf(file.path("datasets", "intclust_dfs_discoveryset.pdf"))
plot(sf, col=mycol) 
legend("bottomleft", legend=c(sprintf("IntClust%s %i(%i)", levels(intclust), intclust.n, intclust.nd), sprintf("Logrank p-value = %.2E", lrp)), col=c(mycol, "white"), lty=1, bty="n")
dev.off()

## -> Supplementary Figure S27d
## validation set
message("Validation set")
sf <- survfit(Surv(mysurv[[1]], mysurv[[2]])[!myx] ~ 1)
plot(sf)
## median survival time
iix <- max(which(sf$surv > 0.5))
message(sprintf("median survival time (months): %.3g", sf$time[iix]))

## -> Supplementary Figure S35c
## validation set
message("Validation set")
sf <- survfit(Surv(mysurv[[1]], mysurv[[2]])[!myx] ~ 1)
plot(sf)
## median survival time
iix <- max(which(sf$surv > 0.5))
message(sprintf("median survival time (months): %.3g", sf$time[iix]))

tumor.info <- data.frame(tumor.info, "t.os"=t.os[rownames(tumor.info)], "e.os"=e.os[rownames(tumor.info)], "t.dfs"=t.dfs[rownames(tumor.info)], "e.dfs"=e.dfs[rownames(tumor.info)], "t.rfs"=NA, "e.rfs"=NA, "t.dmfs"=NA, "e.dmfs"=NA)

## match normal-tumor pair
dd <- read.csv(file.path("README", "tumour_normal_pair.txt"), sep="\t", stringsAsFactors=FALSE)
dd[ ,1] <- toupper(gsub(badchars, "_", dd[ ,1]))
dd[ ,2] <- toupper(gsub(badchars, "_", dd[ ,2]))
## tumor
iix <- is.element(dd[ ,1], rownames(tumor.info))
match.pair <- rep(NA, nrow(tumor.info))
names(match.pair) <- rownames(tumor.info)
match.pair[dd[iix,1]] <- dd[iix,2]
tumor.info <- data.frame(tumor.info, "MatchNormal"=match.pair)
tumor.info[ , "MatchNormal"] <- as.character(tumor.info[ , "MatchNormal"])
## normal
iix <- is.element(dd[ ,2], rownames(normal.info))
match.pair <- rep(NA, nrow(normal.info))
names(match.pair) <- rownames(normal.info)
match.pair[dd[iix,2]] <- dd[iix,1]
normal.info <- data.frame(normal.info, "MatchTumor"=match.pair)
normal.info[!is.na(normal.info[ , "MatchTumor"]), "age"] <- tumor.info[normal.info[!is.na(normal.info[ , "MatchTumor"]), "MatchTumor"], "age"]

## reorder patients
tumor.ge <- tumor.ge[rownames(tumor.info), rownames(annot.ge), drop=FALSE]

## save dataset
message("Save full dataset ...")
save(list=c("tumor.ge", "normal.ge", "annot.ge", "tumor.info", "normal.info"), compress=TRUE, file=file.path(saveres, "METABRIC_GE_FULL.RData"))

## merge normal and tumor
gn <- colnames(tumor.ge)
sample.ge <- rbind(tumor.ge[ , gn, drop=FALSE], normal.ge[ , gn, drop=FALSE])
cn <- union(colnames(tumor.info), colnames(normal.info))
class.col <- sapply(cn, function (x, f1, f2) {
  if(is.element(x, names(f1)) && !is.element(x, names(f2))) { return(f1[x]) }
  if(!is.element(x, names(f1)) && is.element(x, names(f2))) { return(f2[x]) }
  if(!is.element(x, names(f1)) && !is.element(x, names(f2))) { return("character") }
  if(is.element(x, names(f1)) && is.element(x, names(f2))) {
    if(f1[x] != f2[x]) { warning(sprintf("different class for %s: 5S vs. %s", x, f1[x], f2[x])) }
    return(f1[x])
  }
}, f1=sapply(tumor.info, class), f2=sapply(normal.info, class))
class.lev <- sapply(cn, function (x, f1, f2) {
  if(is.element(x, names(f1)) && !is.element(x, names(f2))) { return(f1[[x]]) }
  if(!is.element(x, names(f1)) && is.element(x, names(f2))) { return(f2[[x]]) }
  if(!is.element(x, names(f1)) && !is.element(x, names(f2))) { return(NULL) }
  if(is.element(x, names(f1)) && is.element(x, names(f2))) {
    return(sort(unique(c(f1[[x]], f2[[x]]))))
  }
}, f1=sapply(tumor.info, levels), f2=sapply(normal.info, levels))

sample.info <- data.frame(matrix(NA, nrow=nrow(tumor.info) + nrow(normal.info), ncol=length(cn), dimnames=list(c(rownames(tumor.info), rownames(normal.info)), cn)))
sample.info <- setcolclass.df(df=sample.info, colclass=class.col, factor.levels=class.lev)
sample.info[rownames(tumor.info), colnames(tumor.info)] <- tumor.info
sample.info[rownames(normal.info), colnames(normal.info)] <- normal.info

save(list=c("sample.ge", "annot.ge", "sample.info"), compress=TRUE, file=file.path(saveres, "METABRIC_GE_MERGED_FULL.RData"))

## concatenate the adjacent normal samples to the tumor datasets

## select the most variant probe for each EntrezGene ID
gid <- as.character(annot.ge[ ,"EntrezID"])
names(gid) <- rownames(annot.ge)
ugid <- sort(unique(gid))
names(ugid) <- paste("geneid", ugid, sep=".")
rr <- genefu::geneid.map(geneid1=gid, data1=tumor.ge, geneid2=ugid)
## filter and rename probes
tumor.ge <- tumor.ge[ , names(rr$geneid1), drop=FALSE]
normal.ge <- normal.ge[ , names(rr$geneid1), drop=FALSE]
annot.ge <- annot.ge[names(rr$geneid1), , drop=FALSE]
colnames(tumor.ge) <- colnames(normal.ge) <- rownames(annot.ge) <- names(rr$geneid2)
message("Save ENTREZ dataset ...")
save(list=c("tumor.ge", "normal.ge", "annot.ge", "tumor.info", "normal.info"), compress=TRUE, file=file.path(saveres, "METABRIC_GE_ENTREZ.RData"))

setwd(original.dir)
write.csv(sample.info, file="../curation/breast/uncurated/METABRIC.csv")
