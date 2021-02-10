#####################################
#### Get public data from TCGA ######
#####################################

###1. Install stable release from CRAN###
#install.packages("UCSCXenaTools")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("GEOquery")
library("UCSCXenaTools")
library("GEOquery")

### 2.Select datasets from TGCA Xena browser####
#### To define which type of data download. 
# A. Select HostName (gdcHub or tcgaHub)
# B.Select Dataset
#data(XenaData)
#XenaData$XenaDatasets
XenaGenerate(subset = XenaHostNames=="gdcHub") %>% 
  XenaFilter(filterDatasets = "methylation|phenotype") %>% 
  XenaFilter(filterDatasets = "BRCA") -> info
options(timeout=1000) # avoid time out
XenaQuery(info) %>%
  XenaDownload() -> xe_download
#Prepare data into R for analysis.
dat = XenaPrepare(xe_download)
#Get data and filter them to work on Age clock code
#Filter in probes used in Age clock pipeline

#### METHYLATION 450K########
met450=dat$TCGA.BRCA.methylation450.tsv.gz
colnames(met450)[1] <- "ProbeID"
met450=met450[, -grep("-06A$", colnames(met450))]###Exclude methastasic samples
#write.csv(met450, file=gzfile("Methylation450_BRCA.csv.gz"), row.names = FALSE) ### All delta matrix
###Filter in probes used in Age clock analysis
Age_markers=read.table("Age_clock_probes.txt")
Age_markers.ma=as.matrix(Age_markers)
Age_clock.450=met450[met450$ProbeID%in%Age_markers.ma,]
dim(Age_clock.450)
#write.csv(Age_clock.450, file=gzfile("AC_Methylation450_BRCA.csv.gz"), row.names = FALSE)
#####Filter samples with normal adjacent-tumor tissue
samples.tab=read.table("Samples_to_Use_TCGA.txt") ### 
samples.tab=as.matrix(samples.tab)
names.use <- names(Age_clock.450)[(names(Age_clock.450)%in%samples.tab)]
Age_clock.450.subset <- Age_clock.450[, names.use]
dim(Age_clock.450.subset)
write.csv(Age_clock.450.subset, file=gzfile("AC_NT_Meth450_BRCA.csv.gz"), row.names = FALSE)

#### METHYLATION 27K########
met27=dat$TCGA.BRCA.methylation27.tsv.gz
colnames(met27)[1] <- "ProbeID"
met27=met27[, -grep("-06A$", colnames(met27))]###Exclude methastasic samples
Age_markers=read.table("Age_clock_probes.txt")
Age_markers.ma=as.matrix(Age_markers)
Age_clock.27=met27[met27$ProbeID%in%Age_markers.ma,]
dim(Age_clock.27)
#write.csv(Age_clock.27, file=gzfile("Methylation27_BRCA.csv.gz"), row.names = FALSE)

###Filter in probes used in Age clock analysis
Age_markers=read.table("Age_clock_probes.txt")
Age_markers.ma=as.matrix(Age_markers)
Age_clock.27=met27[met27$ProbeID%in%Age_markers.ma,]
dim(Age_clock.27)
#write.csv(Age_clock.27, file=gzfile("AC_Methylation27_BRCA.csv.gz"), row.names = FALSE)
#####Filter samples with normal adjacent-tumor tissue
samples.tab=read.table("Samples_to_Use_TCGA.txt") ### 
samples.tab=as.matrix(samples.tab)
names.use <- names(Age_clock.27)[(names(Age_clock.27)%in%samples.tab)]
Age_clock.27.subset <- Age_clock.27[, names.use]
dim(Age_clock.27.subset)
write.csv(Age_clock.27.subset, file=gzfile("AC_NT_Meth27_BRCA.csv.gz"), row.names = FALSE)

####Clinical information#####
pheno=dat$TCGA.BRCA.GDC_phenotype.tsv.gz

### 3.Select datasets from GEO omnibus####
write.csv(pheno, "Pheno_data_BRCA.csv", row.names = FALSE)

library(GEOquery)
library(data.table)
#####Download Series matrix
#rm(list=ls())
Series=c("GSE67919")
gse1 <- getGEO(Series,GSEMatrix=TRUE)
datMet1 = exprs (gse1[[1]])

Age_markers=read.table("Age_clock_probes.txt")
Age_markers.ma=as.matrix(Age_markers)
datMet1=as.data.frame(datMet1)
setDT(datMet1, keep.rownames = "ProbeID")
Age_clock.datMet=datMet1[datMet1$ProbeID%in%Age_markers.ma,]
dim(Age_clock.datMet)
write.csv(Age_clock.datMet, file=paste("AC.Work",Series,"series.csv"))
Pheno=pData(phenoData(gse1[[1]]))
write.csv(Pheno, file=paste(Series,"Pheno.csv"))

### Analysing from raw data processing I recommend minfi package
