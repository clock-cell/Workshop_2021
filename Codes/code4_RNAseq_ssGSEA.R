############# Immnue characterization ##########
###1. Install stable release from CRAN###
library("UCSCXenaTools")
### 2.Select datasets from TGCA Xena browser####
XenaGenerate(subset = XenaHostNames=="gdcHub") %>% 
  XenaFilter(filterDatasets = "htseq_counts") %>% 
  XenaFilter(filterDatasets = "BRCA") -> info
options(timeout=1000) # avoid time out
XenaQuery(info) %>%
  XenaDownload() -> xe_download
#Prepare data into R for analysis.
dat = XenaPrepare(xe_download)

###Filter in only mRNA genes
mRNA=read.table("Anotations_mRNA_Final.txt")
mRNA.ma=as.matrix(mRNA$V1)
mRNA.counts=dat[dat$Ensembl_ID%in%mRNA.ma,]
dim(mRNA.counts)
#####Filter samples with normal adjacent-tumor tissue
samples.tab=read.table("Samples_to_Use_TCGA_RNAseq.txt") ### 
ens="Ensembl_ID"
samples.tab=rbind(samples.tab,ens)
samples.tab=as.matrix(samples.tab)
names.use <- names(mRNA.counts)[(names(mRNA.counts)%in%samples.tab)]
mRNA.counts <- mRNA.counts[, names.use]
dim(mRNA.counts)
write.csv(mRNA.counts, "Counts_NT_BRCA.csv", row.names = FALSE)

###Normalized 
library(DESeq2)
###Raw count matrix input
cts <- read.csv("Counts_NT_BRCA.csv",header=TRUE)
###integers values
cts[,-1]<- lapply(lapply(cts[,-1],round),as.integer) 
#str(cts) 
cts1<- as.matrix(cts[,-1]) 
row.names(cts1)<- cts[,1] 
#Order count matrix and sample annotation
cts1 <- cts1[,order(colnames(cts1))]
Pheno=read.csv("SampleAnnotation_27_450_TN_BRCA.csv",header=TRUE)
Pheno=Pheno[order(Pheno$Sample),]
Pheno=as.data.frame(Pheno)
# PhenoData matrix
Pheno=DataFrame(Pheno)
Pheno$Class_IHC=as.factor(Pheno$Class_IHC)
dds <- DESeqDataSetFromMatrix(countData = cts1,
                              colData = Pheno,
                              design = ~ Class_IHC)
# filter the data to remove genes with few counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
nrow(dds)
# Run DESeq
dds$Class_IHC <- relevel(dds$Class_IHC, ref = "LumA")
dds <- DESeq(dds)
res <- results(dds)
####normalized data
norm=counts(dds, normalized=T)
norm=as.data.frame(norm)
norm
library(dplyr)
library(tibble)
norm <- norm %>% 
  rownames_to_column(var = "Ensembl_ID")
Annotation=read.table("Anotations_mRNA_Final.txt", header=TRUE)
Annotation=Annotation[,1:2]
Norm.data <- merge(norm, Annotation, by.x = "Ensembl_ID", by.y = "Ensembl_ID" )
max=ncol(Norm.data)
max2=max-1
Norm.data= Norm.data[,c(1,max,2:max2)]
Norm.data=subset(Norm.data, select = -c(Ensembl_ID))
write.csv(Norm.data,"Normalized_T_BRCA_matched.csv", row.names = FALSE)


####Eliminate duplicates
delta=read.csv("Normalized_T_BRCA_matched.csv", header=T)
dim(delta)
maxCol=ncol(delta)
delta1=transform(delta, IQR=apply(delta[2:maxCol],1, IQR, na.rm = TRUE ))
gene=delta1[,1] 
IQR_delta=delta1$IQR
IQR_delta=as.numeric(as.character(IQR_delta))
mydat.max=do.call(rbind,lapply(split(delta1,gene),function(chunk) chunk[which.max(chunk$IQR),]))
##Delete columns added by IQR column
mydata=subset(mydat.max, select = -c(IQR))
dim(mydata)
library(dplyr)
library(tibble)
mydata <-  mydata%>% 
  rownames_to_column(var = "Gene")
mydata=subset(mydata, select = -c(Symbol))
write.csv(mydata, "Normalized_T_BRCA_matched_SD.csv", row.names = FALSE)


####Adjust to z.score to compere with other data from different source like microarrays
tabla=read.csv("Normalized_T_BRCA_matched_SD.csv", header=TRUE, row.names=1)
tabla=as.matrix(tabla)
z.score=scale(tabla,scale=TRUE)
write.csv(z.score,"Normalized_T_BRCA_matched_SD_z.score.csv")


library(GSEABase)
library(GSVAdata)
library(parallel)
library(GSVA)
expr=read.csv("Normalized_T_BRCA_matched_SD_z.score.csv", header=T, row.names = 1)
expr=as.matrix(expr)
names(dimnames(expr)) <- c("Gene", "")
expr
gsm <- getGmt("IMMUNE.gmt")
gsva=gsva(expr, gsm,
          method=c("ssgsea"),
          kcdf=c("Gaussian"), ##"Poisson" for RPKM or FPKM or transformed RNA-seq data back into integer counts
          #Gaussian for microarrays and RNA-seq normalized continuous values (voom,DESeq or EDGE)
          abs.ranking=FALSE,####
          min.sz=10, #!!!!!!!!!!!!Default 10 dependerá de cuantos genes hay por geneSet
          max.sz=Inf,
          parallel.sz=0,
          #parallel.type="SOCK",
          mx.diff=TRUE,
          tau=0.25,
          ssgsea.norm=TRUE, #####normalizing the scores by the absolute difference
          #between the minimum and the maximum
          verbose=TRUE)
write.csv(gsva,"ssGSEA_Immune_BRCA_NT.csv")
rm(list=ls())

