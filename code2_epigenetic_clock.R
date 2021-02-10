############################
#####Age clock pipiline#####
############################

###Code 2###
##1. Load library
library(WGCNA)
#library(sqldf)
library(impute)
library(RPMM)

#2. Normalization step
#BMIQ function from Teschendorff 2013 (Bioinformatics. 2013 Jan 15;29(2):189-96) 
# adjusts for the type-2 bias in Illumina Infinium 450k data by Michael Hovert
## All files must be on the same folder or indicate their location
source("NORMALIZATION.R")

###3. Age transformation and probe annotation functions
## All files must be on the same folder or indicate their location
trafo= function(x,adult.age=20) { x=(x+1)/(1+adult.age); y=ifelse(x<=1, log( x),x-1);y }
anti.trafo= function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }
probeAnnotation21kdatMethUsed=read.csv("probeAnnotation21kdatMethUsed.csv")
probeAnnotation27k=read.csv("datMiniAnnotation27k.csv")
datClock=read.csv("AdditionalFile3.csv")

# 4.Read in the DNA methylation data (beta values)
##To imporve computatinal time Beta matrix filtered for
# probes annotated datMiniAnnotation27k.csv and samples
dat0=read.csv("AC_NT_Meth450_BRCA.csv.gz") 
nSamples=dim(dat0)[[2]]-1
nProbes= dim(dat0)[[1]]

#########################################
#####Run step 5,6, 7 and 8 together!!!!#####
#########################################

#5. Create a log file which will be output into your directory
# The code serves to create a log file (for error checks etc).
# It will automatically create a log file.
file.remove("LogFile.txt")
file.create("LogFile.txt")
DoNotProceed=FALSE
cat(paste( "The methylation data set contains", nSamples, "samples (e.g. arrays) and ", nProbes, " probes."),file="LogFile.txt")
if (nSamples==0) {DoNotProceed=TRUE; cat(paste( "\n ERROR: There must be a data input error since there seem to be no samples.\n Make sure that you input a comma delimited file (.csv file)\n that can be read using the R command read.csv.sql . Samples correspond to columns in that file  ."), file="LogFile.txt",append=TRUE) } 
if (nProbes==0) {DoNotProceed=TRUE; cat(paste( "\n ERROR: There must be a data input error since there seem to be zero probes.\n Make sure that you input a comma delimited file (.csv file)\n that can be read using the R command read.csv.sql  CpGs correspond to rows.")   , file="LogFile.txt",append=TRUE) } 
if (  nSamples > nProbes  ) { cat(paste( "\n MAJOR WARNING: It worries me a lot that there are more samples than CpG probes.\n Make sure that probes correspond to rows and samples to columns.\n I wonder whether you want to first transpose the data and then resubmit them? In any event, I will proceed with the analysis."),file="LogFile.txt",append=TRUE) }
if (  is.numeric(dat0[,1]) ) { DoNotProceed=TRUE; cat(paste( "\n Error: The first column does not seem to contain probe identifiers (cg numbers from Illumina) since these entries are numeric values. Make sure that the first column of the file contains probe identifiers such as cg00000292. Instead it contains ", dat0[1:3,1]  ),file="LogFile.txt",append=TRUE)  } 
if (  !is.character(dat0[,1]) ) {  cat(paste( "\n Major Warning: The first column does not seem to contain probe identifiers (cg numbers from Illumina) since these entries are numeric values. Make sure that the first column of the file contains CpG probe identifiers such as cg00000292. Instead it contains ", dat0[1:3,1]  ),file="LogFile.txt",append=TRUE)  } 
datout=data.frame(Error=c("Input error. Please check the log file for details","Please read the instructions carefully."), Comment=c("", "email Steve Horvath."))
if ( ! DoNotProceed ) {
  nonNumericColumn=rep(FALSE, dim(dat0)[[2]]-1)
  for (i in 2:dim(dat0)[[2]] ){ nonNumericColumn[i-1]=! is.numeric(dat0[,i]) }
  if (  sum(nonNumericColumn) >0 ) { cat(paste( "\n MAJOR WARNING: Possible input error. The following samples contain non-numeric beta values: ", colnames(dat0)[-1][ nonNumericColumn], "\n Hint: Maybe you use the wrong symbols for missing data. Make sure to code missing values as NA in the Excel file. To proceed, I will force the entries into numeric values but make sure this makes sense.\n" ),file="LogFile.txt",append=TRUE)  } 
  XchromosomalCpGs=as.character(probeAnnotation27k$Name[probeAnnotation27k$Chr=="X"])
  selectXchromosome=is.element(dat0[,1], XchromosomalCpGs )
  selectXchromosome[is.na(selectXchromosome)]=FALSE
  meanXchromosome=rep(NA, dim(dat0)[[2]]-1)
  if (   sum(selectXchromosome) >=500 )  {
    meanXchromosome= as.numeric(apply( as.matrix(dat0[selectXchromosome,-1]),2,mean,na.rm=TRUE)) }
  if (  sum(is.na(meanXchromosome)) >0 ) { cat(paste( "\n \n Comment: There are lots of missing values for X chromosomal probes for some of the samples. This is not a problem when it comes to estimating age but I cannot predict the gender of these samples.\n " ),file="LogFile.txt",append=TRUE)  } 
  ## filter by match probes
  match1=match(probeAnnotation21kdatMethUsed$Name , dat0[,1])
  if  ( sum( is.na(match1))>0 ) { 
    missingProbes= probeAnnotation21kdatMethUsed$Name[!is.element( probeAnnotation21kdatMethUsed$Name , dat0[,1])]    
    DoNotProceed=TRUE; cat(paste( "\n \n Input error: You forgot to include the following ", length(missingProbes), " CpG probes (or probe names):\n ", 
                                  paste( missingProbes, sep="",collapse=", ")),file="LogFile.txt",append=TRUE)  }  
  
  ####6. Filter STEP: Restrict the data to 21k probes and ensure they are numeric
  match1=match(probeAnnotation21kdatMethUsed$Name , dat0[,1])
  if  ( sum( is.na(match1))>0 ) stop(paste(sum( is.na(match1)), "CpG probes cannot be matched"))
  dat1= dat0[match1,]
  asnumeric1=function(x) {as.numeric(as.character(x))}
  dat1[,-1]=apply(as.matrix(dat1[,-1]),2,asnumeric1)
  
  #7. Create the output file called datout
  set.seed(1)
  # Do you want to normalize the data (recommended)? 
  #Just run normalizeData=FALSE in order to verify R code is not crashing, 
  #otherwise the results are not going to be realiable. 
  normalizeData=TRUE
  source("StepwiseAnalysis.txt")
  # 8. Output the results 
  if (  sum(  datout$Comment  != "" )   ==0 ) { cat(paste( "\n The individual samples appear to be fine. "),file="LogFile.txt",append=TRUE)  } 
  if (  sum(  datout$Comment != "" )   >0 ) { cat(paste( "\n Warnings were generated for the following samples.\n", datout[,1][datout$Comment != ""], "\n Hint: Check the output file for more details."),file="LogFile.txt",append=TRUE)  } 
} 

#####End run here!!!!###
datout1=datout
#9. Relate DNAm age to chronological age
#To address this task, we read in the sample annotation data that contain the chronological ages.
#Make sure it includes columns: Tissue, Age, Female (Female=1, Male=0). Note the spelling (capitalization).
datSample=read.csv("SampleAnnotation_450_NT_BRCA.csv")
#Order samples in the same way in the output file and sample annotation file
datSample <- datSample[order(datSample$Sample),]
datout1 <- datout1[order(datout1$SampleID),]
DNAmAge=datout1$DNAmAge
medianAbsDev=function(x,y) median(abs(x-y),na.rm=TRUE)
medianAbsDev1=signif(medianAbsDev(DNAmAge, datSample$Age),2)
#par(mfrow=c(1,1))
#verboseScatterplot(DNAmAge, datSample$Age,xlab="DNAm Age", ylab="Chronological Age",main=paste("All, err=", medianAbsDev1) );abline(0,1) 

#10. Compute DNAm age acceleration in samples
# FIRST ACCELERATION SCORE: The first acceleration meausure is based on the difference.
AgeAccelerationDiff=DNAmAge- datSample$Age
# SECOND ACCELERATION SCORE: residual resulting from regressing DNAmAge on chronological age
Age=datSample$Age
restNonMissing= !is.na(DNAmAge) & !is.na(Age)
AgeAccelerationResidual=rep(NA, length(Age) )
if (sum(restNonMissing,na.rm=TRUE) >3 ){
  AgeAccelerationResidual[restNonMissing]=residuals(lm(as.numeric(datout1$DNAmAge)~as.numeric(datSample$Age), 
                                                       subset= restNonMissing))
}
datout1$AgeAccelerationResidual=AgeAccelerationResidual
# Replace - for . in sample names 
datSample$Sample=gsub("-", ".", datSample$Sample) 
dataout1 <- merge(datout1, datSample, by.x = "SampleID", by.y = "Sample" )

###################################
###########Methylation 27K#########
###################################

# 4.Read in the DNA methylation data (beta values)
##To imporve computatinal time Beta matrix filtered for
# probes annotated datMiniAnnotation27k.csv and samples
dat0=read.csv("AC_NT_Meth27_BRCA.csv.gz") 
nSamples=dim(dat0)[[2]]-1
nProbes= dim(dat0)[[1]]
#########################################
#####Run step 5,6, 7 and 8 together!!!!#####
#########################################

#5. Create a log file which will be output into your directory
# The code serves to create a log file (for error checks etc).
# It will automatically create a log file.
file.remove("LogFile.txt")
file.create("LogFile.txt")
DoNotProceed=FALSE
cat(paste( "The methylation data set contains", nSamples, "samples (e.g. arrays) and ", nProbes, " probes."),file="LogFile.txt")
if (nSamples==0) {DoNotProceed=TRUE; cat(paste( "\n ERROR: There must be a data input error since there seem to be no samples.\n Make sure that you input a comma delimited file (.csv file)\n that can be read using the R command read.csv.sql . Samples correspond to columns in that file  ."), file="LogFile.txt",append=TRUE) } 
if (nProbes==0) {DoNotProceed=TRUE; cat(paste( "\n ERROR: There must be a data input error since there seem to be zero probes.\n Make sure that you input a comma delimited file (.csv file)\n that can be read using the R command read.csv.sql  CpGs correspond to rows.")   , file="LogFile.txt",append=TRUE) } 
if (  nSamples > nProbes  ) { cat(paste( "\n MAJOR WARNING: It worries me a lot that there are more samples than CpG probes.\n Make sure that probes correspond to rows and samples to columns.\n I wonder whether you want to first transpose the data and then resubmit them? In any event, I will proceed with the analysis."),file="LogFile.txt",append=TRUE) }
if (  is.numeric(dat0[,1]) ) { DoNotProceed=TRUE; cat(paste( "\n Error: The first column does not seem to contain probe identifiers (cg numbers from Illumina) since these entries are numeric values. Make sure that the first column of the file contains probe identifiers such as cg00000292. Instead it contains ", dat0[1:3,1]  ),file="LogFile.txt",append=TRUE)  } 
if (  !is.character(dat0[,1]) ) {  cat(paste( "\n Major Warning: The first column does not seem to contain probe identifiers (cg numbers from Illumina) since these entries are numeric values. Make sure that the first column of the file contains CpG probe identifiers such as cg00000292. Instead it contains ", dat0[1:3,1]  ),file="LogFile.txt",append=TRUE)  } 
datout=data.frame(Error=c("Input error. Please check the log file for details","Please read the instructions carefully."), Comment=c("", "email Steve Horvath."))
if ( ! DoNotProceed ) {
  nonNumericColumn=rep(FALSE, dim(dat0)[[2]]-1)
  for (i in 2:dim(dat0)[[2]] ){ nonNumericColumn[i-1]=! is.numeric(dat0[,i]) }
  if (  sum(nonNumericColumn) >0 ) { cat(paste( "\n MAJOR WARNING: Possible input error. The following samples contain non-numeric beta values: ", colnames(dat0)[-1][ nonNumericColumn], "\n Hint: Maybe you use the wrong symbols for missing data. Make sure to code missing values as NA in the Excel file. To proceed, I will force the entries into numeric values but make sure this makes sense.\n" ),file="LogFile.txt",append=TRUE)  } 
  XchromosomalCpGs=as.character(probeAnnotation27k$Name[probeAnnotation27k$Chr=="X"])
  selectXchromosome=is.element(dat0[,1], XchromosomalCpGs )
  selectXchromosome[is.na(selectXchromosome)]=FALSE
  meanXchromosome=rep(NA, dim(dat0)[[2]]-1)
  if (   sum(selectXchromosome) >=500 )  {
    meanXchromosome= as.numeric(apply( as.matrix(dat0[selectXchromosome,-1]),2,mean,na.rm=TRUE)) }
  if (  sum(is.na(meanXchromosome)) >0 ) { cat(paste( "\n \n Comment: There are lots of missing values for X chromosomal probes for some of the samples. This is not a problem when it comes to estimating age but I cannot predict the gender of these samples.\n " ),file="LogFile.txt",append=TRUE)  } 
  ## filter by match probes
  match1=match(probeAnnotation21kdatMethUsed$Name , dat0[,1])
  if  ( sum( is.na(match1))>0 ) { 
    missingProbes= probeAnnotation21kdatMethUsed$Name[!is.element( probeAnnotation21kdatMethUsed$Name , dat0[,1])]    
    DoNotProceed=TRUE; cat(paste( "\n \n Input error: You forgot to include the following ", length(missingProbes), " CpG probes (or probe names):\n ", 
                                  paste( missingProbes, sep="",collapse=", ")),file="LogFile.txt",append=TRUE)  }  
  
  ####6. Filter STEP: Restrict the data to 21k probes and ensure they are numeric
  match1=match(probeAnnotation21kdatMethUsed$Name , dat0[,1])
  if  ( sum( is.na(match1))>0 ) stop(paste(sum( is.na(match1)), "CpG probes cannot be matched"))
  dat1= dat0[match1,]
  asnumeric1=function(x) {as.numeric(as.character(x))}
  dat1[,-1]=apply(as.matrix(dat1[,-1]),2,asnumeric1)
  
  #7. Create the output file called datout
  set.seed(1)
  # Do you want to normalize the data (recommended)? 
  #Just run normalizeData=FALSE in order to verify R code is not crashing, 
  #otherwise the results are not going to be realiable. 
  normalizeData=TRUE
  source("StepwiseAnalysis.txt")
  # 8. Output the results 
  if (  sum(  datout$Comment  != "" )   ==0 ) { cat(paste( "\n The individual samples appear to be fine. "),file="LogFile.txt",append=TRUE)  } 
  if (  sum(  datout$Comment != "" )   >0 ) { cat(paste( "\n Warnings were generated for the following samples.\n", datout[,1][datout$Comment != ""], "\n Hint: Check the output file for more details."),file="LogFile.txt",append=TRUE)  } 
} 
#####End run here###
datout2 = datout

#9. Relate DNAm age to chronological age
#To address this task, we read in the sample annotation data that contain the chronological ages.
#Make sure it includes columns: Tissue, Age, Female (Female=1, Male=0). Note the spelling (capitalization).
datSample=read.csv("SampleAnnotation_27_NT_BRCA.csv")
#Order samples in the same way in the output file and sample annotation file
datSample <- datSample[order(datSample$Sample),]
datout2 <- datout2[order(datout2$SampleID),]
DNAmAge=datout2$DNAmAge
medianAbsDev=function(x,y) median(abs(x-y),na.rm=TRUE)
medianAbsDev1=signif(medianAbsDev(DNAmAge, datSample$Age),2)

#10. Does DNAm age acceleration is related to Pam50 or IHC subtypes?
# FIRST ACCELERATION SCORE: The first acceleration meausure is based on the difference.
AgeAccelerationDiff=DNAmAge- datSample$Age
# SECOND ACCELERATION SCORE: residual resulting from regressing DNAmAge on chronological age
Age=datSample$Age
restNonMissing= !is.na(DNAmAge) & !is.na(Age)
AgeAccelerationResidual=rep(NA, length(Age) )
if (sum(restNonMissing,na.rm=TRUE) >3 ){
  AgeAccelerationResidual[restNonMissing]=residuals(lm(as.numeric(datout2$DNAmAge)~as.numeric(datSample$Age), subset= restNonMissing))
}
datout2$AgeAccelerationResidual=AgeAccelerationResidual
# Replace - for . in sample names 
datSample$Sample=gsub("-", ".", datSample$Sample) 
dataout2 <- merge(datout2, datSample, by.x = "SampleID", by.y = "Sample" )
colnames(dataout1) <- colnames(dataout2) 
dataout3=rbind(dataout1, dataout2)
# output the results into the directory
write.table(dataout3,"Output_NT_27_450K.csv", row.names=F, sep="," )


########################
####Plots##############
library(ggplot2)
library(ggpubr)
library(extrafont)
library(grid)
library(EnvStats)
#Load SampleAnnotation file
datSample=read.csv("Output_NT_27_450K.csv",header=TRUE)
### Boxplot of DNAm Age and Age Acceleration between biological groups
##### DNAm AGE plot
y.text=c("DNAm Age \n (Horvath's epigenetic clock)")
p=ggboxplot(datSample, x = "Class", y ="DNAmAge",
            color = "Class", palette = "Set1", #jco, Set1, Spectral
            add = "point", xlab = FALSE, ylab= paste(y.text)) +
  #facet_wrap(~Class_IHC) +
  theme(legend.position="none") 
p= p + stat_compare_means(label = "..p.signif..", method = "wilcox.test", ref.group = "Tumor", hide.ns = TRUE, size= 15)     
####Include total number of samples per group
p= p + stat_n_text(size =5,position = "identity") 
p
ggsave(file=paste(y.text,".tiff"), p + 
         theme(axis.text=element_text(size=15,family="ArialMT"), 
               axis.title=element_text(size=15,face="bold",family="ArialMT"), 
               plot.title = element_text(size=15,family="ArialMT")))

##### Age Acceleration Residual plot
y.text=c("Age Acceleration Residual \n (Horvath's epigenetic clock)")
p=ggboxplot(datSample, x = "Class", y ="AgeAccelerationResidual",
            color = "Class", palette = "Set1", #jco, Set1, Spectral
            add = "point", xlab = FALSE, ylab= paste(y.text)) +
  #facet_wrap(~Class_IHC) +
  theme(legend.position="none") 
p= p + stat_compare_means(label = "..p.signif..", method = "wilcox.test", ref.group = "Tumor", hide.ns = TRUE, size= 15)     
####Agregar N= 3
p= p + stat_n_text(size =5,position = "identity") 
p
ggsave(file=paste(y.text,".tiff"), p + 
         theme(axis.text=element_text(size=15,family="ArialMT"), 
               axis.title=element_text(size=15,face="bold",family="ArialMT"), 
               plot.title = element_text(size=15,family="ArialMT")))

######Correlation plot between Chronological Age vs DNAm Age
library(plyr)
library(reshape2)
library(cowplot)
Name=c("Correlation")
bp= ggplot(datSample, aes(x=DNAmAge, y=Age,color=Class), size = 15)+ geom_point(size = 1.5)+
  scale_colour_manual(values = c( "#928AAC","#BBB9D3","#7089C3","#4F7083"))+
  xlab("DNAm Age") + ylab("Chronological age") +
  theme(legend.position = "none",show.legend.text = F)+
  #facet_wrap(~Class_IHC) +
  theme_classic()
bp

#Get correlation for class all analyzed samples
bp=bp + stat_cor(method = "spearman",col="black", label.x.npc = 0.05, 
                 label.y.npc = 0, hjust = 0)
###Get correlation for class (Normal and Tumor)
p2=bp + stat_cor(aes(color=Class),method = "pearson", show.legend = FALSE)
p2
ggsave(file=paste(Name,".pdf"), p2 + theme(axis.text=element_text(size=25,family="ArialMT"), axis.title=element_text(size=25,face="bold",family="ArialMT"), plot.title = element_text(size=25,family="ArialMT")))

#### Get differences between matched normal adjacent and tumoral tissue
#Load SampleAnnotation file
Y=read.csv("Output_NT_27_450K.csv",header=TRUE)
Y=Y[,1:2]
Y=as.data.frame(Y)
Y_new <- NULL
for (i in 1:nrow(Y))
{
  if(i%%2==1){ 
    
    Y_temp <- Y[i:(i +1), ]
    Y_temp$SampleID <- Y_temp$SampleID[1]
    Y_temp[2, 2] <-  -Y_temp[2, 2]
    Y_sum <- aggregate(Y_temp$DNAmAge, 
                       by=list(a = Y_temp$SampleID),
                       FUN=sum)
    
    Y_new <- rbind(Y_new, Y_sum)
    
    
  }else{'NULL'}
}
colnames(Y_new) <- colnames(Y)
Y_new <- Y_new[order(Y_new$DNAmAge),] 
max=nrow(Y_new)
Y_new$rank=1:96
write.csv(Y_new,"AgeAcelaration_Results.csv", row.names = FALSE)