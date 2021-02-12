library(UCSCXenaTools)
options(timeout=1000)
xena_default_hosts()
host <- "https://gdc.xenahubs.net"
#data=XenaData$XenaDatasets
#write(data,"List.data.csv")
dataset <- "TCGA-BRCA.methylation27.tsv"
samples.tab=read.table("Samples_to_Use_TCGA.txt", sep = "\t")
samples=samples.tab[-1,]
samples <- unlist(samples)
probes=read.table("Age_clock_probes.txt")
probes=probes[-1,]
probes <- unlist(probes)
# Fetch samples
fetch_dataset_samples(host, dataset, 2)
# Fetch identifiers
fetch_dataset_identifiers(host, dataset)
# Fetch expression value by probes
data=fetch_dense_values(host, dataset, probes, samples, check = FALSE)
write.table(data,"Meth27_BRCA.txt",sep = "\t")
