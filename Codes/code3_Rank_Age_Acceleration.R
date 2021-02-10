#Segmented regression: define the cutoff
library(segmented)
delta=read.csv("AgeAcelaration_Results.csv", header=T)
#### Filter outliers by quiantile (<5% and >95%)
DNAmAge=delta$DNAmAge
quantile(DNAmAge, seq(0, 1, by=0.005))
### Define class based on quantile
dati=delta
dati=dati[which(dati$DNAmAge< 27.02570202),]
os<-segmented(lm(DNAmAge~rank,data=dati),seg.Z=~rank, 
              control=seg.control(n.boot=0,it.max = 1000,K = 10))
os
#plot.segmented(os)
summary(delta)
#### Assign age accelerated classes based on segmented regression and mean value 
delta$Class_Acel_Class <- cut(delta$DNAmAge, breaks = c(-Inf, 4.337 , 14.59 , Inf), 
                         labels = c("Low", "Medium", "High"))
write.csv(delta,"AgeAcelaration_Results.csv", row.names = FALSE)
