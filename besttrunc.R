library(stringr)
library(tidyr)
library(ggplot2)
filelist<-dir("./03-qiime/05-dada2_trial",pattern="\\.tsv$",full.names=T,recursive=TRUE)

sumstats<-data.frame("mean_input_passed"=NA, "median_input_passed"=NA, "mean_input_merged"=NA, "median_input_merged"=NA,"mean_non_chimeric"=NA, "median_non_chimeric"=NA, "trunc"=NA )

k=1
for (file in filelist){
data<-read.table(file, skip=2, sep="\t")
names(data)<-read.table(file, nrow=1, sep="\t", comment.char="")[1,]
stats<-data.frame("mean_input_passed"=NA, "median_input_passed"=NA, "mean_input_merged"=NA, "median_input_merged"=NA,"mean_non_chimeric"=NA, "median_non_chimeric"=NA )
j=1
for (i in c(4,7,9)){
  stats[1,j]<-mean(data[,i])
  stats[1,j+1]<-median(data[,i])
  j=j+2
}
name<-str_remove(file,"./03-qiime/05-dada2_trial/stats-dada2_")
name<-str_remove(name, "/metadata.tsv")
sumstats[k,]<-stats[1,]
sumstats[k,7]<-name
k=k+1
}

datalong<-pivot_longer(sumstats,1:6,names_to="feature", values_to="percentage")

ggplot(datalong, aes(trunc, percentage))+
  geom_boxplot()+
  theme(axis.text.x=element_text(angle=90))
  
