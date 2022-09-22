#sample C4 (reads=439) is excluded
datared<-data[,!grepl("C4",colnames(data))]#removal of sample Veggie20

#rarefaction by size of 2463
if(file.exists("./06-outputfiles/02rarefied-data-ByCoverage.txt"))
{
  bysize<-read.table("./06-outputfiles/02rarefied-data-BySize.txt",header=T, sep=" ")
}   else   {
  bysize<-estimateD(datared[,2:ncol(datared)],base="size", level=2463)
  bysize$method="bysize"
  write.table(bysize,"./06-outputfiles/02rarefied-data-BySize.txt")
}

#rarefaction by coverage of 0.995
if(file.exists("./06-outputfiles/02rarefied-data-ByCoverage.txt"))
{
  bycoverage<-read.table("./06-outputfiles/02rarefied-data-ByCoverage.txt",header=T, sep=" ")
}   else   {
  bycoverage<-estimateD(data[,2:ncol(data)],base="coverage", level=0.995)
  bycoverage$method="bycoverage"
  write.table(bycoverage,"./06-outputfiles/02rarefied-data-ByCoverage.txt")
}

hilldiv<-rbind(bysize,bycoverage)
colnames(hilldiv)[1]<-"sample.id"
hilldiv_m<-join(hilldiv, meta, by="sample.id", type="left")

#alpha-diversity----
#comparison of alpha diversity among the predefined groups----
for (i in 0:2){
  dat<-hilldiv_m%>%
    filter(order==i, method=="bysize")
  aov<-aov(qD~group,data=dat)
  print(summary(aov))
  print(qqnorm(aov$residuals))
  kw<-kruskal.test(qD~group,data=dat)
  print(kw)
  print(dunn.test(dat$qD,dat$group, method="bonferroni"))
  
}

for (i in 0:2){
  dat<-hilldiv_m%>%
    filter(order==i, method=="bycoverage")
  aov<-aov(qD~group,data=dat)
  print(summary(aov))
  print(qqnorm(aov$residuals))
  kw<-kruskal.test(qD~group,data=dat)
  print(kw)
  print(dunn.test(dat$qD,dat$group, method="bonferroni"))
}