#sample C4 (reads=439) is excluded
data<-readRDS("./06-outputfiles/C_data.rds")
datared<-data[,!grepl("C4",colnames(data))]#removal of sample Veggie20

#rarefaction by size of 2463
if(file.exists("./06-outputfiles/02rarefied-data-BySize.txt"))
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
  bycoverage<-estimateD(datared[,2:ncol(datared)],base="coverage", level=0.995)
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
    filter(Order.q==i, method=="bysize")
  aov<-aov(qD~group,data=dat)
  print(summary(aov))
  print(qqnorm(aov$residuals))
  kw<-kruskal.test(qD~group,data=dat)
  print(kw)
  print(dunn.test(dat$qD,dat$group, method="bonferroni"))
  
}

kwpvals<-data.frame("method"=NA, "pval"=NA)
dunnlist<-list()

for (i in 0:2){
  dat<-hilldiv_m%>%
    filter(Order.q==i, method=="bycoverage")
  aov<-aov(qD~group,data=dat)
  print(summary(aov))
  print(qqnorm(aov$residuals))
  kw<-kruskal.test(qD~group,data=dat)
  print(kw)
  kwpvals[i+1,1]=i
  kwpvals[i+1,2]=kw$p.value
  print(dunn.test(dat$qD,dat$group, method="bonferroni"))
  dunnpvals<-data.frame("method"=i,"pair"=dunn.test(dat$qD,dat$group, method="bonferroni")$comparisons, "adjpval"=dunn.test(dat$qD,dat$group, method="bonferroni")$P.adjusted )
  dunnlist[[i+1]]<-dunnpvals
}
kwpvals
saveRDS(kwpvals, "./06-outputfiles/E1_alphadivkw.rds")
saveRDS(dunnlist, "./06-outputfiles/E1_alphadivdunn.rds")
