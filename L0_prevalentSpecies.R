#additional: detect and describe prevalent genera----
colnames(relgen)[3]<-"relab"
mean<-aggregate(relgen$relab,by=list(Genus=relgen$Genus),FUN=function(x) {round(mean(x,na.rm=T),digits=2)})
mean[mean$x>=0.05,]
commons<-mean[mean$x>=0.05,1]

for (i in 1:length(commons)){
  mins<-min(relgen%>%
              filter(Genus==commons[i],relab>0)%>%select(relab))
  print(paste(commons[i],mins,sep="-"))
} 

for (i in 1:length(commons)){
  maxs<-max(relgen%>%
              filter(Genus==commons[i])%>%select(relab))
  print(paste(commons[i],maxs,sep="-"))
} 
for (i in 1:length(commons)){
  selection<-ohnenull%>%filter(Genus==commons[i])%>%select(Genus, sample.id)
  nrows<-nrow(unique(selection))
  print(paste(commons[i],nrows,sep="-"))
} 
allgen<-mean[,1]
df<-data.frame(Genus=NULL,n=NULL)
for (i in 1:length(allgen)){
  selection<-ohnenull%>%filter(Genus==allgen[i])%>%select(Genus, sample.id)
  nrows<-nrow(unique(selection))
  df[i,1]<-allgen[i]
  df[i,2]<-nrows
} 
df[order(-df$V2),]

#excurs some stats for results
relabstats<-data.frame("feature"=NA, "value"=NA)
relabstats[1,]<-c("no. of samples", length(sample_names(relab_po)))
relabstats[2,]<-c("overall reads", sum(otu_table(relab_po)))
relabstats[3,]<-c("min reads", min(sample_sums(relab_po)))
relabstats[4,]<-c("max reads", max(sample_sums(relab_po)))
relabstats[5,]<-c("median reads", median(sample_sums(relab_po)))
relabstats[6,]<-c("no. of asvs", ntaxa(relab_po))
relabstats[7,]<-c("no. of Genera", length(unique(data.frame(tax_table(relab_po))$Genus)))
relabstats[8,]<-c("no. of Phyla", length(unique(data.frame(tax_table(relab_po))$Phylum)))

reltable<-data.frame(otu_table(relab_po))
reltable$sseqid<-row.names(reltable)
reltable$max<-apply(reltable[,-ncol(reltable)],1,max)
reltable$min<-apply(reltable[,-ncol(reltable)],1,min)
min3perc<-reltable[,(ncol(reltable)-2):ncol(reltable)]%>%
  filter(max>=0.03)
reltax<-data.frame(tax_table(phylo_object))
reltax$sseqid<-rownames(reltax)
min3perc<-plyr::join(min3perc, reltax, by="sseqid", type="left")

relabstats[9,]<-c("no. of asvs >3perc", length(unique(min3perc$sseqid)))
relabstats[10,]<-c("no. of Genera >3perc", length(unique(min3perc$Genus)))
relabstats[11,]<-c("no. of Phyla >3perc", length(unique(min3perc$Phylum)))
relabstats[,2]<-as.numeric(relabstats[,2])
relabstats
saveRDS(relabstats, file ="./06-outputfiles/relabstats.rds")

reltable2<-reltable
reltable2$sseqid<-rownames(reltable2)
reltable2<-join(reltable2, reltax, by="sseqid", type="left")
larger3percentphyla<-unique(min3perc$Phylum)
rell3pp<-data.frame("Phylum"=NA, "max"=NA, "min"=NA)
for ( i in 1:length(larger3percentphyla)){
  redtab<-reltable2%>%
    filter(Phylum==larger3percentphyla[i])
  rell3pp[i,1]=larger3percentphyla[i]
  rell3pp[i,2]=max(redtab$max)
  rell3pp[i,3]=min(redtab$min)
  
}
rell3pp
saveRDS(rell3pp, file ="./06-outputfiles/relrangedominantphyla.rds")

#stats of the most common genera median >5%
median<-aggregate(relgen$relab,by=list(Genus=relgen$Genus),FUN=function(x) {round(median(x,na.rm=T),digits=2)})
median[median$x>=0.05,]
commongenera<-mean[mean$x>=0.05,1]
larger3percentphyla<-unique(min3perc$Phylum)
rell3pp<-data.frame("Phylum"=NA, "max"=NA, "min"=NA)
for ( i in 1:length(larger3percentphyla)){
  redtab<-reltable2%>%
    filter(Phylum==larger3percentphyla[i])
  rell3pp[i,1]=larger3percentphyla[i]
  rell3pp[i,2]=max(redtab$max)
  rell3pp[i,3]=min(redtab$min)
  
}
rell3pp
rell3pp$min<-as.numeric(rell3pp$min)
saveRDS(rell3pp, file ="./06-outputfiles/relrangedominantphyla.rds")

#genus prevalence
prevlist<-data.frame(relgen%>%
  group_by(Genus)%>%
  summarise(n=n())
)
prevlist<-prevlist[order(-prevlist$n),]
prevlist
prevlist$min<-NA
prevlist$max<-NA

for (i in 1:length(prevlist$Genus)){
  relrange=relgen%>%filter(Genus==prevlist$Genus[i])
  prevlist$min[i]=min(relrange$relab)*100
  prevlist$max[i]=max(relrange$relab)*100
}

prevlist


saveRDS(prevlist, file ="./06-outputfiles/prevlist.rds")

maxpersample<-relgen%>%
  group_by(sample.id)%>%
  summarise(max=max(relab, na.rm=TRUE))

dominantgenera<-data.frame("Genus"=NA, "sample.id"=NA, "relab"=NA)
sampleids<-unique(relgen$sample.id)
for (i in 1:length(sampleids)){
  sampleset<-relgen%>%
    filter(sample.id==sampleids[i])
  maxset<-max(sampleset$relab)
  dominantgenera[i,]<-sampleset%>%
    filter(round(relab,7)==round(maxset,7))
  
}

domgen<-data.frame(dominantgenera%>%
  group_by(Genus)%>%
  summarise(n=n()))
domgen<-domgen[order(-domgen$n),]
saveRDS(domgen,file ="./06-outputfiles/domgen.rds" )

#phylum prevalence
relphyl<-aggregate(ohnenull$relAb,by=list(Phylum=ohnenull$Phylum,sample.id=ohnenull$sample.id),FUN=sum)
colnames(relphyl)[3]<-"relab"
phlyprevlist50<-data.frame(relphyl%>%
                           filter(relab>=0.5)%>%
                         group_by(Phylum)%>%
                           summarise(n=n()))
phlyprevlist50
saveRDS(phlyprevlist50, file ="./06-outputfiles/phylprevlist50.rds")

phlyprevlist90<-data.frame(relphyl%>%
                             filter(relab>=0.9)%>%
                             group_by(Phylum)%>%
                             summarise(n=n()))
phlyprevlist90
saveRDS(phlyprevlist90, file ="./06-outputfiles/phylprevlist90.rds")

