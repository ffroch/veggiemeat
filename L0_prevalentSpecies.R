#additional: detect and describe prevalent genera----
#import df of relative abundances per genus and sample.id called relgen
relgen<-readRDS("./06-outputfiles/K0_relgen.rds")
relab_po<-readRDS("./06-outputfiles/E4_relab_po.rds")
#rename column x
colnames(relgen)[3]<-"relab"
#sum relative abundances over all samples for each genus, and filter for those genera with at least a sum of 5% 
mean<-aggregate(relgen$relab,by=list(Genus=relgen$Genus),FUN=function(x) {round(mean(x,na.rm=T),digits=2)})
mean[mean$x>=0.05,]
#generate a list of most "common" genera
commons<-mean[mean$x>=0.05,1]

#exkurs: some general stats for results chapter
relabstats<-data.frame("feature"=NA, "value"=NA)
relabstats[1,]<-c("no. of samples", length(sample_names(relab_po)))
relabstats[2,]<-c("overall reads", sum(otu_table(relab_po)))
relabstats[3,]<-c("min reads", min(sample_sums(relab_po)))
relabstats[4,]<-c("max reads", max(sample_sums(relab_po)))
relabstats[5,]<-c("median reads", median(sample_sums(relab_po)))
relabstats[6,]<-c("no. of asvs", ntaxa(relab_po))
relabstats[7,]<-c("no. of Genera", length(unique(data.frame(tax_table(relab_po))$Genus)))
relabstats[8,]<-c("no. of Phyla", length(unique(data.frame(tax_table(relab_po))$Phylum)))

#find the min an the max relative abundance for each genus among all samples
reltable<-data.frame(otu_table(relab_po))
reltable$sseqid<-row.names(reltable)
reltable$max<-apply(reltable[,-ncol(reltable)],1,max)
reltable$min<-apply(reltable[,-ncol(reltable)],1,min)
#filter for genera with at least 3% relative abundance in at least one sample
min3perc<-reltable[,(ncol(reltable)-2):ncol(reltable)]%>%
  filter(max>=0.03)
#connect this table with taxonomy
reltax<-data.frame(tax_table(phylo_object_red))
reltax$sseqid<-rownames(reltax)
min3perc<-plyr::join(min3perc, reltax, by="sseqid", type="left")
#calculate number of asvs, genera and phyla with at least 3% relative abundance in at least one sample
relabstats[9,]<-c("no. of asvs >3perc", length(unique(min3perc$sseqid)))
relabstats[10,]<-c("no. of Genera >3perc", length(unique(min3perc$Genus)))
relabstats[11,]<-c("no. of Phyla >3perc", length(unique(min3perc$Phylum)))
relabstats[,2]<-as.numeric(relabstats[,2])
relabstats
saveRDS(relabstats, file ="./06-outputfiles/L0_relabstats.rds")

#extract phyla with at least 3% relative abundance in at least 1 sample and calculate min and max relative abundances among the samples
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
rell3pp$min<-as.numeric(rell3pp$min)
saveRDS(rell3pp, file ="./06-outputfiles/L0_relrangedominantphyla.rds")

#calculate genus prevalence (in how many samples is the genus present) and the min and max relative abundances for each genus among the samples
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


saveRDS(prevlist, file ="./06-outputfiles/L0_prevlist.rds")

#find the dominant genus per sample
#1. find the maximum relative abundance
maxpersample<-relgen%>%
  group_by(sample.id)%>%
  summarise(max=max(relab, na.rm=TRUE))
#2. connect it to the corresponding genus
dominantgenera<-data.frame("Genus"=NA, "sample.id"=NA, "relab"=NA)
sampleids<-unique(relgen$sample.id)
for (i in 1:length(sampleids)){
  sampleset<-relgen%>%
    filter(sample.id==sampleids[i])
  maxset<-max(sampleset$relab)
  dominantgenera[i,]<-sampleset%>%
    filter(round(relab,7)==round(maxset,7))
  
}
#count in how many samples the genus is the dominating one
domgen<-data.frame(dominantgenera%>%
  group_by(Genus)%>%
  summarise(n=n()))
domgen<-domgen[order(-domgen$n),]
saveRDS(domgen,file ="./06-outputfiles/L0_domgen.rds" )

#phylum prevalence - in how many sample the phylum has over 50% relative abundance
#import a df with relative abundance per genus and sample (called ohnenull)
ohnenull<-readRDS("./06-outputfiles/K0_ohnenull.rds")
relphyl<-aggregate(ohnenull$relAb,by=list(Phylum=ohnenull$Phylum,sample.id=ohnenull$sample.id),FUN=sum)
colnames(relphyl)[3]<-"relab"
phlyprevlist50<-data.frame(relphyl%>%
                           filter(relab>=0.5)%>%
                         group_by(Phylum)%>%
                           summarise(n=n()))
phlyprevlist50
saveRDS(phlyprevlist50, file ="./06-outputfiles/L0_phylprevlist50.rds")

#phylum prevalence - in how many sample the phylum has over 90% relative abundance
phlyprevlist90<-data.frame(relphyl%>%
                             filter(relab>=0.9)%>%
                             group_by(Phylum)%>%
                             summarise(n=n()))
phlyprevlist90
saveRDS(phlyprevlist90, file ="./06-outputfiles/L0_phylprevlist90.rds")

