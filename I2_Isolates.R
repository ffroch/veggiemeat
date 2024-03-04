relab_po<-readRDS("./06-outputfiles/E4_relab_po.rds")
isometa<-readRDS("./06-outputfiles/C_isometa.rds")
meta<-readRDS("./06-outputfiles/C_meta.rds")
isolates<-read.table("./06-outputfiles/prefiltered_blastout.txt", header=T, sep=" ")
isolates$corrseqid=isolates$sseqid
isolatesn<-isolates[,c(1,ncol(isolates))]
names(isolatesn)<-c("isolate.id","sseqid")
isolatesn1<-plyr::join(isolatesn,isometa, by="isolate.id",type="left")
isolatesn2<-plyr::join(isolatesn1,meta,by="sample.id",type="left")
tax<-data.frame(tax_table(relab_po))
tax$sseqid<-row.names(tax)
saveRDS(tax, "./06-outputfiles/I2_tax.rds")
isolatesn3<-plyr::join(isolatesn2,tax,by="sseqid",type="left")
oturel<-readRDS("./06-outputfiles/I1_oturel.rds")
oturel$sseqid<-row.names(oturel)
isolates<-plyr::join(isolatesn3,oturel,by="sseqid",type="left")

write.table(isolates, "./06-outputfiles/complete-cleaned-filtered-isolates.txt",sep="\t")

#take only isolates with at least Family level
assignedisolates<-isolates%>%
  filter(!is.na(Family))
#generate fasta with only isolates with at least Family level
isolate_ids<-assignedisolates$isolate.id
seqstrimclean<-readDNAStringSet("./06-outputfiles/seqstrimclean_40_15_100bp.fa")
assignedisolatesfasta<-seqstrimclean[isolate_ids]
writeXStringSet(assignedisolatesfasta,"./06-outputfiles/ValidBacterialIsolates.fa")

counts<-data.frame("type"=NA, "counts"=NA)
counts[1,]<-c("overall bacterial isolates", nrow(isolates)+20)#fungi isolates added
counts[2,]<-c("unassigned", nrow(isolates)-nrow(assignedisolates))
counts[3,]<-c("assigned on genus level", nrow(assignedisolates%>%filter(!is.na(Genus))))
counts[4,]<-c("only family level", nrow(assignedisolates%>%filter(is.na(Genus)))-nrow(assignedisolates%>%filter(is.na(Family))))
counts[5,]<-c("number of genera", length(unique(assignedisolates$Genus)))
counts[6,]<-c("number of phyla", length(unique(assignedisolates$Phylum)))
counts[7,]<-c("number of Enterobacteriacea only family level", nrow(assignedisolates%>%filter(is.na(Genus),Family=="Enterobacteriaceae")))
saveRDS(counts, file="./06-outputfiles/I2_counts1.rds")


#count isolates per genus, class and family
length(unique(assignedisolates$Genus))
length(unique(assignedisolates$Class))
length(unique(assignedisolates$Family))


#count isolates per Genus
isolatecounts<-assignedisolates%>%
  group_by(Genus)%>%
  dplyr::summarise(n=n())%>%
  arrange(desc(n))
isolatecounts<-data.frame(isolatecounts)
isolatecounts
write.table(isolatecounts,"./06-outputfiles/isolatespergenus.txt",sep="\t")
prevalence<-data.frame("Genus"=NA,"samplecounts"=NA)
for (i in 1:length(isolatecounts$Genus)){
  check<-assignedisolates%>%
    filter(Genus==isolatecounts$Genus[i])%>%
    group_by(sample.id)
  prevalence[i,1]=isolatecounts$Genus[i]
  prevalence[i,2]=length(unique(check$sample.id))
}
prevalence<-prevalence[order(-prevalence$samplecounts),]
saveRDS(prevalence,"./06-outputfiles/I2_isolateprevalence.rds")

groupprevalence<-data.frame("Genus"=NA,"samplecounts"=NA)
for (i in 1:length(isolatecounts$Genus)){
  check<-assignedisolates%>%
    filter(Genus==isolatecounts$Genus[i])%>%
    group_by(sample.id)
  groupprevalence[i,1]=isolatecounts$Genus[i]
  groupprevalence[i,2]=length(unique(check$group))
}
groupprevalence<-groupprevalence[order(-groupprevalence$samplecounts),]
groupprevalence
saveRDS(groupprevalence,"./06-outputfiles/I2_isolategroupprevalence.rds")




#count isolates per Family
isolatecounts<-assignedisolates%>%
  group_by(Family)%>%
  dplyr::summarise(n=n())%>%
  arrange(desc(n))
isolatecounts<-data.frame(isolatecounts)
isolatecounts

write.table(isolatecounts,"./06-outputfiles/isolatesperfamily.txt",sep="\t")

manufacturerprevalence<-data.frame("Family"=NA,"samplecounts"=NA)
for (i in 1:length(isolatecounts$Family)){
  check<-assignedisolates%>%
    filter(Family==isolatecounts$Family[i])%>%
    group_by(sample.id)
  manufacturerprevalence[i,1]=isolatecounts$Family[i]
  manufacturerprevalence[i,2]=length(unique(check$pseudoprod))
}
manufacturerprevalence<-manufacturerprevalence[order(-manufacturerprevalence$samplecounts),]
manufacturerprevalence


isolatecounts<-assignedisolates%>%
  group_by(Phylum)%>%
  dplyr::summarise(n=n())%>%
  arrange(desc(n))
isolatecounts

write.table(isolatecounts,"./06-outputfiles/isolatesperphylum.txt",sep="\t")

assignedisolates%>%
  filter(Family=="Enterobacteriaceae")

isolatecounts<-assignedisolates%>%
  group_by(Genus,sample.id)%>%
  dplyr::summarise(n=n())%>%
  arrange(desc(n))

saveRDS(isolates, "./06-outputfiles/I2_isolates.rds")