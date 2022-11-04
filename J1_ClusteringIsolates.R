placeholder<-data.frame(isolate.id=names(assignedisolatesfasta))
placeholder<-join(placeholder, isometa, by="isolate.id", type="left")
placeholder<-join(placeholder, meta, by="sample.id", type="left")


taxtab<-read.table("./06-outputfiles/complete-cleaned-filtered-isolates.txt",sep="\t")
rownames(taxtab)<-taxtab$isolate.id
taxmet<-taxtab[complete.cases(taxtab$sample.id),]


samplelist<-unique(placeholder$sample.id[!is.na(placeholder$sample.id)])
emptylist<-list()

for (i in 1:length(samplelist)){

adjustment<-data.frame(isolate.id=placeholder$isolate.id[grepl(samplelist[i],placeholder$sample.id)])
ssa1<-assignedisolatesfasta[adjustment$isolate.id]
ssa1

aligned<-AlignSeqs(ssa1)

writeXStringSet(aligned, file = "./06-outputfiles/tempal.fasta")

dna<- read.alignment("./06-outputfiles/tempal.fasta", format = "fasta")
D <- dist.alignment(dna, matrix = "identity")
labs<-attr(D, "Labels")
Dn<-as.dist(additive(D))
attr(Dn, "Labels") <-labs
clus<-IdClusters(Dn,cutoff=0.06,showPlot=T)
cpc<-data.frame(clus)

ccpc<-table(cpc)
cpc$isolate.id<-rownames(cpc)
colnames(cpc)[1]<-"counts"


comb<-join(cpc, taxmet, by="isolate.id", type="inner")
combred<-data.frame(comb$counts, comb$Genus)
comuni<-data.frame(unique(combred))
gencounts<-comuni%>%
  group_by(comuni[,2])%>%
  tally()
gencounts<-data.frame(gencounts)
gencounts$sample.id<-samplelist[i]
colnames(gencounts)[1]<-"Genus"

emptylist[[i]]<-gencounts

}

completelist<-ldply(emptylist,data.frame)

write.table(completelist, "./06-outputfiles/listofuniqueisolatesperproduct.txt")


#calculate unique isolates per sample group----
grouplist<-unique(placeholder$group[!is.na(placeholder$group)])


emptylist<-list()

for (i in 1:length(grouplist)){
  
 adjustment<-data.frame(isolate.id=placeholder$isolate.id[grepl(grouplist[i],placeholder$group)])
  ssa1<-assignedisolatesfasta[adjustment$isolate.id]
  ssa1
  
  aligned<-AlignSeqs(ssa1)
  
  writeXStringSet(aligned, file = "./06-outputfiles/tempal.fasta")
  
  dna<- read.alignment("./06-outputfiles/tempal.fasta", format = "fasta")
  D <- dist.alignment(dna, matrix = "identity")
  labs<-attr(D, "Labels")
  Dn<-as.dist(additive(D))
  attr(Dn, "Labels") <-labs
  clus<-IdClusters(Dn,cutoff=0.06,showPlot=T)
  cpc<-data.frame(clus)
  
  ccpc<-table(cpc)
  cpc$isolate.id<-rownames(cpc)
  colnames(cpc)[1]<-"counts"
  
  comb<-join(cpc, taxmet, by="isolate.id", type="inner")
  combred<-data.frame(comb$counts, comb$Genus)
  comuni<-data.frame(unique(combred))
  gencounts<-comuni%>%
    group_by(comuni[,2])%>%
    tally()
  gencounts<-data.frame(gencounts)
  gencounts$samplegroup<-grouplist[i]
  colnames(gencounts)[1]<-"Genus"
  
  emptylist[[i]]<-gencounts
  
}

completelist<-ldply(emptylist,data.frame)

write.table(completelist, "./06-outputfiles/listofuniqueisolatespergroup.txt")


#calculate unique isolates per producer----
producerlist<-unique(placeholder$pseudoprod[!is.na(placeholder$pseudoprod)])
emptylist<-list()

for (i in 1:length(producerlist)){
  
  adjustment<-data.frame(isolate.id=placeholder$isolate.id[grepl(producerlist[i],placeholder$pseudoprod)])
  ssa1<-assignedisolatesfasta[adjustment$isolate.id]
  ssa1
  
  aligned<-AlignSeqs(ssa1)
  
  writeXStringSet(aligned, file = "./06-outputfiles/tempal.fasta")
  
  dna<- read.alignment("./06-outputfiles/tempal.fasta", format = "fasta")
  D <- dist.alignment(dna, matrix = "identity")
  labs<-attr(D, "Labels")
  Dn<-as.dist(additive(D))
  attr(Dn, "Labels") <-labs
  clus<-IdClusters(Dn,cutoff=0.06,showPlot=T)
  cpc<-data.frame(clus)
  
  ccpc<-table(cpc)
  cpc$isolate.id<-rownames(cpc)
  colnames(cpc)[1]<-"counts"
  
  comb<-join(cpc, taxmet, by="isolate.id", type="inner")
  combred<-data.frame(comb$counts, comb$Genus)
  comuni<-data.frame(unique(combred))
  gencounts<-comuni%>%
    group_by(comuni[,2])%>%
    tally()
  gencounts<-data.frame(gencounts)
  gencounts$pseudoprod<-producerlist[i]
  colnames(gencounts)[1]<-"Genus"
  
  emptylist[[i]]<-gencounts
  
}

completelist<-ldply(emptylist,data.frame)

write.table(completelist, "./06-outputfiles/listofuniqueisolatesperproducer.txt")
