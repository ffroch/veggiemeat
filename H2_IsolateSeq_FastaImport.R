emptylist2<-readRDS("./06-outputfiles/SangerReadsList.RData")

#make a list of all the fasta files
filelist2<-list.files(path="./06-outputfiles/trimmedSangerfasta/",pattern="*.fa",full.names=T)
#generate namelist for the DNA strings
emptynamelist=vector(mode="list", length=length(filelist2))
j=0
for (file in filelist2) {
  j=j+1
  temporary<-str_remove(file,"./06-outputfiles/trimmedSangerfasta/")
  emptynamelist[[j]]<-str_remove(temporary,".fa")
}
emptynamelist

# generate DNA stringsset from the sequence list and name them with isolate ID
seqstrim<-DNAStringSet(unlist(emptylist2))
seqstrim
names(seqstrim)<-emptynamelist
seqstrim
seqstrim<-OrientNucleotides(seqstrim)


#filter for trimmed sequences >100 bp
seqstrimclean<-seqstrim[width(seqstrim)>=100]

writeXStringSet(seqstrimclean,"./06-outputfiles/seqstrimclean_40_15_100bp.fa")

