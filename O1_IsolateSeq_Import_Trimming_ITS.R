#find all ab1 files from bacterial isolates (27F) and make list
filelist<-intersect(list.files(path="./02-rawdata/SangerSeqData",pattern="*.ab1",full.names=T),
                    list.files(path="./02-rawdata/SangerSeqData",pattern="*ITS",full.names=T))


#generate DNA strings from all the ab1 files incl. trimming and generate fasta files 
emptylist=vector(mode="list", length=length(filelist))
i=0
for (file in filelist){
  vegprod1 <- SangerRead(readFeature = "Forward Read", readFileName=file, TrimmingMethod = "M2", M2CutoffQualityScore = 40, M2SlidingWindowSize = 15, M1TrimmingCutoff = NULL)
  if (vegprod1@QualityReport@trimmedMeanQualityScore>35){
  string<-vegprod1@primarySeqRaw
  trims<-vegprod1@QualityReport@trimmedStartPos
  trime<-vegprod1@QualityReport@trimmedFinishPos
  trimseq<-subseq(string, start=trims, end=trime)
  print(trimseq)
  i=i+1
  emptylist[[i]]<-trimseq
  writeFastaSR(vegprod1, outputDir= "./06-outputfiles/trimmedSangerfasta_ITS/")
  }
}
emptylist
emptylist2<-emptylist[-which(sapply(emptylist, is.null))]

saveRDS(emptylist2, file="./06-outputfiles/SangerReadsList_ITS.RData")
