#a database with the sequence data from MiSeq run was genereated with 08-CombineIsolatesAndMiSeq/makeblastsdb.sh
#sequences from the isolates were blasted to this database with 08-CombineIsolatesAndMiSeq/blastisolates.sh

#the generated blastout file was corrected for 
completeblast<-read.table("./08-CombineIsolatesAndMiSeq/comparison.blastout")
colnames(completeblast)<-c("isolate.id", "sseqid", "pident", "qcovs", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
bestblast<-read.table("./08-CombineIsolatesAndMiSeq/best.blastout")
colnames(bestblast)<-c("isolate.id", "sseqid", "pident", "qcovs", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
taxaqiime<-data.frame(read_tsv("./08-CombineIsolatesAndMiSeq/taxonomy/taxonomy.tsv",col_names=T))
colnames(taxaqiime)[1]<-"sseqid"
bestfits<-completeblast[completeblast$pident==100,]
n_occur<-data.frame(table(bestfits$isolate.id))
n_occur[n_occur$Freq>1,]
tocheck<-bestfits[bestfits$isolate.id %in% n_occur$Var1[n_occur$Freq >1],]

tocheck<-join(tocheck, taxaqiime, type="left")
isometa<-readRDS("./06-outputfiles/C_isometa.rds")

tochecksamples<-join(tocheck, isometa, by="isolate.id",type="left")
tochecksamples$relab<-0
tochecksamples<-tochecksamples[complete.cases(tochecksamples),]

##chekck, is there already somewhere a relative abundance table
oturel<-data.frame(otu_table(relab_po))

#find relative abundances per asv and sample for the isolates
for (i in 1:nrow(tochecksamples)){
  rowind=tochecksamples[i,2]
  colind=tochecksamples[i,grepl("sample.id",colnames(tochecksamples))]
  tochecksamples$relab[i]=ifelse(is.numeric(oturel[grepl(rowind,rownames(oturel)),grepl(colind,colnames(oturel))]),oturel[grepl(rowind,rownames(oturel)),grepl(colind,colnames(oturel))],NA)
  }

#chose asv with highest rel ab
bestchoice<-data.frame(tochecksamples%>%
    group_by(isolate.id)%>%
    top_n(1,relab))
  #remove lines with 0 percent rel ab
  bestchoice<-bestchoice[bestchoice$relab>0,]
  bestchoice<-bestchoice[complete.cases(bestchoice),]
  bestchoicered<-bestchoice[,1:13]
  bestchoicered$overrule<-1


  #ownbestblast<-data.frame(completeblast%>%
                           #group_by(isolate.id)%>%
                           #top_n(1,pident))
  #ownbestblast$overrule<-0
  bestblast$overrule<-0

  combi<-rbind(bestchoicered, bestblast)
  prefilter<-data.frame(combi%>%
                         group_by(isolate.id)%>%
                         top_n(1,overrule))

  #check whole data for rel abs and remove 0 valued isolates
  prefiltersamples<-join(prefilter, isometa, by="isolate.id",type="left")
  prefiltersamples$relab<-0
  prefiltersamples<-prefiltersamples[complete.cases(prefiltersamples),]
  for (i in 1:nrow(prefiltersamples)){
    rowind=prefiltersamples[i,2]
    colind=prefiltersamples[i,grepl("sample.id",colnames(prefiltersamples))]
    prefiltersamples$relab[i]=ifelse(is.numeric(oturel[grepl(rowind,rownames(oturel)),grepl(colind,colnames(oturel))]),oturel[grepl(rowind,rownames(oturel)),grepl(colind,colnames(oturel))],NA)
  }
  cleanprefilter<-prefiltersamples[prefiltersamples$relab>0,]
  cleanprefilter<-cleanprefilter[complete.cases(cleanprefilter),]

  #remove pident < 98
  final<-cleanprefilter[cleanprefilter$pident>=98,]

write.table(final, "./06-outputfiles/cleaned_blastout.txt")
  
  
write.table(prefiltersamples, "./06-outputfiles/prefiltered_blastout.txt")
  

saveRDS(oturel,"./06-outputfiles/I1_oturel.rds")

