meta<-readRDS("./06-outputfiles/C_meta.rds")
taxmet<-readRDS("./06-outputfiles/J1_taxmet.rds")
#isolate per product
genprod<-taxmet%>%
  group_by(sample.id)%>%
  distinct(Genus,.keep_all=T)
ctgp<-genprod%>%
  group_by(sample.id)%>%
  tally()
genprod$cat<-"product"
gp<-genprod[,c(grep("Kingdom",colnames(genprod)):grep("Genus",colnames(genprod)),grep("sample.id",colnames(genprod)),ncol(genprod))]
col<-colnames(gp)
ccol<-colnames(ctgp)
ctgp$cat<-"product"


#isolate per producer
genfi<-taxmet%>%
  group_by(pseudoprod)%>%
  distinct(Genus,.keep_all=T)
ctgf<-genfi%>%
  group_by(pseudoprod)%>%
  tally()
genfi$cat<-"manufacturer"
gf<-genfi[,c(grep("Kingdom",colnames(genfi)):grep("Genus",colnames(genfi)),grep("pseudoprod",colnames(genfi)),ncol(genfi))]
colnames(gf)<-col
colnames(ctgf)<-ccol
ctgf$cat<-"manufacturer"

#isolate per protein and texture
genpt<-taxmet%>%
  group_by(group)%>%
  distinct(Genus,.keep_all=T)
ctgpt<-genpt%>%
  group_by(group)%>%
  tally()
genpt$cat<-"protein/texture"
gpt<-genpt[,c(grep("Kingdom",colnames(genpt)):grep("Genus",colnames(genpt)),grep("group",colnames(genpt)),ncol(genpt))]
colnames(gpt)<-col
colnames(ctgpt)<-ccol
ctgpt$cat<-"protein/texture"


presentgen<-rbind(gp,gf,gpt)
gencounts<-data.frame(rbind(ctgp,ctgf,ctgpt))

#get isolates per sample
ct<-gmodels::CrossTable(presentgen$Genus,presentgen$sample.id)
ctlist<-data.frame(ct[[1]])
colnames(ctlist)<-c("Genus", "sample.id", "present")
#make present absent 
ctlist$present<-ifelse(ctlist$present==0,0,1)

#add phylodata
ctlist1<-plyr::join(ctlist, presentgen[,1:6],by="Genus", type="left")
#add meta data
ctlist2<-plyr::join(ctlist1, presentgen[,c(7:8)],by="sample.id", type="left")
ctlist3<-unique(ctlist2)
ctlist4<-plyr::join(ctlist3,gencounts[,1:3],by="sample.id", type="left")
ctlist4$comblab<-paste0(ctlist4$sample.id," (",ctlist4$n,")")

#replace unnecessary strings in taxonomy
toclean<-c("Kingdom","Phylum","Class", "Order", "Family", "Genus")
for (i in toclean){
ctlist4[,i]<-sapply(ctlist4[,i],function(x)gsub(".*\\_\\_","",as.character(x)))
}
colnames(ctlist4)[11]<-"cat_2"

#add information for point size which indicates different "species" within a genus
prodspec<-read.table("./06-outputfiles/listofuniqueisolatesperproduct.txt")
prodspec$cat<-"sample.id"
colnames(prodspec)[3]<-"sample.id"
#for manufacturer
compspec<-read.table("./06-outputfiles/listofuniqueisolatesperproducer.txt")
compspec$cat<-"manufacturer"
colnames(compspec)[3]<-"sample.id"
#for group
prottexspec<-read.table("./06-outputfiles/listofuniqueisolatespergroup.txt")
prottexspec$cat<-"protein/texture"
colnames(prottexspec)[3]<-"sample.id"

uniqueisolates<-rbind(prodspec,compspec,prottexspec)
uniqueisolates$Genus<-sapply(uniqueisolates$Genus,function(x)gsub(".*\\_\\_","",as.character(x)))
colnames(uniqueisolates)[2]<-"specct"

ctlist5<-plyr::join(ctlist4,uniqueisolates,by=c("Genus","sample.id"), type="left")



#reduce to unique dataset
ctlist7<-unique(ctlist5[,c(1:2,4:12)])
#count number of genera per ID or group
speccts<- aggregate(ctlist5$specct~ctlist5$Genus+ctlist5$comblab, FUN=sum)
colnames(speccts)<-c("Genus","comblab","specct")
ctlist7<-plyr::join(ctlist7,speccts,by=c("Genus","comblab"),type="left")
#add present absent information
pres<-aggregate(ctlist5$present~ctlist5$Genus+ctlist5$comblab, FUN=max)
colnames(pres)<-c("Genus","comblab","present")
#combine data
ctlist7<-plyr::join(ctlist7,pres,by=c("Genus","comblab"),type="left")

#add relative abundances
##get latest folder for input data for relative abundance data
relAb<-read.table("./06-outputfiles/relAb_per_ASV_product_new.txt",sep="\t")
relAb$Genus<-ifelse(relAb$Family=="Enterobacteriaceae", 
                     ifelse(is.na(relAb$Genus), "unclassified Enterobacteriaceae", relAb$Genus),relAb$Genus )
relAbgen<-aggregate(relAb$relAb~relAb$Genus+relAb$sample.id,FUN=sum)
colnames(relAbgen)<-c("Genus", "sample.id", "relAb")
relAbgen$relAb<-relAbgen$relAb*100
#add rel ab data
ctlist8<-plyr::join(ctlist7,relAbgen,by=c("Genus", "sample.id"),type="left")
ctlist8$relAb<-ifelse(is.na(ctlist8$relAb),0,ctlist8$relAb)



#generate a dataset of genera, presented in MiSeq with at least 10% but without isolates
relAb$relAb<-relAb$relAb*100
relAbfil<-relAb[relAb$relAb>=10,]
comp<-unique(ctlist8$Genus)
newGen<-relAbfil$Genus[!grepl(paste(comp,collapse="|"),relAbfil$Genus)]
relAbmissing<-relAb[grepl(paste(newGen,collapse="|"),relAb$Genus),]
relAbsum<-aggregate(relAbmissing$relAb~relAbmissing$Genus+relAbmissing$sample.id,FUN=sum)
colnames(relAbsum)<-c("Genus", "sample.id", "relAb")
products<-unique(ctlist8$sample.id)
ctlist8a<-expand.grid(Genus=newGen,sample.id=products)
tax<-unique(relAbmissing[,3:8])
ctlist8a<-plyr::join(ctlist8a,tax,type="left")
label<-unique(ctlist8[,c(2,8:10)])
ctlist8a<-plyr::join(ctlist8a,label,type="left")
ctlist8a$comblab<-paste0(ctlist8a$sample.id," (",ctlist8a$n,")")
ctlist8a$specct<-0
ctlist8a$present<-0
ctlist8a<-plyr::join(ctlist8a,relAbsum,type="left")
ctlist8a$relAb<-ifelse(is.na(ctlist8a$relAb),0,ctlist8a$relAb)

#dcapitalize genera without isoaltes
ctlist8a$Genus<-toupper(ctlist8a$Genus)
ctlist8a<-unique(ctlist8a)


#combine both datasets
ctlist9<-rbind(ctlist8,ctlist8a)


#calculate mean relative abundances per group and manufacturer
ctlist9_m<-plyr::join(ctlist9[ctlist9$cat=="product",],meta,by="sample.id",type="left")
meanrelabproducer<-aggregate(ctlist9_m$relAb,list(sample.id=ctlist9_m$pseudoprod,Genus=ctlist9_m$Genus),mean)
meanrelabid<-aggregate(ctlist9_m$relAb,list(sample.id=ctlist9_m$sample.id,Genus=ctlist9_m$Genus),mean)
meanrelabprottex<-aggregate(ctlist9_m$relAb, list(sample.id=ctlist9_m$group,Genus=ctlist9_m$Genus),mean)
meanrelab<-rbind(meanrelabid,meanrelabproducer,meanrelabprottex)
colnames(meanrelab)<-c("sample.id","Genus","meanrelab")

#combine present absent lists ctlist9 with mean relative abundances meanrelab
ctlist10<-plyr::join(ctlist9,meanrelab,by=c("sample.id","Genus"),type="left")
#adapt color palette for relative abundance data with own breaks
ctlist10$logrelAb<-ifelse(is.na(ctlist10$relAb),0,
                         ifelse(ctlist10$relAb==0,0,log10(ctlist10$relAb+1)))
ctlist10$meanrelabclus<-discretize(ctlist10$meanrelab,method="fixed", breaks=c(0,0.00001,0.1,1,10,20,40,60,80,100.1))
#add color palette
pal<-brewer.pal(n=9,"YlOrRd")
pal2<-c("white",pal)

#change wrap title order and labels
ctlist10$cat_f<-factor(ctlist10$cat, levels=c("product", "protein/texture", "manufacturer"))
catnames<-c("product"="product","protein/texture"="protein-\ntexture", "manufacturer"="manufacturer")

#change x axis labels
original<-c("producer", "peacut", "peaminced","soycut", "soyminced")
replacement<-c("M","pea-fibrous", "pea-minced", "soy-fibrous","soy-minced")
for (i in 1:5){
  ctlist10$comblab<-str_replace_all(ctlist10$comblab,original[i],replacement[i])
}




#reroder genera by phylo and change class labels
classlist<-c("Actinobacteria", "Bacteroidia","Deinococci",
             "Bacilli","Alphaproteobacteria","Betaproteobacteria",
             "Gammaproteobacteria")
ctlist10$index2<-ctlist10$Class
reclass<-as.character(c(1:7))
for (i in 1:7){
  ctlist10$index2<-str_replace_all(ctlist10$index2,classlist[i],reclass[i])
}
ctlist10$index2<-as.numeric(ctlist10$index2)
ctlist10<-ctlist10[order(ctlist10$Kingdom,ctlist10$Phylum,ctlist10$index2,ctlist10$Order,ctlist10$Family,ctlist10$Genus),]
ctlist10$index<-c(1:nrow(ctlist10))
class.labs<-c("Actinobacteria\n(Actinomycetota)", "Bacteroidia\n(Bacteroidota)","Deinococci\n(Deinococcota)",
               "Bacilli\n(Bacillota)","Alphaproteobacteria\n(Pseudomonadota)","Betaproteobacteria\n(Pseudomonadota)",
               "Gammaproteobacteria\n(Pseudomonadota)" )
names(class.labs)<-classlist
ctlist10<-ctlist10[complete.cases(ctlist10[,"Class"]),]

#change sample order fitting to new grouporder
meta<-meta[order(meta$protein.source,meta$texture,meta$pseudoprod,meta$sample.id),]
meta$index3<-c(1:nrow(meta))
newindex<-data.frame(sample.id=meta$sample.id,index3=meta$index3)
nip2<-data.frame(sample.id=c("pea-fibrous", "pea-minced", "soy-fibrous", "soy-minced", "producer 1", "producer 2", "producer 3",
                              "producer 4", "producer 5", "producer 6", "producer 7", "producer 8", "producer 9"),index3=c(33,34,35,36,37,38,39,40,41,42,43,44,45))
newindex<-rbind(newindex, nip2)
ctlist10<-plyr::join(ctlist10,newindex,by="sample.id",type="left")


#plot
ctlist10%>%
  filter(cat!="group")%>%
  mutate(Genus=fct_reorder(Genus,index),Class=fct_reorder(Class,index2),comblab=fct_reorder(comblab,index3))%>%
  
  ggplot(aes(comblab,Genus))+
  geom_tile(aes(fill=meanrelabclus),color="lightgrey")+
  geom_point(aes(alpha=factor(present),size=specct),col="#5C5C5C")+
  #geom_text(aes(label=specct),size=2)+
  facet_grid(Class~cat_f,scales="free",space="free",labeller=labeller(cat_f=as_labeller(catnames),Class=class.labs))+
  scale_alpha_manual(values=c(0,1),guide="none")+
  scale_fill_manual( values=pal2,name="relative abundances\nmiseq %")+
  scale_size_continuous(name="estimated different\nisolated strains")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust=0,vjust=0.25),
        axis.text.y=element_text(face="italic"),
        strip.text.y = element_text(angle=0,size = 8), strip.background = element_blank(),
        legend.position = "bottom",
        legend.direction="horizontal", legend.box="horizontal",
        legend.text=element_text(size=8), legend.title=element_text(size=8),
        legend.key.size=unit(3,"mm"))+
  xlab("")

#save plot
ggsave("./04-finalRplots/Fig1inclmeanrelabnewnom20240222.tiff", width=9.53,height=12, units="in",dpi=2000)
ggsave("./VeggieMeat/Fig1inclmeanrelabnewnom20240222.png", width=9.53,height=12, units="in",dpi=300)



