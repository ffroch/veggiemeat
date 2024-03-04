#join relative abundance data with taxonomy
oturel<-readRDS("./06-outputfiles/I1_oturel.rds")
oturel$sseqid<-row.names(oturel)
tax<-readRDS("./06-outputfiles/I2_tax.rds")
oturel_tax<-join(oturel, tax, by="sseqid", type="left")
meta<-readRDS("./06-outputfiles/C_meta.rds")
relab_po<-readRDS("./06-outputfiles/E4_relab_po.rds")
isolates<-read.table("./06-outputfiles/complete-cleaned-filtered-isolates.txt",sep="\t")



#generate data structure for plotting, each sample needs the information if there are isolates or not and of the relative abundances
list<-as.character(sample_names(relab_po))
dflist<-vector(mode="list", length=length(list))
i=1  
  
for (i in 1:length(list)){
    set<-isolates[grepl(list[i],isolates$sample.id),
                 c(1:grep("Species",colnames(isolates)),grep(list[i],colnames(isolates)))]
    #set<-na.omit(set)
    colnames(set)[ncol(set)]<-"relAb"
    set<-set%>%add_count(sseqid)
    setsel<-data.frame(isolate.id=set$isolate.id, sseqid=set$sseqid,
                       Kingdom=set$Kingdom,Phylum=set$Phylum,Class=set$Class,
                       Order=set$Order, Family=set$Family, Genus=set$Genus,
                       sample.id=set$sample.id, relAb=set$relAb,n=set$n)
    
    res<-data.frame(isolate.id="no isolate", sseqid=oturel_tax$sseqid,
                    Kingdom=oturel_tax$Kingdom,Phylum=oturel_tax$Phylum,Class=oturel_tax$Class,
                   Order=oturel_tax$Order, Family=oturel_tax$Family, Genus=oturel_tax$Genus,
                  sample.id=list[i], relAb=oturel_tax[,grep(list[i],names(oturel_tax))],n=0)
    resclean<-res[res$relAb>0,]
    comp<-rbind(setsel,resclean)
    
    compclean<-distinct(comp,sseqid,.keep_all=T)
    compclean$iso<-ifelse(compclean$isolate.id=="no isolate","noi","i")
    dflist[[i]]<-compclean
    i=i+1
    #print(sum(compclean$n)) # just to check sums
  }

dflist

completelist<-ldply(dflist,data.frame)
completelist<-completelist[complete.cases(completelist[,"iso"]),]
cl<-completelist
clno<-cl[order(cl$Kingdom,cl$Phylum,cl$Class,cl$Order,cl$Family,cl$Genus),]

#for Fig 1
clno1<-sapply(clno,function(x)gsub(".*\\_\\_","",as.character(x)))


relgen<-aggregate(clno$relAb,by=list(Genus=clno$Genus,clno$sample.id),FUN=sum,drop=F)
relgen$relab<-ifelse(is.na(relgen$x),0,relgen$x)
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



#remove relAb 0 for cleaner presentation, isolates without miseq pendant are disturbing
ohnenull<-clno[clno$relAb>0,]

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

#classes with >3 percent in total per product
relgen<-aggregate(ohnenull$relAb,by=list(Genus=ohnenull$Genus,sample.id=ohnenull$sample.id),FUN=sum)
relgenred<-relgen[relgen$x>0.1,]
genlist<-paste(unique(relgenred$Genus),collapse="|")
relclass<-aggregate(ohnenull$relAb,by=list(class=ohnenull$Class,sample.id=ohnenull$sample.id),FUN=sum)
relclassred<-relclass[relclass$x>0.1,]
classgen<-unique(ohnenull[grepl(genlist, ohnenull$Genus),]$Class)
classgenlist<-paste(classgen,collapse="|")
classlist<-paste(unique(relclassred$class),collapse="|")
relphyl<-aggregate(ohnenull$relAb,by=list(phylum=ohnenull$Phylum,sample.id=ohnenull$sample.id),FUN=sum)
relphylred<-relphyl[relphyl$x>0.1,]
phylgen<-unique(ohnenull[grepl(classlist,ohnenull$Class),]$Phylum)
phylgenlist<-paste(phylgen,collapse="|")
phyllist<-paste(unique(relphylred$phylum),collapse="|")

#rename features if they are low abundant asvs
ohnenull$newtax<-ifelse (grepl(genlist,ohnenull$Genus),ohnenull$Genus,
                         ifelse(grepl(classlist, ohnenull$Class),
                                ifelse(grepl(classgenlist,ohnenull$Class),paste0("other ", ohnenull$Class),ohnenull$Class),
                                ifelse(grepl(phyllist,ohnenull$Phylum),
                                       ifelse(grepl(phylgenlist,ohnenull$Phylum),paste0("other ",ohnenull$Phylum),ohnenull$Phylum),
                                       ifelse(is.na(ohnenull$Phylum),"unassigned", "others"))))


aggregate(ohnenull$relAb,by=list(newtax=ohnenull$newtax),FUN=sum)




#order by taxonomy 
ohnenull<-ohnenull[order(ohnenull$Kingdom,ohnenull$Phylum,ohnenull$Class,ohnenull$Family,ohnenull$newtax),]
ohnenull$index2<-c(1:nrow(ohnenull))


#correct order of genera presentation
indexcorr1<-max(ohnenull[ohnenull$Phylum=="Actinobacteriota",]$index2,na.rm=T)
indexcorr2<-max(ohnenull[ohnenull$Class=="Deinococci",]$index2,na.rm=T)
indexcorr3<-max(ohnenull[ohnenull$Phylum=="Firmicutes",]$index2,na.rm=T)
indexcorr4<-max(ohnenull[ohnenull$Class=="Bacilli",]$index2,na.rm=T)
indexcorr5<-max(ohnenull[ohnenull$Class=="Alphaproteobacteria",]$index2,na.rm=T)
indexcorr6<-max(ohnenull[ohnenull$Class=="Gammaproteobacteria",]$index2,na.rm=T)

#put other ... to end of each class
ohnenull$index3<-ifelse(ohnenull$newtax=="other Actinobacteriota",indexcorr1+1,
                        ifelse(ohnenull$newtax=="other Bacteroidia",indexcorr2,
                               ifelse(ohnenull$newtax=="other Firmicutes",indexcorr3+1,
                                      ifelse(ohnenull$newtax=="other Bacilli",indexcorr4+1,
                                             ifelse(ohnenull$newtax=="other Alphaproteobacteria",indexcorr5+1,
                                                    ifelse(ohnenull$newtax=="other Gammaproteobacteria",indexcorr6+1,
                        ifelse(ohnenull$newtax=="others",ohnenull$index2+10000,
                               ifelse(ohnenull$newtax=="unassigned", ohnenull$index2+100000,ohnenull$index2))))))))







library(RColorBrewer)
palegrey<-adjustcolor("yellow", alpha.f=0.1)

#choose colors for presentation 



pal<-c(brewer.pal(n=11,name="PiYG")[c(1,3)],
       brewer.pal(n=11,name="RdBu")[c(7,9,11)],
       brewer.pal(n=11,name="PRGn")[c(9,10,11)],
       brewer.pal(n=11,name="BrBG")[c(7,8,10)],
       brewer.pal(n=11,name="PRGn")[c(5,3,1)],
       brewer.pal(n=11,name="BrBG")[c(1:4)],
       brewer.pal(n=11,name="RdBu")[c(4:1)],
       "#5C5C5C",
       "lightgrey")

colorlist<-aggregate(relAb~newtax+sample.id,ohnenull, sum)
colorlist<-colorlist[order(-colorlist$relAb),]
colorlist<-distinct(colorlist,newtax,.keep_all=T)
speclist<-ohnenull%>%
  select("Kingdom","Phylum","Class","newtax")
speclist<-unique(speclist)
speclist2<-speclist

speclist2$Kingdom<-ifelse(speclist2$newtax=="others","Z_others", speclist2$Kingdom)
speclist2$Phylum<-ifelse(speclist2$newtax=="others","Z_others", speclist2$Phylum)
speclist2$Class<-ifelse(speclist2$newtax=="others","Z_others", speclist2$Class)
speclist2$Phylum<-ifelse(grepl("p__",speclist2$newtax),speclist2$newtax,speclist2$Phylum)
speclist2$Class<-ifelse(grepl("p__",speclist2$newtax),speclist2$newtax,speclist2$Class)
speclist2$Class<-ifelse(grepl("c__",speclist2$newtax),speclist2$newtax,speclist2$Class)
speclist2<-unique(speclist2)

colorlist<-join(colorlist, speclist2,type="left",by="newtax")
colorlist<-colorlist[order(colorlist$Kingdom,colorlist$Phylum,colorlist$Class,colorlist$newtax),]

i=1
j=1
for (i in 1:nrow(colorlist)){
  colorlist$col[i]<-ifelse(colorlist$relAb[i]<=0,"grey", pal[j] ) 
  if(colorlist$relAb[i]>=0){j=j+1}
  print(colorlist$relAb[i])
  print(j)
}
colorlist
colorlist<-colorlist[,c(1,7)]






ohnenull_m<-join(ohnenull,meta,by="sample.id",type="left")
#index for sorting on x axis
ohnenull_m<-ohnenull_m[order(ohnenull_m$protein.source,ohnenull_m$texture,ohnenull_m$pseudoprod,ohnenull_m$sample.id),]
ohnenull_m$index4<-c(1:nrow(ohnenull_m))



#sorter<-join(sorter,sample_names(relab_po),by="sample.id", type="left")

ohnenull_m$filllabelp1<-ifelse(grepl("other ",ohnenull_m$newtax),"other","")
ohnenull_m$filllabelp2<-ifelse(grepl("other",ohnenull_m$newtax),
                               sub(".* ","",ohnenull_m$newtax),
                               ohnenull_m$newtax)
ohnenull_m$filllabelp2<-gsub("Firmicutes", "Bacillota", ohnenull_m$filllabelp2)


legendlabel<-mixedFontLabel(ohnenull_m$filllabelp1,ohnenull_m$filllabelp2,sep=" ",italic=1:2,always.upright = c("other","others", "Bacteroidia", "Bacilli", "Bacillota", "Alphaproteobacteria", "Gammaproteobacteria", "unassigned"))
names(legendlabel)<-ohnenull_m$newtax

ohnenull_m<-ohnenull_m[complete.cases(ohnenull_m[,"iso"]),]
plot1d<-ohnenull_m%>%
  mutate(sample.id=fct_reorder(sample.id,index4),newtax=fct_reorder(newtax,index3))%>%
  ggplot(aes(sample.id,relAb, fill=newtax, width=0.85),color=palegrey)+
  geom_col(color=ifelse(ohnenull_m$iso=="i", "red",palegrey),
           size=ifelse(ohnenull_m$iso=="i", 1,0.0001),
           alpha=ifelse(ohnenull_m$iso=="i", 1,1))+
  geom_col(color=ifelse(ohnenull_m$iso=="i", "red",palegrey),
           size=ifelse(ohnenull_m$iso=="i", 1,0.0001),
           linetype=ifelse(ohnenull_m$iso=="i", "solid", "blank"),
           alpha=ifelse(ohnenull_m$iso=="i", 1,0))+
  #annotate(geom="text",x=sorter$sample.id,y=-5,label=sorter$reads, size=3,vjust=sorter$vjust)+
  scale_fill_manual(name=NULL,labels=legendlabel,values=colorlist$col)+
  xlab("product")+
  ylab("relative abundances [%]")+
  theme_bw()+
  guides(fill=guide_legend(ncol=4))+
  theme(legend.position="bottom",
        legend.direction="horizontal", legend.box="horizontal",
        legend.text=element_text(size=8,hjust=0), legend.title=element_text(size=8),
        legend.key.size=unit(3,"mm"),
        legend.key.width=unit(2*3,"mm"),
        plot.margin=unit(c(0,0.01,0,0.1),"null"))


#make data structure for plot product group header
ohnenulladd<-data.frame(sample.id=ohnenull_m$sample.id,protein.source=ohnenull_m$protein.source,group=ohnenull_m$group,
                        producer=paste0("M",substr(ohnenull_m$pseudoprod,10,10)),index4=ohnenull_m$index4)
ohnenulladd_l<-data.frame(pivot_longer(ohnenulladd,2:4,names_to="rowclassifier",values_to="colorclassifier"))
ohnenulladd_l$index5<-ifelse(ohnenulladd_l$rowclassifier=="producer",1,
                             ifelse(ohnenulladd_l$rowclassifier=="group",2,3))
#nicepal<-RColorBrewer::brewer.pal(11, "BrBG")[c(11,8,4,1)]
nicepal<-c("#918D30", "#5BD39E", "#AEB733", "#73832D", "#4DCF72", "#68D6C9")
colorlistadd<-c(brewer.pal(11,"BrBG")[c(1:4,7:11)],nicepal[3:6],nicepal[1:2])

#addtext<-c("pea","soy","fibrous","minced","fibrous","minced","M9","M2","M1","M3","M9","M2","M7",
#         "M5","M6","M8","M3","M4","M5")
#yvals<-c(3,3,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1)
#xvals<-c("A4","C6","D6","A5","B4","B3","D6","D3","A1","A5","A6","D1","C3","B4","C5",
#     "C8","B8","B3","B2")
#hjustval<-c(0.5,0.5,0,0,0,0,0.5,0.5,0,0,0.5,0.5,0.5,0,0,0.5,0.5,0.5,0)

addtext<-c("pea","soy","fibrous","minced","fibrous","minced","M9","M2","M1","M3","M9","M2","M7",
           "M5","M6","M8","M3","M4","M5")
yvals<-c(3,3,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1)
xvals<-c("A4","C6","D6","A5","B5","B3","D3","D6","A1","A6","D1","A5","C5","C1","B4",
         "C8","B8","B3","B2")
hjustval<-c(0.5,0.5,0,0.5,0.5,0,0.5,0.5,0,0.5,0.5,0,0,0,0,0.5,0.5,0.5,0)


sorter<-ohnenulladd_l%>%
  mutate(sample.id=fct_reorder(sample.id,index4),colorclassifier=fct_reorder(colorclassifier,index5))
sorter<-sorter[complete.cases(sorter[,"sample.id"]),]
sorter<-data.frame(sample.id=unique(sorter$sample.id))

sorter$vjust<-c(rep(0:1,13),0)
#sorter<-join(sorter,sample_names(relab_po),by="sample.id", type="left")



ohnenulladd_l<-ohnenulladd_l[complete.cases(ohnenulladd_l[,"sample.id"]),]



p1<-ohnenulladd_l%>%
  mutate(sample.id=fct_reorder(sample.id,index4),colorclassifier=fct_reorder(colorclassifier,index5))%>%
  ggplot(aes(sample.id,factor(index5), col=colorclassifier, width=0.85),color=palegrey)+
  geom_point(size=11,shape=15)+
  theme_bw()+
  scale_color_manual(values=colorlistadd)+
  scale_y_discrete(name="attributes",labels=c("manufacturer","texture","main protein"))+
  annotate(geom="text",x=xvals,y=yvals,label=addtext,col="white",hjust=hjustval)+
  theme(legend.position="none",axis.title.x=element_blank(),axis.text.x=element_blank(),
        plot.margin=unit(c(0.01,0.01,0,0.038),"null"))


plot_grid(p1,plot1d,nrow=2,rel_heights = c(1,5))
ggsave("VeggieMeat/Fig2_10percentnewnom20240222.png", width=9.53,height=11, units="in",dpi=300)

