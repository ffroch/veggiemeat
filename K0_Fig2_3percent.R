#join relative abundance data with taxonomy
oturel_tax<-join(oturel, tax, by="sseqid", type="left")

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

write.table(clno,"./06-outputfiles/relAb_per_ASV_product_new.txt",sep="\t")#for Figure 1


#aggregate relative abundances per sample and genus
relgen<-aggregate(clno$relAb,by=list(Genus=clno$Genus,clno$sample.id),FUN=sum,drop=F)
relgen$relab<-ifelse(is.na(relgen$x),0,relgen$x)

#remove relAb 0 for cleaner presentation, isolates without miseq pendant are disturbing
ohnenull<-clno[clno$relAb>0,]

#classes with >3 percent in total per product for a nicer presentation
relgen<-aggregate(ohnenull$relAb,by=list(Genus=ohnenull$Genus,sample.id=ohnenull$sample.id),FUN=sum)
relgenred<-relgen[relgen$x>0.03,]
genlist<-paste(unique(relgenred$Genus),collapse="|")
relclass<-aggregate(ohnenull$relAb,by=list(class=ohnenull$Class,sample.id=ohnenull$sample.id),FUN=sum)
relclassred<-relclass[relclass$x>0.03,]
classgen<-unique(ohnenull[grepl(genlist, ohnenull$Genus),]$Class)
classgenlist<-paste(classgen,collapse="|")
classlist<-paste(unique(relclassred$class),collapse="|")
relphyl<-aggregate(ohnenull$relAb,by=list(phylum=ohnenull$Phylum,sample.id=ohnenull$sample.id),FUN=sum)
relphylred<-relphyl[relphyl$x>0.03,]
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

#put "other"to end of each class
ohnenull$index3<-ifelse(ohnenull$newtax=="other Actinobacteriota",indexcorr1+1,
                        ifelse(ohnenull$newtax=="other Bacteroidia",indexcorr2,
                               ifelse(ohnenull$newtax=="other Firmicutes",indexcorr3+1,
                                      ifelse(ohnenull$newtax=="other Bacilli",indexcorr4+1,
                                             ifelse(ohnenull$newtax=="other Alphaproteobacteria",indexcorr5+1,
                                                    ifelse(ohnenull$newtax=="other Gammaproteobacteria",indexcorr6+1,
                        ifelse(ohnenull$newtax=="others",ohnenull$index2+10000,
                               ifelse(ohnenull$newtax=="unassigned", ohnenull$index2+100000,ohnenull$index2))))))))


#prepare color scheme for nice presentation
library(RColorBrewer)
palegrey<-adjustcolor("yellow", alpha.f=0.1)

#choose colors for presentation 
pal<-c(brewer.pal(n=11,name="PiYG")[c(1:5)],
       brewer.pal(n=11,name="RdBu")[c(7:11)],
       brewer.pal(n=11,name="BrBG")[c(7:11)],
       brewer.pal(n=11,name="PRGn")[c(8:11)],
       brewer.pal(n=11,name="RdYlBu")[c(8:11)],
       brewer.pal(n=11,name="PRGn")[c(1,3)],
       brewer.pal(n=11,name="PuRd")[c(4:6)],
       brewer.pal(n=11,name="BrBG")[c(1:5)],
       brewer.pal(n=11,name="RdBu")[c(5:1)],
       brewer.pal(n=11,name="Spectral")[c(4:7)],
       "black")

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


#add meta data to plot data matrix
ohnenull_m<-join(ohnenull,meta,by="sample.id",type="left")
#index for sorting on x axis first protein source, then texture, then producer and sample id
ohnenull_m<-ohnenull_m[order(ohnenull_m$protein.source,ohnenull_m$texture,ohnenull_m$pseudoprod,ohnenull_m$sample.id),]
ohnenull_m$index4<-c(1:nrow(ohnenull_m))



ohnenull_m$filllabelp1<-ifelse(grepl("other ",ohnenull_m$newtax),"other","")
ohnenull_m$filllabelp2<-ifelse(grepl("other",ohnenull_m$newtax),
                               sub(".* ","",ohnenull_m$newtax),
                               ohnenull_m$newtax)

legendlabel<-mixedFontLabel(ohnenull_m$filllabelp1,ohnenull_m$filllabelp2,sep=" ",italic=1:2,always.upright = c("other","others"))
names(legendlabel)<-ohnenull_m$newtax

ohnenull_m<-ohnenull_m[complete.cases(ohnenull_m[,"iso"]),]
plot1d<-ohnenull_m%>%
  mutate(sample.id=fct_reorder(sample.id,index4),newtax=fct_reorder(newtax,index3))%>%
  ggplot(aes(sample.id,relAb*100, fill=newtax, width=0.85),color=palegrey)+
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
  ylab("relative abundancies [%]")+
  theme_bw()+
  theme(legend.position="bottom",
        legend.direction="horizontal", legend.box="horizontal",
        legend.text=element_text(size=8,hjust=0), legend.title=element_text(size=8),
        legend.key.size=unit(3,"mm"),
        plot.margin=unit(c(0,0.01,0,0.1),"null"))







#make data structure for plot product group header
ohnenulladd<-data.frame(sample.id=ohnenull_m$sample.id,protein.source=ohnenull_m$protein.source,texture=ohnenull_m$texture,
                        producer=paste0("M",substr(ohnenull_m$pseudoprod,10,10)),index4=ohnenull_m$index4)
ohnenulladd_l<-data.frame(pivot_longer(ohnenulladd,2:4,names_to="rowclassifier",values_to="colorclassifier"))
ohnenulladd_l$index5<-ifelse(ohnenulladd_l$rowclassifier=="producer",1,
                             ifelse(ohnenulladd_l$rowclassifier=="texture",2,3))
nicepal<-RColorBrewer::brewer.pal(11, "BrBG")[c(11,8,4,1)]
colorlistadd<-c(brewer.pal(11,"BrBG")[c(1:4,7:11)],nicepal[c(1,2)],nicepal[c(3,4)])

addtext<-c("pea","soy","fibrous","minced","fibrous","minced","M9","M2","M1","M3","M9","M2","M7",
           "M5","M6","M8","M3","M4","M5")
yvals<-c(3,3,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1)
xvals<-c("A4","C6","D6","A5","B4","B3","D6","D3","A1","A5","A6","D1","C3","B4","C5",
        "C8","B8","B3","B2")
hjustval<-c(0.5,0.5,0,0,0,0,0.5,0.5,0,0,0.5,0.5,0.5,0,0,0.5,0.5,0.5,0)

sorter<-ohnenulladd_l%>%
      mutate(sample.id=fct_reorder(sample.id,index4),colorclassifier=fct_reorder(colorclassifier,index5))
sorter<-sorter[complete.cases(sorter[,"sample.id"]),]
sorter<-data.frame(sample.id=unique(sorter$sample.id))

sorter$vjust<-c(rep(0:1,14))
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


ggsave("04-finalRplots/Fig2_3percent.tiff", width=9.53,height=11, units="in",dpi=2000)



#knitr::stitch_rhtml(paste0(scriptname,".R"),output=paste0(outputdir,scriptname,"_",Sys.Date(),".html"))



#for poster
clus<-read.table("H:/0300LMM/0302VeggieProducts/06Rnew/00_MetaData/tSNEcluster2.txt",sep="\t", header=T)

ohnenullm1<-join(ohnenull_m,clus,by="ID", type="left")

ohnenullm1$index5<-ifelse(ohnenullm1$Cluster=="Leu", ohnenullm1$index4+1000000,0)
ohnenullm1$index5<-ifelse(ohnenullm1$Cluster=="Lac", ohnenullm1$index4+3000000,ohnenullm1$index5)
ohnenullm1$index5<-ifelse(ohnenullm1$Cluster=="Pro", ohnenullm1$index4+2000000,ohnenullm1$index5)

sorter2<-data.frame(ID=unique(ohnenullm1$ID[order(ohnenullm1$index5)]))
sorter2$vjust<-c(rep(0:1,13),0)
sorter2<-join(sorter2,samplereads,by="ID", type="left")

pd2<-ohnenullm1%>%
  mutate(ID=fct_reorder(ID,index5),newtax=fct_reorder(newtax,index3))%>%
  ggplot(aes(ID,relAb, fill=newtax, width=0.85),color=palegrey)+
  geom_col(color=ifelse(ohnenullm1$iso=="i", "red",palegrey),
           size=ifelse(ohnenullm1$iso=="i", 1,0.0001),
           alpha=ifelse(ohnenullm1$iso=="i", 1,1))+
  geom_col(color=ifelse(ohnenullm1$iso=="i", "red",palegrey),
           size=ifelse(ohnenullm1$iso=="i", 1,0.0001),
           linetype=ifelse(ohnenullm1$iso=="i", "solid", "blank"),
           alpha=ifelse(ohnenullm1$iso=="i", 1,0))+
  annotate(geom="text",x=sorter2$ID,y=-5,label=sorter2$reads, size=3,vjust=sorter2$vjust)+
  scale_fill_manual(name=NULL,labels=legendlabel,values=colorlist$col)+
  xlab("product")+
  ylab("relative abundancies [%]")+
  theme_bw()+
  theme(legend.position="bottom",
        legend.direction="horizontal", legend.box="horizontal",
        legend.text=element_text(size=8,hjust=0), legend.title=element_text(size=8),
        legend.key.size=unit(3,"mm"),
        legend.background = element_rect(fill="white"),
        plot.margin=unit(c(0,0.01,0,0.1),"null"),
        plot.background=element_rect(fill="white"),
        panel.background = element_rect(fill="white"))
pd2

ohnenulladd1<-data.frame(ID=ohnenullm1$ID,main_prot=ohnenullm1$main_prot,category=ohnenullm1$category,
                        producer=paste0("M",substr(ohnenullm1$pseudoprod,10,10)),index5=ohnenullm1$index5)
ohnenulladd1_l<-pivot_longer(ohnenulladd1,2:4,names_to="rowclassifier",values_to="colorclassifier")
ohnenulladd1_l$index6<-ifelse(ohnenulladd1_l$rowclassifier=="producer",1,
                             ifelse(ohnenulladd1_l$rowclassifier=="category",2,3))
nicepal<-RColorBrewer::brewer.pal(11, "BrBG")[c(11,8,4,1)]
colorlistadd<-c(brewer.pal(11,"BrBG")[c(1:4,7:11)],nicepal[c(1,3)],nicepal[c(2,4)])

addtext<-c("pea","soy","pea","soy",
           "fibrous","minced","fibrous","minced","fib.","min.","fib.","min.",
           "M9","M2","M7","M6","M1","M3","M9","M2","M5","M6","M3","M4","M5","M8","M5")
yvals<-c(3,3,3,3,
         2,2,2,2,2,2,2,2,
         1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
xvals<-c("D8","C6","A6","B3",
         "D4","A5","C6","A4","C1","B8","C3","B2",
         "D3","A3","C5","B5","A1","A6","D1","A2","C1","B4","B8","B3","C3","C8","B2")
hjustval<-c(0.5,0.5,0,0,
            0,0,0.5,0,0,0,0,0,
            0.5,0.5,0,0.5,0,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0)

nicepal<-c("#006957","#A55300","#E98219","#00392F","#10937D","#147591")
colorlistadd<-c(brewer.pal(9,"PuBuGn")[5:9],brewer.pal(9,"YlOrBr")[6:9],nicepal[c(3,2)],nicepal[c(1,4)])

p2<-ohnenulladd1_l%>%
  mutate(ID=fct_reorder(ID,index5),colorclassifier=fct_reorder(colorclassifier,index6))%>%
  ggplot(aes(ID,factor(index6), col=colorclassifier, width=0.85),color=palegrey)+
  geom_point(size=11,shape=15)+
  theme_bw()+
  scale_color_manual(values=colorlistadd)+
  scale_y_discrete(name="attributes",labels=c("manufacturer","texture","main protein"))+
  annotate(geom="text",x=xvals,y=yvals,label=addtext,col="white",hjust=hjustval)+
  theme(legend.position="none",axis.title.x=element_blank(),axis.text.x=element_blank(),
        plot.margin=unit(c(0.01,0.01,0,0.038),"null"),
        plot.background=element_rect(fill="white"),
        panel.background = element_rect(fill="white"))
p2


#ohnenulladd1_l%>%
 # mutate(ID=fct_reorder(ID,index5),colorclassifier=fct_reorder(colorclassifier,index6))%>%
  #ggplot(aes(ID,factor(index6), col=colorclassifier, width=0.85),color=palegrey)+
  #geom_point(size=11,shape=15)+
  #theme_bw()+
  #scale_color_manual(values=colorlistadd)+
  #scale_y_discrete(name="attributes",labels=c("manufacturer","texture","main protein"))+
  #annotate(geom="text",x=xvals,y=yvals,label=addtext,col="white",hjust=hjustval)+
  #theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
   #     plot.margin=unit(c(0.01,0.01,0,0.038),"null"),
    #    legend.background=element_rect("white"))

plot_grid(p2,pd2,nrow=2,rel_heights = c(1,5))
ggsave(paste0(outputdir,"Fig2.tiff"), width=45,height=52, units="cm",dpi=1200)

ragg::agg_png(paste0(outputdir,"Fig2ragg.png"), width=14.43,height=12.53, units="in",res=1200,scaling=1.5)
plot_grid(p2,pd2,nrow=2,rel_heights = c(1,5))
dev.off()
