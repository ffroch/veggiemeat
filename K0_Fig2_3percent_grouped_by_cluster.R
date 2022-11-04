#for poster
ohnenull_mp<-ohnenull_m
ohnenull_mp$index5<-ifelse(ohnenull_mp$profile=="Le", ohnenull_mp$index4+1000000,0)
ohnenull_mp$index5<-ifelse(ohnenull_mp$profile=="La", ohnenull_mp$index4+3000000,ohnenull_mp$index5)
ohnenull_mp$index5<-ifelse(ohnenull_mp$profile=="Pro", ohnenull_mp$index4+2000000,ohnenull_mp$index5)


pd2<-ohnenull_mp%>%
  mutate(sample.id=fct_reorder(sample.id,index5),newtax=fct_reorder(newtax,index3))%>%
  ggplot(aes(sample.id,relAb, fill=newtax, width=0.85),color=palegrey)+
  geom_col(color=ifelse(ohnenull_mp$iso=="i", "red",palegrey),
           size=ifelse(ohnenull_mp$iso=="i", 1,0.0001),
           alpha=ifelse(ohnenull_mp$iso=="i", 1,1))+
  geom_col(color=ifelse(ohnenull_mp$iso=="i", "red",palegrey),
           size=ifelse(ohnenull_mp$iso=="i", 1,0.0001),
           linetype=ifelse(ohnenull_mp$iso=="i", "solid", "blank"),
           alpha=ifelse(ohnenull_mp$iso=="i", 1,0))+
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


#unfinished
ohnenulladd1<-data.frame(ID=ohnenull_mp$ID,main_prot=ohnenull_mp$main_prot,category=ohnenull_mp$category,
                        producer=paste0("M",substr(ohnenull_mp$pseudoprod,10,10)),index5=ohnenull_mp$index5)
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
