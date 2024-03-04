#rarefaction by coverage of 0.995
if(file.exists("./06-outputfiles/02rarefied-data-ByCoverage.txt"))
  {
  bycoverage<-read.table("./06-outputfiles/02rarefied-data-ByCoverage.txt",header=T, sep=" ")
   }   else   {
     bycoverage<-estimateD(datared[,2:ncol(datared)],base="coverage", level=0.995)
     bycoverage$method="bycoverage"
     write.table(bycoverage,"./06-outputfiles/02rarefied-data-ByCoverage.txt")
     }



hilldiv<-bycoverage
colnames(hilldiv)[1]<-"sample.id"
hilldiv_m<-join(hilldiv, meta, by="sample.id", type="left")

saveRDS(hilldiv_m, "./06-outputfiles/E2_hilldiv.RDS")

#alpha-diversity----
#comparison of alpha diversity among the predefined groups----
pal<-c("#AEB733", "#73832D", "#4DCF72", "#68D6C9")

for (i in 0:2){
  dat<-hilldiv_m%>%
    filter(Order.q==i, method=="bycoverage")
  aov<-aov(qD~group,data=dat)
  print(summary(aov))
  print(qqnorm(aov$residuals))
  kw<-kruskal.test(qD~group,data=dat)
  print(kw)
  print(dunn.test(dat$qD,dat$group, method="bonferroni"))
}

p1<-hilldiv_m%>%
  filter(Order.q==0, method=="bycoverage")%>%
  ggplot(aes(group,qD,col=group))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(size=6,alpha=0.3,position=position_jitter(width=0.3,seed=4))+
  scale_color_brewer(palette="Paired")+
  theme_classic()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())+
  ylab("species richness")+
  labs(color="group")+
  geom_segment(aes(x=1,y=max(qD)*0.99,xend=2,yend=max(qD)*0.99),col="black")+
  
  geom_text(aes(label="*",x=1.5,y=max(qD)),col="black")+

  geom_text(aes(group,qD,label=sample.id),position=position_jitter(width=0.3,seed=4),size=2,color="black")
p1

p2<-hilldiv_m%>%
  filter(Order.q==1, method=="bycoverage")%>%
  ggplot(aes(group,qD, color=group))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(size=6,alpha=0.3,position=position_jitter(width=0.3,seed=4))+
  scale_color_manual(values=pal)+
  theme_classic()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())+
  ylab("Hill-Shannon")+
  labs(color="group")+
  geom_segment(aes(x=1,y=max(qD)*0.99,xend=2,yend=max(qD)*0.99),col="black")+

  geom_text(aes(label="*",x=1.5,y=max(qD)),col="black")+

  geom_text(aes(group,qD,label=sample.id),position=position_jitter(width=0.3,seed=4),size=2,color="black")
p2

p3<-hilldiv_m%>%
  filter(Order.q==2, method=="bycoverage")%>%
  ggplot(aes(group,qD, color=group))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(size=6,alpha=0.3,position=position_jitter(width=0.3,seed=4))+
  scale_color_manual(values=pal)+
  scale_x_discrete(labels=c("pea-fibrous"="pea\nfibrous","pea-minced"="pea\nminced",
                            "soy-fibrous"="soybean\nfibrous","soy-minced"="soybean\nminced"))+
  theme_classic()+
  ylab("Hill-Simpson")+
  #labs(color="group")+
  geom_text(aes(group,qD,label=sample.id),position=position_jitter(width=0.3,seed=4),size=2,color="black")
p3
plot_grid(p2,p3,nrow=2,ncol=1)
ggsave("./04-finalRplots/PlotDivByCov.tiff",width=4,height=10, units="in",dpi=2000)
ggsave("./VeggieMeat/PlotDivByCov.png", width=4,height=10, units="in",dpi=300)


#alpha diversity correlated to product properties----

facetlabs<-c('0'="richness", '1'="Hill-Shannon", '2'="Hill-Simpson")
hilldiv_m%>%
  filter(method=="bycoverage")%>%
  ggplot(aes(no.ingredients, qD))+
  geom_point(aes(col=colony_density,size=1/(daysleft+1),alpha=1/(daysleft+1)))+
  geom_text(aes(label=sample.id,col=colony_density),size=3)+
  stat_smooth(method="lm",col="lightgrey", se=F)+
  facet_grid(Order.q~.,scale="free_y",labeller=as_labeller(facetlabs))+
  scale_color_manual(values=c("high"="#A93226","medium"="#A569BD",
                              "very low/neg"="#7FB3D5"))+
  theme_bw()
ggsave("./04-finalRplots/DivVsIngredients.tiff")



hilldiv_m%>%
  filter(method=="bycoverage")%>%
  ggplot(aes(daysleft, qD))+
  geom_point(aes(col=colony_density,size=1/no.ingredients,alpha=1/no.ingredients))+
  geom_text(aes(label=sample.id,col=colony_density),size=3)+
  stat_smooth(method="lm",col="lightgrey", se=F)+
  facet_grid(Order.q~.,scale="free_y",labeller=as_labeller(facetlabs))+
  scale_color_manual(values=c("high"="#A93226","medium"="#A569BD",
                              "very low/neg"="#7FB3D5"))+
  theme_bw()
ggsave("./04-finalRplots/DivVsexpdays.tiff")


#qPCR genomic equivalents vs diversity

hilldiv_mg<-join(hilldiv_m,ge,type="left")


facetlabs<-c('0'="richness", '1'="Hill-Shannon", '2'="Hill-Simpson")
d0<-hilldiv_mg%>%
  filter(method=="bycoverage")%>%
  filter(daysleft<=100)%>%
  ggplot(aes(log10(ge_g),qD,col=group))+
  geom_point(size=5,alpha=0.5,aes(color=colony_density))+
  scale_color_manual(values=c("high"="#A93226","medium"="#A569BD",
                              "very low/neg"="#7FB3D5"))+
  theme_bw()+
  ylab("species richness")+
  labs(color="group")+
  facet_grid(Order.q~.,scale="free_y",labeller=as_labeller(facetlabs))+
  stat_smooth(method="lm",col="lightgrey", se=F)+
  geom_text(aes(log10(ge_g),qD,label=sample.id,color=colony_density),size=3)

ggsave("./04-finalRplots/DivVsgenomicequivalents.tiff")
