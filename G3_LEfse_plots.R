#hand made LEfse plots
combtnames<-data.frame("sample.id"=rownames(combt))
combtn<-cbind(combtnames, combt)



combtnlong<-pivot_longer(combtn, 3:ncol(combtn),names_to="feature",values_to="relab")
combtnlong$newfeature<-sub(".*\\.","",combtnlong$feature)


lefseresults<-read.table("./07-lefse/lefse_data.res",header=F,sep="\t")
colnames(lefseresults)<-c("feature", "LDA","lefsegroup","adjLDA","pval")
lefsefiltered<-lefseresults[!is.na(lefseresults$adjLDA),]

longjoined<-join(combtnlong,lefsefiltered,by="feature", type="left")
longfiltered<-longjoined[!is.na(longjoined$lefsegroup),]

p1<-longfiltered%>%
  mutate(neworderfeature=reorder(newfeature,adjLDA))%>%
  ggplot(aes(relab,neworderfeature,col=group))+
  #geom_jitter()+
  stat_summary(fun.data=median_hilow,geom="pointrange",
               fun.args=list(conf.int=0.5),
               position=position_dodge(width=0.8),
               size=0.25)+
  facet_grid(lefsegroup~.,scale="free",space="free")+
  theme_bw()+
  theme(legend.position = "none",strip.text.y=element_blank(),
        plot.margin=unit(c(0.5,0.00001,0.5,0.02),"null"),
        axis.title.y= element_blank())+
  xlab("relative abundance")+
  scale_color_brewer(palette="Paired")

p2<-longfiltered%>%
  mutate(neworderfeature=reorder(newfeature,adjLDA))%>%
  ggplot(aes(adjLDA,neworderfeature,fill=lefsegroup))+
  geom_bar(position="dodge", stat="identity")+
  facet_grid(lefsegroup~.,scale="free",space="free")+
  theme_bw()+
  theme(axis.text.y=element_blank(),axis.title.y=element_blank(),axis.ticks.y=element_blank(),
        plot.margin=unit(c(0.5,0.5,0.5,0),"null"))+
  labs(x="LDA score (log 10)",fill="group")+
  scale_fill_brewer(palette="Paired")


plot_grid(p1,p2,rel_widths=c(3,2),nrow=1)
ggsave("./04-finalRplots/lefsestatplot.tiff")

nicepal<-RColorBrewer::brewer.pal(4,"Paired")

p3<-longfiltered%>%
  filter(lefsegroup=="pea-fibrous"| lefsegroup=="pea-minced")%>%
  mutate(neworderfeature=reorder(newfeature,adjLDA))%>%
  ggplot(aes(relab,neworderfeature,col=group))+
  stat_summary(fun.data=median_hilow,geom="pointrange",
               fun.args=list(conf.int=0.5),
               position=position_dodge(width=0.8),
               size=0.25)+
  facet_grid(lefsegroup~.,scale="free",space="free")+
  theme_bw()+
  theme(legend.position = "none",strip.text.y=element_blank(),axis.title.x=element_blank(),
        plot.margin=unit(c(0.5,0.01,0.5,0.02),"null"),
        axis.title.y=element_blank())+
  scale_color_manual(values=nicepal)


p4<-longfiltered%>%
  filter(lefsegroup=="pea-fibrous"| lefsegroup=="pea-minced")%>%
  mutate(neworderfeature=reorder(newfeature,adjLDA))%>%
  ggplot(aes(adjLDA,neworderfeature,fill=lefsegroup))+
  geom_bar(position="dodge", stat="identity")+
  facet_grid(lefsegroup~.,scale="free",space="free")+
  theme_bw()+
  theme(axis.text.y=element_blank(),axis.title.y=element_blank(),axis.ticks.y=element_blank(),
        legend.position = "none",axis.title.x=element_blank(),
        plot.margin=unit(c(0.5,0.5,0.5,0),"null"))+
  scale_fill_manual(values=nicepal[1:2])

p5<-longfiltered%>%
  filter(lefsegroup=="soy-fibrous"| lefsegroup=="soy-minced")%>%
  mutate(neworderfeature=reorder(newfeature,adjLDA))%>%
  ggplot(aes(relab,neworderfeature,col=group))+
  stat_summary(fun.data=median_hilow,geom="pointrange",
               fun.args=list(conf.int=0.5),
               position=position_dodge(width=0.8),
               size=0.25)+
  facet_grid(lefsegroup~.,scale="free",space="free")+
  theme_bw()+
  theme(legend.position = "none",strip.text.y=element_blank(),
        plot.margin=unit(c(0.5,0.001,0.5,0.004),"null"),
        axis.title.y=element_blank())+
  xlab("relative abundance")+
  scale_color_manual(values=nicepal)

p6<-longfiltered%>%
  filter(lefsegroup=="soy-fibrous"| lefsegroup=="soy-minced")%>%
  mutate(neworderfeature=reorder(newfeature,adjLDA))%>%
  ggplot(aes(adjLDA,neworderfeature,fill=lefsegroup))+
  geom_bar(position="dodge", stat="identity")+
  facet_grid(lefsegroup~.,scale="free",space="free")+
  theme_bw()+
  theme(axis.text.y=element_blank(),axis.title.y=element_blank(),axis.ticks.y=element_blank(),legend.position = "none",
        plot.margin=unit(c(0.5,0.5,0.5,0),"null"))+
  scale_fill_manual(values=nicepal[3:4])+
  labs(x="LDA score (log 10)",fill="group")

plot_grid(p3,p4,p5,p6,rel_widths=c(3,1),nrow=2)
ggsave("./04-finalRplots/lefsestatplotdiffscales.tiff")
