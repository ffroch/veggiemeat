#hand made LEfse plots
combt<-readRDS("./06-outputfiles/G2_combt.rds")
combtnames<-data.frame("sample.id"=rownames(combt))
combtn<-cbind(combtnames, combt)


combtnlong<-pivot_longer(combtn, 3:ncol(combtn),names_to="feature",values_to="relab")
combtnlong$newfeature<-sub(".*\\.","",combtnlong$feature)


lefseresults<-read.table("./07-lefse/lefse_data.res",header=F,sep="\t")
colnames(lefseresults)<-c("feature", "LDA","lefsegroup","adjLDA","pval")
lefsefiltered<-lefseresults[!is.na(lefseresults$adjLDA),]

longjoined<-join(combtnlong,lefsefiltered,by="feature", type="left")
longfiltered<-longjoined[!is.na(longjoined$lefsegroup),]

longfiltered$newfeature<-gsub("Firmicutes", "Bacillota", longfiltered$newfeature)
longfiltered$newfeature<-gsub("Proteobacteria", "Pseudomonadota", longfiltered$newfeature)


nicepal<-c("#AEB733", "#73832D", "#4DCF72", "#68D6C9")

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
  scale_color_manual(values=nicepal)+
    scale_y_discrete(labels = function(x) {
    italic_labels <- c("Leuconostoc", "Lactobacillaceae", "Enterococcaceae", "Enterococcus", "Pseudomonas", "Pseudomonadaceae", "Xanthomonadaceae", "Sphingomonadaceae", "Sphingomonas", "Photobacterium", "Dellaglioa", "Flavobacteriaceae") # Define the features you want in italics
    lbls <- ifelse(x %in% italic_labels, parse(text = paste0("italic(", x, ")")), x)
    return(lbls)
  })


p2<-longfiltered%>%
  mutate(neworderfeature=reorder(newfeature,adjLDA))%>%
  ggplot(aes(adjLDA,neworderfeature,fill=lefsegroup))+
  geom_bar(position="dodge", stat="identity")+
  facet_grid(lefsegroup~.,scale="free",space="free")+
  theme_bw()+
  theme(axis.text.y=element_blank(),axis.title.y=element_blank(),axis.ticks.y=element_blank(),
        plot.margin=unit(c(0.5,0.5,0.5,0),"null"))+
  labs(x="LDA score (log 10)",fill="group")+
  scale_fill_manual(values=nicepal)



plot_grid(p1,p2,rel_widths=c(3,2),nrow=1)
ggsave("./04-finalRplots/lefsestatplotnewnom20240222.tiff", width=7,height=8, units="in",dpi=2000)

ggsave("VeggieMeat/lefsestatplotnewnom20240222.png", width=7,height=8, units="in",dpi=300)



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
  scale_color_manual(values=nicepal)+
    scale_y_discrete(labels = function(x) {
    italic_labels <- c("Leuconostoc", "Lactobacillaceae", "Enterococcaceae", "Enterococcus", "Pseudomonas", "Pseudomonadaceae", "Xanthomonadaceae", "Sphingomonadaceae", "Sphingomonas", "Photobacterium", "Dellaglioa", "Flavobacteriaceae") # Define the features you want in italics
    lbls <- ifelse(x %in% italic_labels, parse(text = paste0("italic(", x, ")")), x)
    return(lbls)
  })


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
  scale_color_manual(values=nicepal)+
    scale_y_discrete(labels = function(x) {
    italic_labels <- c("Leuconostoc", "Lactobacillaceae", "Enterococcaceae", "Enterococcus", "Pseudomonas", "Pseudomonadaceae", "Xanthomonadaceae", "Sphingomonadaceae", "Sphingomonas", "Photobacterium", "Dellaglioa", "Flavobacteriaceae") # Define the features you want in italics
    lbls <- ifelse(x %in% italic_labels, parse(text = paste0("italic(", x, ")")), x)
    return(lbls)
  })

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
ggsave("./04-finalRplots/lefsestatplotdiffscalesnewnom20240222.tiff", width=7,height=8, units="in",dpi=2000)

ggsave("VeggieMeat/lefsestatplotdiffscalesnewnom20240222.png", width=7,height=8, units="in",dpi=300)
