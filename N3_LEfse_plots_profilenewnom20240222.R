#hand made LEfse plots
combt<-readRDS("./06-outputfiles/N2_combt_profile.rds")
combtnames<-data.frame("sample.id"=rownames(combt))
combtn<-cbind(combtnames, combt)

combtn$profilelong<-str_replace_all(combtn$profilelong, "Lactobacillus", "Latilactobacillus")
combtn$profilelong<-str_replace_all(combtn$profilelong, "Proteobacteria", "Pseudomonadota")

combtnlong<-pivot_longer(combtn, 3:ncol(combtn),names_to="feature",values_to="relab")
combtnlong$newfeature<-sub(".*\\.","",combtnlong$feature)

combtnlong$newfeature<-gsub("Firmicutes", "Bacillota", combtnlong$newfeature)
combtnlong$newfeature<-gsub("Proteobacteria", "Pseudomonadota", combtnlong$newfeature)

lefseresults<-read.table("./07-lefse/lefse_data_profile.res",header=F,sep="\t")
colnames(lefseresults)<-c("feature", "LDA","lefsegroup","adjLDA","pval")
lefsefiltered<-lefseresults[!is.na(lefseresults$adjLDA),]

longjoined<-join(combtnlong,lefsefiltered,by="feature", type="left")
longfiltered<-longjoined[!is.na(longjoined$lefsegroup),]
longfiltered$lefsegroup<-str_replace_all(longfiltered$lefsegroup,"Lactobacillus", "Latilactobacillus")
longfiltered$lefsegroup<-str_replace_all(longfiltered$lefsegroup,"Proteobacteria", "Pseudomonadota")

nicepal<-c(brewer.pal(n=11,name="RdBu")[11],brewer.pal(n=11,name="PRGn")[9],brewer.pal(n=11,name="BrBG")[3])

p1<-longfiltered%>%
  mutate(neworderfeature=reorder(newfeature,adjLDA))%>%
  ggplot(aes(relab,neworderfeature,col=profilelong))+
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
    italic_labels <- c("Lactobacillus", "Vibrionaceae", "Shewanella", "Shewanellaceae", "Moraxellaceae", "Beijerinckiaceae", "Leuconostoc", "Lactobacillaceae", "Enterococcaceae", "Enterococcus", "Pseudomonas", "Pseudomonadaceae", "Xanthomonadaceae", "Sphingomonadaceae", "Sphingomonas", "Photobacterium", "Dellaglioa", "Flavobacteriaceae", "Psychrobacter", "Acinetobacter", "Weeksellaceae", "Comamonadaceae") # Define the features you want in italics
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
  labs(x="LDA score (log 10)",fill="profile")+
  scale_fill_manual(values=nicepal, breaks=c("Latilactobacillus", "Leuconostocaceae", "Pseudomonadota"), 
                    labels = c(expression(italic("Latilactobacillus")) ,expression(italic("Leuconostocaceae")), "Pseudomonadota"))



plot_grid(p1,p2,rel_widths=c(3,2),nrow=1)
ggsave("./04-finalRplots/lefsestatplot_profilenewnom20240222.tiff", width=7,height=8, units="in",dpi=2000)

ggsave("VeggieMeat/lefsestatplot_profilenewnom20240222.png", width=7,height=8, units="in",dpi=300)

