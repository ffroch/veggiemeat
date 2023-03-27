rare100<-readRDS("./06-outputfiles/F0_rare100.rds")
meta<-readRDS("./06-outputfiles/C_meta.rds")
#determine initial dims-----

PC <- prcomp(t(log10(otu_table(rare100[[1]]) + 1)), center=TRUE, scale=FALSE)
expl_var <- PC$sdev^2/sum(PC$sdev^2)
barplot(expl_var[1:50], ylab="EXPLAINED VARIANCE",
        names.arg=paste0("PC",seq(1:50)), col="darkgreen")
N_perm <- 1000
expl_var_perm <- matrix(NA, ncol = length(PC$sdev), nrow = N_perm)
for(k in 1:N_perm)
{
  expr_perm <- apply(otu_table(rare100[[1]]),2,sample)
  PC_perm <- prcomp(t(log10(expr_perm+1)), center=TRUE, scale=FALSE)
  expl_var_perm[k,] <- PC_perm$sdev^2/sum(PC_perm$sdev^2)
}
plot(expl_var[1:50]~seq(1:50), ylab="EXPLAINED VARIANCE",
     col="darkgreen", type='o', xlab="PRINCIPAL COMPONENTS")
lines(colMeans(expl_var_perm)[1:50]~seq(1:50),col="red")
legend("topright", c("Explained by PCS", "Explained by chance"),
       fill=c("darkgreen","red"), inset=0.02)

pval <- apply(t(expl_var_perm) >= expl_var,1,sum) / N_perm
plot(pval[1:50]~seq(1:50),col="darkred",type='o',
     xlab="PRINCIPAL COMPONENTS",ylab="PVALUE")
optPC<-head(which(pval>=0.05),1)-1
mtext(paste0("OPTIMAL NUMBER OF PRINCIPAL COMPONENTS = ",optPC))



#calculate perplexity
tdatred<-readRDS("./06-outputfiles/E4_tdatred.rds")
perp=round(nrow(tdatred)^0.5)



#plots----
metared<-join(data.frame(sample.id=rownames(tdatred)),meta,type="left")

workinglist<-list()
for (i in 1:length(dmlist)){
  set.seed(1511)
  tSNE<-Rtsne(dmlist[[i]],is_distance=T, dims=2,initial_dims=optPC, perplexity=perp,verbose=T, max_iter=1000)
  tSNEY<-data.frame(tSNE$Y)
  metared$X1<-tSNEY$X1
  metared$X2<-tSNEY$X2
  metared$method<-dmlistnames[i]
  workinglist[[i]]<-metared
}

tSNEdatacomb<-do.call("rbind", workinglist)

tSNEdatacomb<-tSNEdatacomb%>%
  mutate(method=case_when(
    method=="jaccard"~"Jaccard",
    method=="bray"~"Bray-Curtis",
    method=="jsd"~"Jensen-Shannon",
    method=="unifrac"~"UniFrac",
    method=="wunifrac"~"weighted UniFrac")
  )

#analytic plots

tSNEdatacomb%>%
  ggplot(aes(X1,X2,fill=protein.source,color=protein.source,shape=texture))+
  geom_point(size=6,alpha=0.9)+
  geom_text(aes(label=sample.id),size=3,color="black")+
  scale_color_manual(name="main protein",values=pal)+
  scale_fill_manual(name="main protein",values=pal)+
  facet_wrap(~method,ncol = 2,scale="free")+
  theme_bw()+
  scale_shape_manual(name="texture",values=c(15,16,22,21))

nicepal<-c("#AEB733", "#73832D", "#4DCF72", "#68D6C9")

#nice plots
extrapal<-c("#35978f","#01665e","#fddecb")
braytsne<-tSNEdatacomb%>%
  filter(method=="Bray-Curtis",profilelong!="")%>%
  droplevels()%>%
  ggplot(aes(X1,X2))+
  geom_point(aes(shape=group,color=group),size=8,alpha=0.7)+
  geom_text(aes(label=sample.id),size=3,color="black")+
  ggforce::geom_mark_ellipse(expand=0.05,aes(fill=profilelong),alpha=0.2,color="white")+
  scale_shape_manual(name="group", labels = c("pea-fibrous" ,"pea-minced", "soy-fibrous", "soy-minced"),
                     values=c(15,16,15,16))+
  scale_color_manual(name="group", labels = c("pea-fibrous" ,"pea-minced", "soy-fibrous", "soy-minced"),
                     values=c(nicepal[1],nicepal[2],nicepal[3],nicepal[4]))+
  scale_fill_manual(name="profile",breaks=c("Lactobacillus", "Leuconostocaceae", "Proteobacteria"), 
                    labels = c(expression(italic("Latilactobacillus")) ,expression(italic("Leuconostocaceae")), expression(italic("Proteobacteria"))),
                    values=extrapal)+
  theme_bw()+
  geom_text(aes(label="A",x=min(X1),y=max(X2)),col="black",size=6)
braytsne

#additional alpha diversity with profile as group
hilldiv_m<-readRDS("./06-outputfiles/E2_hilldiv.RDS")

alph2<-hilldiv_m%>%
  filter(Order.q==1, method=="bycoverage")%>%
  ggplot(aes(profilelong,qD,col=profilelong))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(size=3,alpha=0.5,position=position_jitter(width=0.3,seed=123))+
  scale_color_manual(name="profile",breaks=c("Lactobacillus", "Leuconostocaceae", "Proteobacteria"), 
                     labels = c(expression(italic("Latilactobacillus")) ,expression(italic("Leuconostocaceae")), expression(italic("Proteobacteria"))),
                     values=extrapal)+
  theme_classic()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        legend.position="none")+
  ylab("Hill-Shannon")+
  labs(color="profile")+
  #geom_segment(aes(x=1,y=max(qD)*0.99,xend=2,yend=max(qD)*0.99),col="black")+
  #geom_text(aes(label="*",x=1.5,y=max(qD)),col="black")+
  scale_x_discrete(labels=c("Lactobacillus"=expression(italic("Latilactobacillus")),"Leuconstocaceae"=expression(italic("Leuconostocaceae")),
                            "Proteobacteria"=expression(italic("Proteobacteria"))))+
  geom_text(aes(profilelong,qD,label=sample.id),position=position_jitter(width=0.3,seed=123),size=2,color="black")+
  geom_text(aes(label="B",x=0.6,y=max(qD)),col="black",size=6)
alph2

alph3<-hilldiv_m%>%
  filter(Order.q==2, method=="bycoverage")%>%
  ggplot(aes(profilelong,qD,col=profilelong))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(size=3,alpha=0.5,position=position_jitter(width=0.3,seed=123))+
  scale_color_manual(name="profile",breaks=c("Lactobacillus", "Leuconostocaceae", "Proteobacteria"), 
                     labels = c(expression(italic("Latilactobacillus")) ,expression(italic("Leuconostocaceae")), expression(italic("Proteobacteria"))),
                     values=extrapal)+
  scale_x_discrete(labels=c("Lactobacillus"=expression(italic("Latilactobacillus")),"Leuconostocaceae"=expression(italic("Leuconostocaceae")),
                            "Proteobacteria"=expression(italic("Proteobacteria"))))+
  theme_classic()+
  theme(legend.position="none")+
  ylab("Hill-Simpson")+
  xlab("profile")+
  labs(color="profile")+
  geom_text(aes(profilelong,qD,label=sample.id),position=position_jitter(width=0.3,seed=123),size=2,color="black")+
  geom_text(aes(label="C",x=0.6,y=max(qD)),col="black",size=6)
alph3


t1<-plot_grid(alph2,alph3,nrow=2,ncol=1)

plot_grid(braytsne,t1,nrow=1,ncol=2,rel_widths = c(1.5,1))
#save plot
ggsave("./04-finalRplots/Fig3profiles.tiff", width=9.53,height=5.96, units="in",dpi=2000)
ggsave("./VeggieMeat/Fig3profiles.png", width=9.53,height=5.96, units="in",dpi=300)



#nice plots supplements
tSNEdatacomb%>%
  filter(method%in%c("Bray-Curtis", "Jaccard", "Jensen-Shannon"),profilelong!="")%>%
  droplevels()%>%
  ggplot(aes(X1,X2))+
  geom_point(aes(shape=group,color=group),size=8,alpha=0.7)+
  geom_text(aes(label=sample.id),size=3,color="black")+
  ggforce::geom_mark_ellipse(expand=0.05,aes(fill=profilelong),alpha=0.2,color="white")+
  scale_shape_manual(name="group", labels = c("pea-fibrous" ,"pea-minced", "soy-fibrous", "soy-minced"),
                     values=c(15,16,15,16))+
  scale_color_manual(name="group", labels = c("pea-fibrous" ,"pea-minced", "soy-fibrous", "soy-minced"),
                     values=c(nicepal[1],nicepal[2],nicepal[3],nicepal[4]))+
  scale_fill_manual(name="profile",breaks=c("Lactobacillus", "Leuconostocaceae", "Proteobacteria"), labels = c("Latilactobacillus" ,"Leuconstocaceae", "Proteobacteria"),
                    values=c("#35978f","#01665e","#fddecb"))+
  facet_wrap(~method,ncol = 1,scale="free")+
  theme_bw()

#save plot
ggsave("./04-finalRplots/FigSM4.tiff", width=5.96,height=9.53, units="in",dpi=2000)
ggsave("./VeggieMeat/FigSM4.png", width=5.96,height=9.53, units="in",dpi=300)
