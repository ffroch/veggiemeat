#generate NMDS plots with average distance matrices
dmlist<-readRDS("./06-outputfiles/F0_dmlist.rds")

for (i in 1:length(dmlist)){
  NMDS<-metaMDS(dmlist[[i]], k=2)
  print(stressplot(NMDS, main=dmlistnames[i]))
  
  NMDSpoints<-data.frame(NMDS$points)
  NMDSpoints$sample.id<-rownames(NMDSpoints)
  NMDS_m<-join(NMDSpoints, meta, type="left")
  
  print(
    ggplot(NMDS_m, aes(x=MDS1, y=MDS2, col=group))+
      geom_point(size = 5, alpha =0.5)+
      ggforce::geom_mark_ellipse(expand=0.02,aes(color=group,fill=group))+
      ggtitle(dmlistnames[i])+
      theme_bw()
  )
}

#NMDS with 1 random rarefied data set
phylo_object_red<-readRDS("./06-outputfiles/E4_phylo_object_red.rds")
set.seed(12345)
rarecred<-phyloseq_coverage_raref(phylo_object_red,iter=1,coverage=0.995)
rarecred<-readRDS("./06-outputfiles/E3_rarecred.rds")
oturarecred<-otu_table(rarecred)
trarecred<-t(oturarecred)

methodlist<-c("bray", "jaccard")
for (i in methodlist){
NMDS<-metaMDS(trarecred,distance=i,k=2)

print(stressplot(NMDS, main=i))

print(plot(NMDS, main=i))

dscores<-as.data.frame(scores(NMDS, display="sites"))
dscores$sample.id<-rownames(trarecred)


scores_meta<-join(dscores,meta,by="sample.id",type="left")

print(
  ggplot(scores_meta, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes( shape = texture, colour = protein.source))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Time", y = "NMDS2", shape = "Type")  + 
  scale_colour_manual(values = c("#009E73", "#E69F00"))+
  ggtitle(i)
)
met<-data.frame(sample.id=sample_names(phylo_object_red))
met<-join(met, data.frame(sample.id=meta$sample.id,group=meta$group), by="sample.id",type="left")
groups<-data.frame(group=met$group)
set.seed=1234
metfit<-envfit(NMDS,groups,permutations = 999, na.rm = T)

#metfit_coord_cont = as.data.frame(scores(metfit, "vectors")) * ordiArrowMul(metfit)
metfit_coord_cat = as.data.frame(scores(metfit, "factors")) 
print(
  ggplot(data = scores_meta, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = scores_meta, aes(colour = protein.source,shape=texture), size = 3, alpha = 0.5) + 
  scale_colour_manual(values = c("orange", "steelblue")) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) +
  labs(colour = "main protein")
)

}

#generate NMDS plots with average distance matrices for supplements
NMDS<-metaMDS(dmlist[[1]], k=2)
df1<-data.frame(diss=NMDS$diss,dist=NMDS$dist,dhat=NMDS$dhat,method="Bray-Curtis")
NMDS<-metaMDS(dmlist[[2]], k=2)
df2<-data.frame(diss=NMDS$diss,dist=NMDS$dist,dhat=NMDS$dhat,method="weighted UniFrac")
NMDS<-metaMDS(dmlist[[3]], k=2)
df3<-data.frame(diss=NMDS$diss,dist=NMDS$dist,dhat=NMDS$dhat,method="UniFrac")
NMDS<-metaMDS(dmlist[[4]], k=2)
df4<-data.frame(diss=NMDS$diss,dist=NMDS$dist,dhat=NMDS$dhat,method="Jaccard")
NMDS<-metaMDS(dmlist[[5]], k=2)
df5<-data.frame(diss=NMDS$diss,dist=NMDS$dist,dhat=NMDS$dhat,method="Jensen-Shannon")

df<-rbind(df1, df2, df3, df4, df5)

pal<-c("#AEB733", "#73832D", "#4DCF72", "#68D6C9")

p1<-ggplot(df, aes(diss,dist))+
  geom_point(color="blue", alpha=0.5)+
  geom_line(aes(diss, dhat),size=1.5,color="red")+
  theme_bw()+
  xlab("Observed Dissimilarity")+
  ylab("Ordination Distance")+
  facet_wrap(~method, ncol=1, scales="free_x")
  
p1

nmdslist<-list()
for (i in 1:length(dmlist)){
  NMDS<-metaMDS(dmlist[[i]], k=2)
  NMDSpoints<-data.frame(NMDS$points)
  NMDSpoints$sample.id<-rownames(NMDSpoints)
  NMDS_m<-join(NMDSpoints, meta, type="left")
  NMDS_m$method<-dmlistnames[i]
  nmdslist[[i]]<-NMDS_m
}
nmdstotal<-do.call(rbind.data.frame, nmdslist)
 
supp.labs <- c("Bray-Curtis", "Jaccard", "Jensen-Shannon", "UniFrac", "weighted UniFrac")
names(supp.labs) <- c("bray","jaccard", "jsd", "unifrac", "wunifrac")

p2<-ggplot(nmdstotal, aes(x=MDS1, y=MDS2, col=group))+
      geom_point(size = 5, alpha =0.5)+
      ggforce::geom_mark_ellipse(expand=0.02,aes(color=group,fill=group))+
      theme_bw()+
  scale_color_manual(values=pal)+
  scale_fill_manual(values=pal)+
  facet_wrap(~method, ncol=1, scales="free", labeller=labeller(method=supp.labs))

plot_grid(p1, p2, rel_widths = c(1,1.7))

ggsave("./04-finalRplots/NMDSplots.tiff",width=8.9,height=9.8, units="in",dpi=2000)
ggsave("./VeggieMeat/NMDSplots.png", width=8.9,height=9.8, units="in",dpi=300)
