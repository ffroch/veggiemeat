#generate NMDS plots with average distance matrices
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
set.seed(12345)
rarec<-phyloseq_coverage_raref(phylo_object,iter=1,coverage=0.995)
oturarec<-otu_table(rarec)
trarec<-t(oturarec)

methodlist<-c("bray", "jaccard")
for (i in methodlist){
NMDS<-metaMDS(trarec,distance=i,k=2)

print(stressplot(NMDS, main=i))

print(plot(NMDS, main=i))

dscores<-as.data.frame(scores(NMDS, display="sites"))
dscores$sample.id<-rownames(trarec)


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
met<-data.frame(sample.id=sample_names(phylo_object))
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




