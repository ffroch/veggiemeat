#PERMANOVA----
#betadispersion
varlist<-c("group", "protein.source", "texture")
for (i in 1:length(dmlist)){
  for (j in 1:length(varlist)){
    k<-which(colnames(sample_data(phylo_object))==varlist[j])
    vardat<-as.vector(sample_data(phylo_object)[,k][[1]])
    betadisp<-betadisper(dmlist[[i]],vardat)
    print(dmlistnames[i])
    print(varlist[j])
    print(betadisp)
    print(anova(betadisp))
    betadispc<-betadisper(dmlist[[i]],vardat,type="centroid")
    plot(betadispc,hull=F, ellipse=T, main=paste(dmlistnames[i], varlist[j],sep=" - "))
    boxplot(betadispc, main=paste(dmlistnames[i], varlist[j],sep=" - "))
    print(anova(betadispc))
    print(TukeyHSD(betadispc))
  }
}


#adonis
for (i in 1:length(dmlist)){
  print(dmlistnames[i])
  set.seed(1234)
  print(adonis(formula=dmlist[[i]]~protein.source*texture, data=data.frame(sample_data(phylo_object)), perm=999))
  print(adonis(formula=dmlist[[i]]~texture*protein.source, data=data.frame(sample_data(phylo_object)), perm=999))
  }

for (i in 1:length(dmlist)){
  print(dmlistnames[i])
  set.seed(1234)
  print(adonis(formula=dmlist[[i]]~protein.source*texture+pseudoprod, data=data.frame(sample_data(phylo_object)), perm=999))
  print(adonis(formula=dmlist[[i]]~texture*protein.source+pseudoprod, data=data.frame(sample_data(phylo_object)), perm=999))
}


