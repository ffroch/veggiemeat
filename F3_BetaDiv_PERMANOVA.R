#PERMANOVA----
#betadispersion
varlist<-c("group", "protein.source", "texture")
for (i in 1:length(dmlist)){
  for (j in 1:length(varlist)){
    k<-which(colnames(sample_data(relab_po))==varlist[j])
    vardat<-as.vector(sample_data(relab_po)[,k][[1]])
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
varexpl<-data.frame("method"=NA, "R2"=NA)
for (i in 1:length(dmlist)){
  print(dmlistnames[i])
  set.seed(1234)
  print(adonis2(formula=dmlist[[i]]~protein.source*texture, data=data.frame(sample_data(relab_po)), perm=999))
  print(adonis2(formula=dmlist[[i]]~texture*protein.source, data=data.frame(sample_data(relab_po)), perm=999))
  varexpl[i,1]<-dmlistnames[i]
  varexpl[i,2]<-(1-adonis2(formula=dmlist[[i]]~protein.source*texture, data=data.frame(sample_data(relab_po)), perm=999)$R2[4])*100
  }


for (i in 1:length(dmlist)){
  print(dmlistnames[i])
  set.seed(1234)
  print(adonis2(formula=dmlist[[i]]~protein.source*texture+pseudoprod, data=data.frame(sample_data(relab_po)), perm=999))
  print(adonis2(formula=dmlist[[i]]~texture*protein.source+pseudoprod, data=data.frame(sample_data(relab_po)), perm=999))
}




#for supplements
dir.create("./06-outputfiles/permanovatables/")

dmlistshort<-list(bd, jd, jsd)
dmlistnamesshort<-c("bray", "jaccard", "jsd")
permanova<-read_docx()
permanova%>%
  print(target="./06-outputfiles/permanovatables/permanova.docx")

for (i in 1:length(dmlistshort)){
  headline<-attributes(adonis2(formula=dmlistshort[[i]]~protein.source*texture, data=data.frame(sample_data(relab_po)), perm=999))$heading
  set.seed(1234)
  flextable(adonis2(formula=dmlistshort[[i]]~protein.source*texture, data=data.frame(sample_data(relab_po)), perm=999))%>%
    add_header_row(values=paste(dmlistnamesshort[i],headline[1], headline[2], sep="\n"), colwidths=5)%>%
    set_table_properties(layout = "autofit")%>%
    align(align="left", part="all")%>%
    save_as_docx(path = "./06-outputfiles/permanovatables/temp1.docx")
  headline<-attributes(adonis2(formula=dmlistshort[[i]]~protein.source*texture, data=data.frame(sample_data(relab_po)), perm=999))$heading
  set.seed(1234)
  flextable(adonis2(formula=dmlistshort[[i]]~texture*protein.source, data=data.frame(sample_data(relab_po)), perm=999))%>%
    add_header_row(values=paste(dmlistnamesshort[i],headline[1], headline[2], sep="\n"), colwidths=5)%>%
    set_table_properties(layout = "autofit")%>%
    align(align="left", part="all")%>%
    save_as_docx(path = "./06-outputfiles/permanovatables/temp2.docx")
  permanova<-read_docx("./06-outputfiles/permanovatables/permanova.docx")
  permanova<-permanova%>%
    body_add_docx(src="./06-outputfiles/permanovatables/temp1.docx")%>%
    body_add_docx(src="./06-outputfiles/permanovatables/temp2.docx")%>%
    print(target = "./06-outputfiles/permanovatables/permanova.docx")
  
}
#with short list
varexpl<-data.frame("method"=NA, "R2"=NA)
for (i in 1:length(dmlistshort)){
  print(dmlistnamesshort[i])
  set.seed(1234)
  print(adonis2(formula=dmlistshort[[i]]~protein.source*texture, data=data.frame(sample_data(relab_po)), perm=999))
  print(adonis2(formula=dmlistshort[[i]]~texture*protein.source, data=data.frame(sample_data(relab_po)), perm=999))
  varexpl[i,1]<-dmlistnamesshort[i]
  varexpl[i,2]<-(1-adonis2(formula=dmlistshort[[i]]~protein.source*texture, data=data.frame(sample_data(relab_po)), perm=999)$R2[4])*100
}

saveRDS(varexpl,"./06-outputfiles/permanovaR2.rds")
