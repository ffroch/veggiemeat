#import single distance matrices
bd<-readRDS("./06-outputfiles/F0_bd.rds")
jd<-readRDS("./06-outputfiles/F0_jd.rds")
jsd<-readRDS("./06-outputfiles/F0_jsd.rds")

#import relative abundance data
relab_po<-readRDS("./06-outputfiles/E4_relab_po.rds")

dmlistshort<-list(bd, jd, jsd)
dmlistnamesshort<-c("Bray-Curtis", "Jaccard", "Jensen-Shannon")

dflist<-list()
j=1
for (i in 1:length(dmlistshort)){
  dflist[[j]]<-data.frame(adonis2(formula=dmlistshort[[i]]~protein.source*texture, data=data.frame(sample_data(relab_po)), perm=999)[,1:5])
  dflist[[j+1]]<-data.frame(adonis2(formula=dmlistshort[[i]]~texture*protein.source, data=data.frame(sample_data(relab_po)), perm=999)[,1:5])
  j=j+2
  }

dflistdf<-do.call(rbind.data.frame, dflist)
dflistdf<-format(dflistdf, digits=3)



dflistdfnames<-c("protein source", "texture", "protein source:texture", "residual", "total")

dflistdf<-data.frame(names=dflistdfnames, dflistdf)
row.names(dflistdf)<-NULL
dflistdf[,5]<-stringr::str_replace_all(dflistdf[,5], "NA", "")
dflistdf[,6]<-stringr::str_replace_all(dflistdf[,6], "NA", "")

headerlist<-list()
j=1
for (i in 1:length(dmlistshort)){
  headerlist[[j]]<-attributes(adonis2(formula=dmlistshort[[i]]~protein.source*texture, data=data.frame(sample_data(relab_po)), perm=999))$heading
  headerlist[[j+1]]<-attributes(adonis2(formula=dmlistshort[[i]]~texture*protein.source, data=data.frame(sample_data(relab_po)), perm=999))$heading
  j=j+2
}

headerlistnew<-list()
for (i in 1:length(headerlist)){
headerlistnew[[i]]<-paste0(headerlist[[i]][1],headerlist[[i]][2])
}

headerdf<-do.call(rbind.data.frame, headerlistnew)
headerdf<-data.frame(header=headerdf[,1])
headerdf<-str_replace_all(headerdf[,1], ", data =", ",\ndata =")

kableExtra::kbl(dflistdf[,1:6], 
                caption = "\\label{tab1}Caption centered above table", longtable=T, booktabs=T, escape=F,valign="b",
                align="l",
                col.names=c("", "Df", "Sum of Sqs", "R2", "F", "Pr. F"))%>%
  kable_styling(latex_options=c("repeat_header")) %>%
  pack_rows("Bray-Curtis", 1,10, bold=T)%>%
  pack_rows("Jaccard", 11,20, bold=T)%>%
  pack_rows("Jensen-Shannon", 21,30, bold=T)%>%
  pack_rows(headerdf[1], 1,5, bold=F, italic=T)%>%
  pack_rows(headerdf[2], 6,10, bold=F, italic=T)%>%
  pack_rows(headerdf[3], 11,15, bold=F, italic=T)%>%
  pack_rows(headerdf[4], 16,20, bold=F, italic=T)%>%
  pack_rows(headerdf[5], 21,25, bold=F, italic=T)%>%
  pack_rows(headerdf[6], 26,30, bold=F, italic=T)
  
saveRDS(dflistdf, "./06-outputfiles/M1_smtabpermanova.rds")
saveRDS(headerdf,"./06-outputfiles/M1_smtabpermanovaheader.rds")
