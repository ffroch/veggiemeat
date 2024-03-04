#use LEfse prepared data for group comparison including correction for muliltple comparisons

lefsedat<-read.table("./07-lefse/lefse_data_profile.txt",header=T,sep="\t", row.names = 1)

#kruskal for all
combt<-data.frame(t(lefsedat))
combt[,2:ncol(combt)]<-lapply(combt[,2:ncol(combt)], function(x) as.numeric(as.character(x)))


sigpval<-data.frame("Var"=names(combt)[1],"KW"="KW","BH"= "BH","BY"= "BY")

pvaldf<-data.frame("Var"= "Var","feature"= "feature","pval"= "pval","adjpval.BH"="adjpval.BH","adjpval.BY"="adjpval.BY")
for (i in 2:ncol(combt)){
  kt<-kruskal.test(combt[,i]~combt[,1],)
  pvaldf[i-1,1]<-colnames(combt)[1]
  pvaldf[i-1,2]<-colnames(combt)[i]
  pvaldf[i-1,3]<-kt$p.value
  pvaldf[i-1,4]<-p.adjust(pvaldf[i-1,3],method="BH", n=nrow(pvaldf))
  pvaldf[i-1,5]<-p.adjust(pvaldf[i-1,3],method="BY", n=nrow(pvaldf))
}
sigpval[1,2]<-sum(pvaldf[,3]<=0.05,na.rm=T)
sigpval[1,3]<-sum(pvaldf[,4]<=0.05,na.rm=T)
sigpval[1,4]<-sum(pvaldf[,5]<=0.05,na.rm=T)

print(sigpval)

filtered_KW<-pvaldf[pvaldf$pval<0.05,]
filtered_BH<-pvaldf[pvaldf$adjpval.BH<0.05,]
filtered_BY<-pvaldf[pvaldf$adjpval.BY<0.05,]

print(filtered_KW)
print(filtered_BH)
print(filtered_BY)


saveRDS(combt,"./06-outputfiles/N2_combt_profile.rds")
