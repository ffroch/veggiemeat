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
perp=round(nrow(tdat)^0.5)



#plots----
metared<-join(data.frame(sample.id=rownames(tdat)),meta,type="left")

workinglist<-list()
for (i in 1:length(dmlist)){
  set.seed(12345)
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

pal<-RColorBrewer::brewer.pal(4,"Paired")[c(1,3)]
tSNEdatacomb%>%
  ggplot(aes(X1,X2,fill=protein.source,color=protein.source,shape=texture))+
  geom_point(size=6,alpha=0.9)+
  geom_text(aes(label=sample.id),size=3,color="black")+
  scale_color_manual(name="main protein",values=pal)+
  scale_fill_manual(name="main protein",values=pal)+
  facet_wrap(~method,ncol = 2,scale="free")+
  theme_bw()+
  scale_shape_manual(name="texture",values=c(15,16,22,21))
