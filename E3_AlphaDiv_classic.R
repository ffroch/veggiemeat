# alphadiversity of not rarefied data
divunrar<-estimate_richness(phylo_object, measures=c("Observed", "Shannon", "Simpson"))
divunrar
plot_richness(phylo_object,x="group", measures=c("Observed", "Shannon", "Simpson"))

#comparison rarefaction methods and sizes
divunrar<-estimate_richness(phylo_object, measures=c("Observed", "Shannon", "Simpson"))
divunrar$meth<-"unrar"
divunrar$sample.id<-rownames(divunrar)
set.seed(12345)
rarec<-phyloseq_coverage_raref(phylo_object,iter=1,coverage=0.995)
divc<-estimate_richness(rarec, measures=c("Observed", "Shannon", "Simpson"))
divc$meth<-"cov"
divc$sample.id<-rownames(divc)
raresize439<-rarefy_even_depth(phylo_object, sample.size=439,rngseed=123)
divsize439<-estimate_richness(raresize439, measures=c("Observed", "Shannon", "Simpson"))
divsize439$meth<-"size439"
divsize439$sample.id<-rownames(divsize439)
raresize2463<-rarefy_even_depth(phylo_object, sample.size=2463,rngseed=123)
divsize2463<-estimate_richness(raresize2463, measures=c("Observed", "Shannon", "Simpson"))
divsize2463$meth<-"size2463"
divsize2463$sample.id<-rownames(divsize2463)

methcomp<-rbind(divunrar, divc, divsize439, divsize2463)
methcomplong<-data.frame(pivot_longer(methcomp, 1:3, names_to="index", values_to="values"))

ggplot(methcomplong, aes(x=sample.id, y=values, col=meth ))+
  geom_jitter(width=0.3, alpha=0.5)+
  facet_grid(index~., scales="free")

methcomplong %>%
  filter(index=="Observed")%>%
  ggplot(aes(x=sample.id, y=values, col=meth ))+
  geom_point()

pdunrar<-estimate_pd(phylo_object)
pdunrar$meth<-"unrar"
pdunrar$sample.id<-rownames(pdunrar)
pdc<-estimate_pd(rarec)
pdc$meth<-"cov"
pdc$sample.id<-rownames(pdc)
pd439<-estimate_pd(raresize439)
pd439$meth<-"size439"
pd439$sample.id<-rownames(pd439)
pd2463<-estimate_pd(raresize2463)
pd2463$meth<-"size2463"
pd2463$sample.id<-rownames(pd2463)

pdbind<-rbind(pdunrar, pdc, pd439, pd2463)
ggplot(pdbind, aes(x=sample.id, y=PD, col=meth ))+
  geom_jitter(width=0.3, alpha=0.5)


pdbindrestruct<-data.frame(meth=pdbind$meth,sample.id=pdbind$sample.id, index="Faith's PD", values=pdbind$PD)
combineall<-rbind(methcomplong, pdbindrestruct)

ggplot(combineall, aes(x=sample.id, y=values, col=meth ))+
  geom_jitter(width=0.3, alpha=0.5)+
  facet_wrap(~index, scales="free", ncol=2)

combineall_m<-join(combineall,meta,type="left")
combineall_m%>%
  filter(meth=="cov")%>%
  ggplot(aes(x=group, y=values))+
  geom_boxplot()+
  geom_jitter(size=3,alpha=0.5,position=position_jitter(width=0.3,seed=4))+
  facet_wrap(~index, scales="free", ncol=2)+
  theme_bw()
