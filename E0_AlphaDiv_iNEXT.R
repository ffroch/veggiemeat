cov<-iNEXT(data[,2:ncol(data)])
cov

tiff("./05-analysisRplots/01rarefaction-diversity_cov-reads.tiff")
ggiNEXT(cov, type=2)+
  scale_shape_manual(values=c(rep(1,length(reads$reads))))+
  xlim(0,median)+
  theme_bw()+
  guides(shape="none")
dev.off()
tiff("./05-analysisRplots/01rarefaction-diversity_div-coverage.tiff")
ggiNEXT(cov, type=3)+
  scale_shape_manual(values=c(rep(1,length(reads$reads))))+
  xlim(0.9,1)+
  theme_bw()+
  guides(shape="none")
dev.off()


#rarefaction curves with minimum sized sample, second minimum sized sample, 10000 and median
covsize<-iNEXT(data[,2:ncol(data)],size=c(min1,min2,10000,median))
covsize


tiff("./05-analysisRplots/01rarefaction-sizerarefied-diversity_cov-reads.tiff")
ggiNEXT(covsize, type=2)+
  scale_shape_manual(values=c(rep(1,length(reads$reads))))+
  xlim(0,median)+
  theme_bw()+
  guides(shape="none")
dev.off()
tiff("./05-analysisRplots/01rarefaction-sizerarefied-diversity_div-coverage.tiff")
ggiNEXT(covsize, type=3)+
  scale_shape_manual(values=c(rep(1,length(reads$reads))))+
  xlim(0.9,1)+
  theme_bw()+
  guides(shape="none")
dev.off()
