
unfil<-read.table("D:/0300LMM/0302VeggieProducts/05Miseq/02-qiime/table.withoutnontargetmitochondriaandchloroplasts-summary/sample-frequency-detail.csv",sep=",")
fil<-read.table("D:/0300LMM/0302VeggieProducts/07MiSeqFerdl/02-qiime/table-no-mitochondria-no-chlorplast/c2dc9899-85c4-4fa0-9976-2914de069dbf/data/sample-frequency-detail.csv",sep=",")

comb<-join(unfil,fil,by="V1",type="left")
colnames(comb)<-c("sample.id","reads_w_chloro","reads_wo_chloro")
comb$leftper<-100/comb$reads_w_chloro*comb$reads_wo_chloro


qpcr<-read.table("D:/0300LMM/0302VeggieProducts/06Rnew/00_MetaData/genomicequivalentsperg.txt",header=T)
meta<-read.table("D:/0300LMM/0302VeggieProducts/06Rnew/00_MetaData/metacat.txt", header=T, sep = "\t") 

comb_qpcr<-join(qpcr, comb,by="sample.id", type="left")
comb_qpcr_m<-join(comb_qpcr,meta,type="left")

ggplot(comb_qpcr_m,aes(leftper,log10(ge_g),col=colony_density))+
    geom_point(size=5, alpha=0.5)+
  theme_bw()

# correct qPCR results by the percentage of reads without chloroplasts and mitochondria
comb_qpcr_m$corr_bce<-comb_qpcr_m$ge_g/100*comb_qpcr_m$leftper
# get a better name for comb_qpcr_m
qpcr <- comb_qpcr_m
# if qpcr$corr_bce is NA, replace it by qpcr$ge_g
qpcr$corr_bce[is.na(qpcr$corr_bce)]<-qpcr$ge_g[is.na(qpcr$corr_bce)]
# max value in the corrected qPCR results
log10(max(qpcr$corr_bce,na.rm=T))
# make it long format for the columns ge_g and corr_bce
long<-data.frame(pivot_longer(qpcr,c(ge_g,corr_bce)))
# replace ge_g by genomic equivalents and corr_bce by genomic equivalents corrected
long$name<-factor(long$name,levels=c("ge_g","corr_bce"),labels=c("genomic equivalents","genomic equivalents corrected"))


paldens<-c("high"="#7A2048","medium"="#408EC6","very low/neg"="#1E2761")
p1 <- ggplot(long,aes(colony_density,log10(value)))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(size=4,alpha=0.5,aes(color=colony_density),position=position_jitter(width=0.25,seed=2))+
  geom_text(size=2.5,aes(label=ID),position=position_jitter(width=0.25,seed=2))+
  scale_color_manual(values=paldens)+
  theme_classic()+
  ylab("log10(genomic equivalents per g food)")+
  xlab("colony density on blood agar")+
  labs(color="colony density\non blood agar") +
  facet_wrap(~name)

p1
ggsave(paste0(outputdir,"S1_colonydensge_g.tiff"), width=20,height=15,units="cm" )

extrapal<-c(brewer.pal(n=11,name="RdBu")[11],brewer.pal(n=11,name="PRGn")[9],brewer.pal(n=11,name="BrBG")[3], "lightgrey")
# if profilelong is empty, replace it by "no 16S data"
qpcr$profilelong[qpcr$profilelong==""]<-"no 16S data"
# replace cut by fibrous in category
qpcr$category[qpcr$category=="cut"]<-"fibrous"
# make groups based on main_prot and category
qpcr$group <- paste0(qpcr$main_prot, "-", qpcr$category)



qpcr$profilelong[is.na(qpcr$profilelong)]<-"no 16S data"
p2<-ggplot(qpcr,aes(group,log10(corr_bce)))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(size=4,alpha=0.5,aes(color=profilelong),position=position_jitter(width=0.25,seed=2))+
  geom_text(size=2.5,aes(label=ID),position=position_jitter(width=0.25,seed=2))+
  scale_color_manual(name="profile",breaks=c("Lactobacillus", "Leuconostocaceae", "Proteobacteria", "no 16S data"), 
                    labels = c(expression(italic("Latilactobacillus")) ,expression(italic("Leuconostocaceae")), expression(italic("Pseudomonadota")), "no 16S data"),
                    values=extrapal)+
  theme_classic()+
  ylab("log10(genomic equivalents per g food)")+
  xlab("group")+
  labs(color="tSNE-cluster")

p2
ggsave(paste0(outputdir,"S2_plotge_g.tiff"), width=20,height=15,units="cm" )


#how good do producers know der products
qpcr$centered_cooking_time_clustered<-factor(qpcr$centered_cooking_time_clustered,
                                             levels=c("<3 min", "3-5 min", "5-8 min", ">8 min", "no rec."))
#replace thorougly by thoroughly
qpcr$additional_recommendation[qpcr$additional_recommendation=="consume only thorougly heated"]<-"consume only thoroughly heated"                                            

                                            
p3<-ggplot(qpcr,aes(centered_cooking_time_clustered,log10(corr_bce)))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(size=5,alpha=0.5,aes(color=colony_density),position=position_jitter(width=0.3,seed=4))+
  geom_text(size=3,aes(label=ID),position=position_jitter(width=0.3,seed=4))+
  scale_color_manual(values=paldens)+
  theme_classic()+
  theme(legend.position="none")+
  ylab("log10(genomic equivalents per g food)")+
  xlab("cooking time (clustered center*)")

p4<-ggplot(qpcr,aes(bbd_exp,log10(corr_bce)))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(size=5,alpha=0.5,aes(color=colony_density),position=position_jitter(width=0.3,seed=2))+
  geom_text(size=3,aes(label=ID),position=position_jitter(width=0.3,seed=2))+
  scale_color_manual(values=paldens)+
  theme_classic()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank())+
  labs(color="colony density\non blood agar")+
  xlab("shelf life information")


p5<-ggplot(qpcr,aes(additional_recommendation,log10(corr_bce)))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(size=5,alpha=0.5,aes(color=colony_density),position=position_jitter(width=0.3,seed=4))+
  geom_text(size=3,aes(label=ID),position=position_jitter(width=0.3,seed=4))+
  scale_color_manual(values=paldens)+
  theme_classic()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank())+ 
   labs(color="colony density\non blood agar")+
  xlab("additional information")


p6<-ggplot(qpcr,aes(after_opening_clustered,log10(corr_bce)))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(size=5,alpha=0.5,aes(color=colony_density),position=position_jitter(width=0.3,seed=2))+
  geom_text(size=3,aes(label=ID),position=position_jitter(width=0.3,seed=2))+
  scale_color_manual(values=paldens)+
  theme_classic()+
  theme(legend.position="none")+
  ylab("log10(genomic equivalents per g food)")+
  labs(color="colony density\non blood agar")+
  xlab("remaining shelf life after opening (clustered)")





cowplot::plot_grid(p3,p4,p6,p5,nrow=2, labels = c("A", "B", "C", "D"), label_size = 15, label_x = 0.07, label_y = 0.95, rel_widths = c(1, 1, 1, 1))

ggsave(paste0(outputdir,"S3_supplementplotge_g.tiff"), width=22.5,height=22.5,units="cm" )

