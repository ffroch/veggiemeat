

#beta-diversity is done with 100 iterations of coverage based rarefaction
if(file.exists("./06-outputfiles/F0_rare100.rds"))
{
  rare100<-readRDS("./06-outputfiles/F0_rare100.rds")
}   else   {
  set.seed(12345)
  rare100<-phyloseq_coverage_raref(phylo_object_red,iter=100,coverage=0.995)
  saveRDS(rare100,"./06-outputfiles/F0_rare100.rds" )
}

# generate several distance matrices
bd<-mult_dissim(rare100, method="bray",average=T)
wuni<-mult_dissim(rare100, method="wunifrac",average=T)
uni<-mult_dissim(rare100, method="unifrac",average=T)
jd<-mult_dissim(rare100, method="jaccard",average=T)
jsd<-mult_dissim(rare100, method="jsd",average=T)



dmlist<-list(bd, wuni, uni, jd, jsd)
dmlistnames<-c("bray", "wunifrac", "unifrac", "jaccard", "jsd")

saveRDS(dmlist,"./06-outputfiles/F0_dmlist.rds")
saveRDS(bd,"./06-outputfiles/F0_bd.rds")
saveRDS(wuni,"./06-outputfiles/F0_wuni.rds")
saveRDS(uni,"./06-outputfiles/F0_uni.rds")
saveRDS(jd,"./06-outputfiles/F0_jd.rds")
saveRDS(jsd,"./06-outputfiles/F0_jsd.rds")