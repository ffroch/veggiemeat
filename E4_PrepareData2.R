#remove 439 size sample from phylo object for later all analysis
phylo_object<-readRDS("./06-outputfiles/C_phylo_object.rds")
phylo_object_red<-subset_samples(phylo_object,sample.id!="C4")
saveRDS(phylo_object_red,"./06-outputfiles/E4_phylo_object_red.rds")

rarecred<-readRDS("./06-outputfiles/E3_rarecred.rds")


rarecredstats<-data.frame("feature"=NA, "value"=NA)
rarecredstats[1,]<-c("no. of samples", length(sample_names(rarecred)))
rarecredstats[2,]<-c("overall reads", sum(otu_table(rarecred)))
rarecredstats[3,]<-c("min reads", min(sample_sums(rarecred)))
rarecredstats[4,]<-c("max reads", max(sample_sums(rarecred)))
rarecredstats[5,]<-c("median reads", median(sample_sums(rarecred)))
rarecredstats[6,]<-c("no. of asvs", ntaxa(rarecred))
rarecredstats[7,]<-c("no. of Genera", length(unique(data.frame(tax_table(rarecred))$Genus)))
rarecredstats[8,]<-c("no. of Phyla", length(unique(data.frame(tax_table(rarecred))$Phylum)))
rarecredstats[,2]<-as.numeric(rarecredstats[,2])

saveRDS(rarecredstats, file ="./06-outputfiles/E4_rarecredstats.rds")

tdat<-readRDS("./06-outputfiles/C_tdat.rds")
tdatred<- tdat[row.names(tdat)!="C4",]
saveRDS(tdatred, "./06-outputfiles/E4_tdatred.rds")


#generate a relative abundance asv table
#export for LEfse tool with bioconda see ./07-lefse/01_lefse_pipe.sh
relab_po<-transform_sample_counts(rarecred, function(x) x / sum(x))

saveRDS(relab_po,"./06-outputfiles/E4_relab_po.rds")