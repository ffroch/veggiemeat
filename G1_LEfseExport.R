#generate a relative abundance asv table
#export for LEfse tool with bioconda see ./07-lefse/01_lefse_pipe.sh
relab_po<-readRDS("./06-outputfiles/E4_relab_po.rds")

phyloseq2lefse(relab_po, covars="group", file.name="./07-lefse/lefse_data.txt",
               taxa.levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
               transpose.otus = TRUE) #https://github.com/kstagaman/phyloseqCompanion/blob/master/R/phyloseq2lefse.R






