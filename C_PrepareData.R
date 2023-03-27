#import and data
data<-read.table("./03-qiime/07-filtered_dada2/exported-table-full-length-uniform-classifier/full-length-uniform-classifier-feature-table.txt",
                 sep="\t",header=T,comment.char = "",skip=1)
meta<-read.table("./01-metadata/samplemeta.txt",sep="\t",header=T)
meta$group<-paste(meta$protein.source, meta$texture, sep="-")
saveRDS(meta, "./06-outputfiles/C_meta.rds")
#rename columns to sample.id
headerold<-data.frame(miseq.id=colnames(data[2:ncol(data)]))
headernew<-meta[,1:2]
headersorted<-join(headerold, headernew)
colnames(data)[2:ncol(data)]<-headersorted$sample.id
saveRDS(data, "./06-outputfiles/C_data.rds")

#loadqPCR data
ge<-read.table("./02-rawdata/genomicequivalentsperg.txt",header=T,sep="\t")
saveRDS(ge, "./06-outputfiles/C_ge.rds")

#data transposition
tdat<-data.frame(t(data[,2:ncol(data)]))
colnames(tdat)<-data[2:length(rownames(data)),1]
rownames(tdat)<-colnames(data[2:ncol(data)])
saveRDS(tdat, "./06-outputfiles/C_tdat.rds")

phylo_object<-qza_to_phyloseq("./03-qiime/07-filtered_dada2/filtered_table-full-length-uniform-classifier.qza",
                              "./03-qiime/09-phylogenetic-tree-full-length-uniform-classifier/rooted-tree.qza",
                              "./03-qiime/06-taxonomy/taxonomy-full-length-uniform-classifier.qza",
                              "./03-qiime/samplemeta_qiime.txt")

sample_data(phylo_object)$group<-paste(sample_data(phylo_object)$protein.source, sample_data(phylo_object)$texture, sep="-")

sample_names(phylo_object)<-sample_data(phylo_object)$sample.id
saveRDS(phylo_object, "./06-outputfiles/C_phylo_object.rds")

phylostats<-data.frame("feature"=NA, "value"=NA)
phylostats[1,]<-c("no. of samples", length(sample_names(phylo_object)))
phylostats[2,]<-c("overall reads", sum(otu_table(phylo_object)))
phylostats[3,]<-c("min reads", min(sample_sums(phylo_object)))
phylostats[4,]<-c("max reads", max(sample_sums(phylo_object)))
phylostats[5,]<-c("median reads", median(sample_sums(phylo_object)))
phylostats[6,]<-c("no. of asvs", ntaxa(phylo_object))
phylostats[7,]<-c("no. of Genera", length(unique(data.frame(tax_table(phylo_object))$Genus)))
phylostats[8,]<-c("no. of Phyla", length(unique(data.frame(tax_table(phylo_object))$Phylum)))
phylostats[,2]<-as.numeric(phylostats[,2])

saveRDS(phylostats, file ="./06-outputfiles/C_phylostats.rds")



#replace d__ in taxa object with empty space
tax <- data.frame(phyloseq::tax_table(phylo_object)[, 1]) %>%
  mutate(Kingdom = stringr::str_replace_all(Kingdom, "d__", ""))
tax_mat <- as.matrix(tax)
rownames(tax_mat) <- phyloseq::taxa_names(phylo_object)
phyloseq::tax_table(phylo_object)[,1] <- tax_mat


isometa<-read.table("./01-metadata/isolatemeta.txt",sep="\t",header=T)
saveRDS(isometa, "./06-outputfiles/C_isometa.rds")