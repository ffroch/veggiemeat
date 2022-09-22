#import and data
data<-read.table("./03-qiime/07-filtered_dada2/exported-table-full-length-uniform-classifier/full-length-uniform-classifier-feature-table.txt",
                 sep="\t",header=T,comment.char = "",skip=1)
meta<-read.table("./01-metadata/samplemeta.txt",sep="\t",header=T)
meta$group<-paste(meta$protein.source, meta$texture, sep="-")
#rename columns to sample.id
headerold<-data.frame(miseq.id=colnames(data[2:ncol(data)]))
headernew<-meta[,1:2]
headersorted<-join(headerold, headernew)
colnames(data)[2:ncol(data)]<-headersorted$sample.id

#loadqPCR data
ge<-read.table("./02-rawdata/genomicequivalentsperg.txt",header=T,sep="\t")


#data transposition
tdat<-data.frame(t(data[,2:ncol(data)]))
colnames(tdat)<-data[2:length(rownames(data)),1]
rownames(tdat)<-colnames(data[2:ncol(data)])

phylo_object<-qza_to_phyloseq("./03-qiime/07-filtered_dada2/filtered_table-full-length-uniform-classifier.qza",
                              "./03-qiime/09-phylogenetic-tree-full-length-uniform-classifier/rooted-tree.qza",
                              "./03-qiime/06-taxonomy/taxonomy-full-length-uniform-classifier.qza",
                              "./03-qiime/samplemeta_qiime.txt")

sample_data(phylo_object)$group<-paste(sample_data(phylo_object)$protein.source, sample_data(phylo_object)$texture, sep="-")

sample_names(phylo_object)<-sample_data(phylo_object)$sample.id

#replace d__ in taxa object with empty space
tax <- data.frame(phyloseq::tax_table(phylo_object)[, 1]) %>%
  mutate(Kingdom = stringr::str_replace_all(Kingdom, "d__", ""))
tax_mat <- as.matrix(tax)
rownames(tax_mat) <- phyloseq::taxa_names(phylo_object)
phyloseq::tax_table(phylo_object)[,1] <- tax_mat
