



#load DNA string set of trimmed isolate sequences
stringset<-readDNAStringSet("./06-outputfiles/seqstrimclean_40_15_100bp.fa","fasta")
length(stringset)

#load reference data
ref_fasta<-"./01-metadata/silva_nr99_v138.1_train_set.fa.gz"


#assign Taxonomy
taxtab<-assignTaxonomy(stringset,refFasta=ref_fasta)
#write taxonomy table
write.table(taxtab,"./06-outputfiles/taxonomy_40_15_100.txt",sep=";")

