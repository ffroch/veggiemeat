#!/bin/bash
conda activate qiime2-2021.4
workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/veggiemeat
mkdir ${workdir}/08-CombineIsolatesAndMiSeq
 
qiime tools export --input-path ${workdir}/03-qiime/07-filtered_dada2/filtered_table-full-length-uniform-classifier.qza --output-path ${workdir}/08-CombineIsolatesAndMiSeq/table
biom convert -i ${workdir}/08-CombineIsolatesAndMiSeq/table/feature-table.biom -o ${workdir}/08-CombineIsolatesAndMiSeq/table/feature-table.txt --to-tsv
tail -n+2 ${workdir}/08-CombineIsolatesAndMiSeq/table/feature-table.txt | wc -l
qiime tools export --input-path ${workdir}/03-qiime/07-filtered_dada2/final-rep-seqs-full-length-uniform-classifier.qza --output-path ${workdir}/08-CombineIsolatesAndMiSeq/repseqs
qiime tools export --input-path ${workdir}/03-qiime/06-taxonomy/taxonomy-full-length-uniform-classifier.qza --output-path ${workdir}/08-CombineIsolatesAndMiSeq/taxonomy

cd ${workdir}/08-CombineIsolatesAndMiSeq/repseqs/
makeblastdb -in dna-sequences.fasta -title rep-seqs -dbtype nucl -hash_index

