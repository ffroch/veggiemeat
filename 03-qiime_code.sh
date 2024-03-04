#!/bin/bash
conda activate qiime2-2021.4

#3 - dada2
mkdir /data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/veggiemeat/03-qiime/addon_270_210
workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/veggiemeat/03-qiime/addon_270_210
cd ${workdir}

mkdir 05-dada2

qiime dada2 denoise-paired \
	--i-demultiplexed-seqs 04-demux/demux-paired-end.qza \
	--p-trim-left-f 15 \
	--p-trim-left-r 15 \
	--p-trunc-len-f 270 \
	--p-trunc-len-r 210 \
	--p-max-ee-f 2 \
	--p-max-ee-r 2 \
	--o-representative-sequences 05-dada2/rep-seqs.qza \
	--o-table 05-dada2/table.qza \
	--o-denoising-stats 05-dada2/stats-dada2.qza \
	
#4 - classification
mkdir 06-taxonomy
for i in silva-138-99-nb-classifier gg-13-8-99-nb-classifier full-length-uniform-classifier; do
qiime feature-classifier classify-sklearn \
	--i-classifier /data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/veggiemeat/03-qiime/00-pretrained_classifier/${i}.qza \
	--i-reads 05-dada2/rep-seqs.qza \
	--o-classification 06-taxonomy/taxonomy-${i}.qza
	done
#5 visualize taxonomy output
for i in silva-138-99-nb-classifier gg-13-8-99-nb-classifier full-length-uniform-classifier; do
qiime metadata tabulate \
	--m-input-file 06-taxonomy/taxonomy-${i}.qza \
	--o-visualization 06-taxonomy/taxonomy-${i}.qzv
	done
	
#6 exclude mitochondrial and chloroplast sequences from table.qza
mkdir 07-filtered_dada2
for i in silva-138-99-nb-classifier gg-13-8-99-nb-classifier full-length-uniform-classifier; do
qiime taxa filter-table \
	--i-table 05-dada2/table.qza \
	--i-taxonomy 06-taxonomy/taxonomy-${i}.qza \
	--p-exclude mitochondria,chloroplast \
	--o-filtered-table 07-filtered_dada2/filtered_table-${i}.qza
done


#7 exclude mitochondrial and chloroplast sequences from rep-seqs.qza
for i in silva-138-99-nb-classifier gg-13-8-99-nb-classifier full-length-uniform-classifier; do
qiime taxa filter-seqs \
	--i-sequences 05-dada2/rep-seqs.qza \
	--i-taxonomy 06-taxonomy/taxonomy-${i}.qza \
	--p-exclude mitochondria,chloroplast \
	--o-filtered-sequences 07-filtered_dada2/final-rep-seqs-${i}.qza
done

#10 final table visualization
for i in silva-138-99-nb-classifier gg-13-8-99-nb-classifier full-length-uniform-classifier; do
qiime feature-table summarize \
	--i-table 07-filtered_dada2/filtered_table-${i}.qza \
	--o-visualization 07-filtered_dada2/final_table-${i}.qzv \
	--m-sample-metadata-file samplemeta_qiime.txt
done

#10a relative abundance asvs
for i in silva-138-99-nb-classifier gg-13-8-99-nb-classifier full-length-uniform-classifier; do

qiime feature-table relative-frequency \
	--i-table 07-filtered_dada2/filtered_table-${i}.qza \
	--o-relative-frequency-table 07-filtered_dada2/filtered_table-${i}-relfreqs.qza

qiime tools export \
	--input-path 07-filtered_dada2/filtered_table-${i}-relfreqs.qza \
	--output-path 07-filtered_dada2/exported-table-${i}-relfreqs

biom convert \
	-i 07-filtered_dada2/exported-table-${i}-relfreqs/feature-table.biom \
	-o 07-filtered_dada2/exported-table-${i}-relfreqs/${i}-feature-table.txt --to-tsv
done

#11 final rep seq visualization
for i in silva-138-99-nb-classifier gg-13-8-99-nb-classifier full-length-uniform-classifier; do
qiime feature-table tabulate-seqs \
	--i-data 07-filtered_dada2/final-rep-seqs-${i}.qza \
	--o-visualization 07-filtered_dada2/final-rep-seqs-${i}.qzv
done	

#12 export
for i in silva-138-99-nb-classifier gg-13-8-99-nb-classifier full-length-uniform-classifier; do
 qiime tools export \
	--input-path 07-filtered_dada2/final-rep-seqs-${i}.qza \
	--output-path 08-asv-seqs-${i}
done

#13 make phylogenetic trees
for i in silva-138-99-nb-classifier gg-13-8-99-nb-classifier full-length-uniform-classifier; do
mkdir 09-phylogenetic-tree-${i}
qiime alignment mafft \
	--i-sequences 07-filtered_dada2/final-rep-seqs-${i}.qza \
	--o-alignment 09-phylogenetic-tree-${i}/aligned-rep-seqs.qza
qiime alignment mask \
	--i-alignment 09-phylogenetic-tree-${i}/aligned-rep-seqs.qza \
	--o-masked-alignment 09-phylogenetic-tree-${i}/masked-aligned-rep-seqs.qza
qiime phylogeny fasttree \
	--i-alignment 09-phylogenetic-tree-${i}/masked-aligned-rep-seqs.qza \
	--o-tree 09-phylogenetic-tree-${i}/unrooted-tree.qza
qiime phylogeny midpoint-root \
	--i-tree 09-phylogenetic-tree-${i}/unrooted-tree.qza \
	--o-rooted-tree 09-phylogenetic-tree-${i}/rooted-tree.qza
done

#14 generating taxa plots
for i in silva-138-99-nb-classifier gg-13-8-99-nb-classifier full-length-uniform-classifier; do
mkdir 10-taxa-plots
qiime taxa barplot \
	--i-table 07-filtered_dada2/filtered_table-${i}.qza \
	--i-taxonomy 06-taxonomy/taxonomy-${i}.qza \
	--m-metadata-file samplemeta_qiime.txt \
	--o-visualization 10-taxa-plots/taxa-bar-plots-${i}.qzv
done

#15 rarefaction
mkdir 11-diversity-rarefaction
for i in silva-138-99-nb-classifier gg-13-8-99-nb-classifier full-length-uniform-classifier; do
qiime diversity alpha-rarefaction \
	--i-table 07-filtered_dada2/filtered_table-${i}.qza \
	--i-phylogeny 09-phylogenetic-tree-${i}/rooted-tree.qza \
	--p-max-depth 24630 \
	--m-metadata-file samplemeta_qiime.txt \
	--o-visualization 11-diversity-rarefaction/alpha-rarefaction-${i}.qzv
done


#16 diversity
for i in silva-138-99-nb-classifier gg-13-8-99-nb-classifier full-length-uniform-classifier; do
qiime diversity core-metrics-phylogenetic \
	--i-phylogeny 09-phylogenetic-tree-${i}/rooted-tree.qza \
	--i-table 07-filtered_dada2/filtered_table-${i}.qza \
	--p-sampling-depth 8000 \
	--m-metadata-file samplemeta_qiime.txt \
	--output-dir 12-diversity-metrics-results-${i}
done

#17 diversity faith-pd
for i in silva-138-99-nb-classifier gg-13-8-99-nb-classifier full-length-uniform-classifier; do
qiime diversity alpha-group-significance \
	--i-alpha-diversity 12-diversity-metrics-results-${i}/faith_pd_vector.qza \
	--m-metadata-file samplemeta_qiime.txt \
	--o-visualization 12-diversity-metrics-results-${i}/faith-pd-group-significance.qzv
done

#18 diversity shannon 
for i in silva-138-99-nb-classifier gg-13-8-99-nb-classifier full-length-uniform-classifier; do
qiime diversity alpha-group-significance \
	--i-alpha-diversity 12-diversity-metrics-results-${i}/shannon_vector.qza \
	--m-metadata-file samplemeta_qiime.txt \
	--o-visualization 12-diversity-metrics-results-${i}/shannon-group-significance.qzv
done	
	
#19 extract taxa bar plots/taxa-bar-plots
for i in silva-138-99-nb-classifier gg-13-8-99-nb-classifier full-length-uniform-classifier; do	
cp 10-taxa-plots/taxa-bar-plots-${i}.qzv 10-taxa-plots/taxa-bar-plots-${i}.zip
unzip 10-taxa-plots/taxa-bar-plots-${i}.zip -d 10-taxa-plots/taxa-bar-plots-${i}
rm 10-taxa-plots/taxa-bar-plots-${i}.zip
done


#20 relative abundances
for i in silva-138-99-nb-classifier gg-13-8-99-nb-classifier full-length-uniform-classifier; do
mkdir 13-relative-abundances-${i}
	for j in 1 2 3 4 5 6; do

qiime taxa collapse \
	--i-table 07-filtered_dada2/filtered_table-${i}.qza \
	--i-taxonomy 06-taxonomy/taxonomy-${i}.qza --p-level ${j} \
	--o-collapsed-table 13-relative-abundances-${i}/table-level-${j}.qza

qiime feature-table relative-frequency \
	--i-table 13-relative-abundances-${i}/table-level-${j}.qza \
	--o-relative-frequency-table 13-relative-abundances-${i}/table-level-${j}-relfreqs.qza

qiime tools export \
	--input-path 13-relative-abundances-${i}/table-level-${j}-relfreqs.qza \
	--output-path 13-relative-abundances-${i}/exported-table-level-${j}-relfreqs

biom convert \
	-i 13-relative-abundances-${i}/exported-table-level-${j}-relfreqs/feature-table.biom \
	-o 13-relative-abundances-${i}/exported-table-level-${j}-relfreqs/${j}-feature-table.txt --to-tsv
	done
done






conda deactivate