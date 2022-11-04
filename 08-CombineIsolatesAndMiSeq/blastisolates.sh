#!/bin/bash
conda activate qiime2-2021.4
workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/veggiemeat

# Outfmt = 6
blastn  -query ${workdir}/06-outputfiles/seqstrimclean_40_15_100bp.fa -db ${workdir}/08-CombineIsolatesAndMiSeq/repseqs/dna-sequences.fasta  -out ${workdir}/08-CombineIsolatesAndMiSeq/comparison.blastout -num_alignments 25 -outfmt "6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore"

# Outfmt = 6 -- best ones
blastn  -query ${workdir}/06-outputfiles/seqstrimclean_40_15_100bp.fa -db ${workdir}/08-CombineIsolatesAndMiSeq/repseqs/dna-sequences.fasta  -out ${workdir}/08-CombineIsolatesAndMiSeq/best.blastout -num_alignments 1 -outfmt "6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore"

conda deactivate