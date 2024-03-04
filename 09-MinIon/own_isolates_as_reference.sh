workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data
mkdir ${workdir}/10_ownref

conda activate blast
for file in /data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/05_rMLST/fasta/*.fasta; do
fl=$(basename "$file")
makeblastdb -in ${file} -out ${workdir}/10_ownref/${fl} -dbtype 'nucl' -hash_index
done

blastn -query /data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/00_ReferenceFasta/cer_gene.fasta -task blastn -db /data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/10_ownref/20220811_42.fasta -out /data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/10_ownref/test.txt -evalue 1e-5 -word_size 4 -num_threads 2 -outfmt '6 qseqid sseqid pident sstart send qstart qend evalue bitscore score qlen length sseq'

tblastx -query /data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/00_ReferenceFasta/cer_gene.fasta -db /data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/10_ownref/20220811_42.fasta -out /data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/10_ownref/test.txt -evalue 0.01 -word_size 4 -num_threads 2 -outfmt '7 qseqid sseqid pident sstart send qstart qend evalue bitscore score qlen length sseq'
