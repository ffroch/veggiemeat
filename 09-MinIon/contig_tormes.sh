#!/bin/bash
conda activate tormes-1.3.0
#get length per contig from assembled fasta files
for i in 6 18 7 19 30 42 54 66 78 90; do
  workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/processed/Run20221223
cat ${workdir}/tormes/*/genomes/${i}.fasta | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }'
done

mkdir ${workdir}/tormes/separated_contigs
for i in 6 18 7 19 30 42 54 66 78 90; do
mkdir ${workdir}/tormes/separated_contigs/barcode_${i}
cat ${workdir}/tormes/*/genomes/${i}.fasta | while read line ; do if [ "${line:0:1}" == ">" ]; then echo -e "\n"$line ; else  echo $line | tr -d '\n' ; fi ; done | tail -n+2 > ${workdir}/tormes/separated_contigs/barcode_${i}/allcontigs.fasta
done

for i in 6 18 7 19 30 42 54 66 78 90; do
  cp ${workdir}/tormes/splitfasta.sh ${workdir}/tormes/separated_contigs/barcode_${i}
  cd ${workdir}/tormes/separated_contigs/barcode_${i}
  ./splitfasta.sh allcontigs.fasta
  rm ./splitfasta.sh
  rm ./allcontigs.fasta
done

for i in 6 18 7 19 30 42 54 66 78 90; do
workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/processed/Run20221223/tormes/separated_contigs/barcode_${i}/
cd ${workdir}
for f in *.fasta; do
  echo $f >> ${workdir}/one
  echo "GENOME" >> ${workdir}/two
  realpath ${workdir}/${f} >> ${workdir}/three
done
done

for i in 6 18 7 19 30 42 54 66 78 90; do
workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/processed/Run20221223/tormes/separated_contigs/barcode_${i}/
mkdir -p ${workdir}/contig_tormes/
paste ${workdir}/one ${workdir}/two ${workdir}/three | sed "1iSamples\tRead1\tRead2" > ${workdir}/contig_tormes/minion-metadata.txt
done


for i in 6 18 7 19 30 42 54 66 78 90; do
  workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/processed/Run20221223/tormes/separated_contigs/barcode_${i}/contig_tormes
  cd ${workdir}
tormes -m minion-metadata.txt -o minion-tormes -t 32
done

conda deactivate

for i in 6 18 7 19 30 42 54 66 78 90; do
  workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/processed/Run20221223/tormes/separated_contigs/barcode_${i}/
  cd ${workdir}
  rm one
  rm two
  rm three
done


for i in 6 18 7 19 30 42 54 66 78 90; do
  workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/processed/Run20221223/tormes/separated_contigs/barcode_${i}/contig_tormes/
  cd ${workdir}
  rm -rf minion-tormes
done
