#!/bin/bash
#run in conda env minionpipe
CONDA_BASE=$(conda info --base)
eval "$(conda shell.bash hook)"
projectdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022
workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/raw/Run20221223
#build project folder and move guppy basecalled fastq files
for i in 06 18 07 19 30 42 54 66 78 90; do
workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/raw/Run20221223
   mkdir -p ${workdir}/fastq_unzipped/barcode${i}
    cp ${workdir}/fastq/barcode${i}/*.gz ${workdir}/fastq_unzipped/barcode${i}
    gunzip ${workdir}/fastq_unzipped/barcode${i}/*.gz
done

# Combine fastq
#concenate fastqs and count reads
    for i in 06 18 07 19 30 42 54 66 78 90; do
    workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/raw/Run20221223
        cat ${workdir}/fastq_unzipped/barcode${i}/*.fastq > ${workdir}/fastq_unzipped/barcode${i}/all_records_bc${i}.fastq
        echo barcode${i}
        cat ${workdir}/fastq_unzipped/barcode${i}/all_records_bc${i}.fastq | wc -l | awk '{print $1/4}' >> ${workdir}/concenated_fastqs.txt
    done

# gzip fastq files

for i in 06 18 07 19 30 42 54 66 78 90; do
workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/raw/Run20221223
        gzip ${workdir}/fastq_unzipped/barcode${i}/all_records_bc${i}.fastq
done

## filtlong
projectdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022
mkdir ${projectdir}/data/processed/Run20221223/filtlong
for i in 06 18 07 19 30 42 54 66 78 90; do
workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/raw/Run20221223
  filtlong --min_length 1000 ${workdir}/fastq_unzipped/barcode${i}/all_records_bc${i}.fastq.gz | gzip > ${projectdir}/data/processed/Run20221223/filtlong/barcode${i}.trim.fastq.gz
done


#flye
projectdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022
  mkdir ${projectdir}/data/processed/Run20221223/flye
  for i in 06 18 07 19 30 42 54 66 78 90; do
      flye --nano-raw ${projectdir}/data/processed/Run20221223/filtlong/barcode${i}.trim.fastq.gz -o ${projectdir}/data/processed/Run20221223/flye/barcode${i} -t 32 --iterations 3
  done

  for i in 06 18 07 19 30 42 54 66 78 90; do
    workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/processed/Run20221223
    contigs="${workdir}/flye/barcode${i}/assembly.fasta"
    threads=32
    reads="${workdir}/filtlong/barcode${i}.trim.fastq.gz"

    minimap2 -x ava-ont -t $threads $contigs $reads > ${workdir}/flye/barcode${i}/overlaps_1.paf
    racon -t $threads $reads ${workdir}/flye/barcode${i}/overlaps_1.paf $contigs > ${workdir}/flye/barcode${i}/racon_1.fasta

    minimap2 -x ava-ont -t $threads ${workdir}/flye/barcode${i}/racon_1.fasta $reads > ${workdir}/flye/barcode${i}/overlaps_2.paf
    racon -m 8 -x -6 -g -8 -w 500 -t $threads $reads ${workdir}/flye/barcode${i}/overlaps_2.paf ${workdir}/flye/barcode${i}/racon_1.fasta > ${workdir}/flye/barcode${i}/racon_2.fasta

    minimap2 -x ava-ont -t $threads ${workdir}/flye/barcode${i}/racon_2.fasta $reads > ${workdir}/flye/barcode${i}/overlaps_3.paf
    racon -m 8 -x -6 -g -8 -w 500 -t $threads $reads ${workdir}/flye/barcode${i}/overlaps_3.paf ${workdir}/flye/barcode${i}/racon_2.fasta > ${workdir}/flye/barcode${i}/racon_3.fasta

    minimap2 -x ava-ont -t $threads ${workdir}/flye/barcode${i}/racon_3.fasta $reads > ${workdir}/flye/barcode${i}/overlaps_4.paf
    racon -m 8 -x -6 -g -8 -w 500 -t $threads $reads ${workdir}/flye/barcode${i}/overlaps_4.paf ${workdir}/flye/barcode${i}/racon_3.fasta > ${workdir}/flye/barcode${i}/racon_4.fasta

    medaka_consensus -i $reads -d ${workdir}/flye/barcode${i}/racon_4.fasta -o ${workdir}/flye/barcode${i}/medaka -t $threads -m r941_min_hac_g507
    rm -rf ${workdir}/flye/barcode${i}/*.paf ${workdir}/flye/barcode${i}/racon_*.fasta*
  done

#TORMES
for i in 06 18 07 19 30 42 54 66 78 90; do
workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/processed/Run20221223
  echo $i >> ${projectdir}/one
  echo "GENOME" >> ${projectdir}/two
  realpath ${workdir}/flye/barcode${i}/medaka/consensus.fasta >> ${projectdir}/three
done

 mkdir -p ${workdir}/tormes/
paste ${projectdir}/one ${projectdir}/two ${projectdir}/three | sed "1iSamples\tRead1\tRead2" > ${workdir}/tormes/minion-metadata.txt

conda deactivate

conda activate tormes-1.3.0
tormes -m ${workdir}/tormes/minion-metadata.txt -o ${workdir}/tormes/minion-tormes-20230103 -t 32
conda deactivate

#get length per contig from assembled fasta files
for i in 6 18 7 19 30 42 54 66 78 90; do
  workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/processed/Run20221223
cat ${workdir}/tormes/*/genomes/${i}.fasta | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }'
done

mkdir ${workdir}/tormes/separated_contigs
for i in 6; do
cat ${workdir}/tormes/*/genomes/${i}.fasta | while read line ; do if [ "${line:0:1}" == ">" ]; then echo -e "\n"$line ; else  echo $line | tr -d '\n' ; fi ; done | tail -n+2 > ${workdir}/tormes/separated_contigs/${i}.fasta
done
