conda activate minionpipe
workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/01_MinIon_RawData/Minion_20220128
mkdir ${workdir}/fastq_unzipped
for i in 02 49 61 73 85; do
mkdir ${workdir}/fastq_unzipped/barcode${i}
        cat ${workdir}/barcode${i}/*/*/*/fastq_pass/*.fastq > ${workdir}/fastq_unzipped/barcode${i}/all_records_bc${i}.fastq
        echo barcode${i}
        cat ${workdir}/fastq_unzipped/barcode${i}/all_records_bc${i}.fastq | wc -l | awk '{print $1/4}' >> ${workdir}/concenated_fastqs.txt
    done

    for i in 02 49 61 73 85; do
            gzip ${workdir}/fastq_unzipped/barcode${i}/all_records_bc${i}.fastq
    done

    ## filtlong
workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data
    for i in 02 49 61 73 85; do
      filtlong --min_length 1000 ${workdir}/01_MinIon_RawData/Minion_20220128/fastq_unzipped/barcode${i}/all_records_bc${i}.fastq.gz | gzip > ${workdir}/02_MinIon_ProcessedData/filtlong/20220128_barcode${i}.trim.fastq.gz
    done

workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data    
    #flye
      for i in 02 49 61 73 85; do
          flye --nano-raw ${workdir}/02_MinIon_ProcessedData/filtlong/20220128_barcode${i}.trim.fastq.gz -o ${workdir}/02_MinIon_ProcessedData/flye/20220128_barcode${i} -t 32 --iterations 3
      done


  workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/02_MinIon_ProcessedData
      for i in 02 49 61 73 85; do
        contigs="${workdir}/flye/20220128_barcode${i}/assembly.fasta"
        threads=32
        reads="${workdir}/filtlong/20220128_barcode${i}.trim.fastq.gz"

        minimap2 -x ava-ont -t $threads $contigs $reads > ${workdir}/flye/20220128_barcode${i}/overlaps_1.paf
        racon -t $threads $reads ${workdir}/flye/20220128_barcode${i}/overlaps_1.paf $contigs > ${workdir}/flye/20220128_barcode${i}/racon_1.fasta

        minimap2 -x ava-ont -t $threads ${workdir}/flye/20220128_barcode${i}/racon_1.fasta $reads > ${workdir}/flye/20220128_barcode${i}/overlaps_2.paf
        racon -m 8 -x -6 -g -8 -w 500 -t $threads $reads ${workdir}/flye/20220128_barcode${i}/overlaps_2.paf ${workdir}/flye/20220128_barcode${i}/racon_1.fasta > ${workdir}/flye/20220128_barcode${i}/racon_2.fasta

        minimap2 -x ava-ont -t $threads ${workdir}/flye/20220128_barcode${i}/racon_2.fasta $reads > ${workdir}/flye/20220128_barcode${i}/overlaps_3.paf
        racon -m 8 -x -6 -g -8 -w 500 -t $threads $reads ${workdir}/flye/20220128_barcode${i}/overlaps_3.paf ${workdir}/flye/20220128_barcode${i}/racon_2.fasta > ${workdir}/flye/20220128_barcode${i}/racon_3.fasta

        minimap2 -x ava-ont -t $threads ${workdir}/flye/20220128_barcode${i}/racon_3.fasta $reads > ${workdir}/flye/20220128_barcode${i}/overlaps_4.paf
        racon -m 8 -x -6 -g -8 -w 500 -t $threads $reads ${workdir}/flye/20220128_barcode${i}/overlaps_4.paf ${workdir}/flye/20220128_barcode${i}/racon_3.fasta > ${workdir}/flye/20220128_barcode${i}/racon_4.fasta

        medaka_consensus -i $reads -d ${workdir}/flye/20220128_barcode${i}/racon_4.fasta -o ${workdir}/flye/20220128_barcode${i}/medaka -t $threads -m r941_min_hac_g507
        rm -rf ${workdir}/flye/barcode${i}/*.paf ${workdir}/flye/20220128_barcode${i}/racon_*.fasta*
      done




    workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data
    #TORMES
    for i in 02 49 61 73 85; do
      echo $i >> ${workdir}/one
      echo "GENOME" >> ${workdir}/two
      realpath ${workdir}/02_MinIon_ProcessedData/flye/20220128_barcode${i}/medaka/consensus.fasta >> ${workdir}/three
    done

     paste ${workdir}/one ${workdir}/two ${workdir}/three | sed "1iSamples\tRead1\tRead2" > ${workdir}/03_MinIon_Tormes/minion-metadata.txt

    conda deactivate

    conda activate tormes-1.3.0

    cd ${workdir}/03_MinIon_Tormes
    tormes -m ${workdir}/03_MinIon_Tormes/minion-metadata.txt -o minion-tormes-20220128 -t 32
    conda deactivate
