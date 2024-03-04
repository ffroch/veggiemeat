#protocol 03.01.Run20221223

#conda create -n antismash antismash # NOT working

conda create -n antismash
conda activate antismash
#conda install -c bioconda antismash # NOT working
#conda install -c "bioconda/label/cf201901" antismash # NOT working


#conda install -c bioconda hmmer2 # NOT working
#conda install -c bioconda hmmer # NOT working
#conda install -c bioconda diamond # NOT working
#conda install -c bioconda fasttree # NOT working
#conda install -c bioconda prodigal # NOT working
#conda install -c bioconda blast # NOT working
#conda install -c bioconda muscle # NOT working
#conda install -c bioconda glimmerhmm # NOT working

conda install -y diamond=2.0.9
conda install -y fasttree=2.1.11
conda install -y GlimmerHMM=3.0.4
conda install -y hmmer2=2.3.2
conda install -y hmmer3=3.1b2
conda activate antismash
conda install -y diamond=2.0.9
conda install -y fasttree=2.1.11
conda install -y GlimmerHMM=3.0.4
conda install -y hmmer2=2.3.2
conda install -y hmmer3
conda install -y hmmer=3.3.2
conda install -y meme=4.11.2
conda install -y muscle=3.8.1551
conda install -y blast=2.10.0
conda install -y prodigal=2.6.3
conda install -y blast=2.6.0
conda install -y blast=2.10.1
conda install -y python=3.8
conda install -y java-jdk


wget https://dl.secondarymetabolites.org/releases/6.1.1/antismash-6.1.1.tar.gz
tar -zxf antismash-6.1.1.tar.gz
pip install ./antisamsh-6.1.1
download-antismash-databases
antismash --check-prereqs


mkdir /data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/processed/Run20221223/antismash
for i in 6 7 18 19 30 42 54 66 78 90; do
  workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/processed/Run20221223
  mkdir ${workdir}/antismash/barcode_${i}
  cd ${workdir}/antismash/barcode_${i}
  antismash ${workdir}/tormes/minion-tormes-20230103/annotation/${i}_annotation/${i}.gbk
done

workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/03_MinIon_Tormes
for i in 6; do
antismash --cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go --smcog-trees ${workdir}/minion-tormes-20220811/annotation/${i}_annotation/${i}.gbk
done

cd /data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/04_MinIon_AntiSmash_detailed
workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/03_MinIon_Tormes
for i in 6 18 42 30; do
  antismash --cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go --smcog-trees ${workdir}/minion-tormes-20221220/annotation/${i}_annotation/${i}.gbk
done

for i in 15 27; do
  antismash --cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go --smcog-trees ${workdir}/minion-tormes-20220522_24/annotation/${i}_annotation/${i}.gbk
done

antismash --cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go --smcog-trees ${workdir}/minion-tormes-20220505/annotation/40_annotation/40.gbk
