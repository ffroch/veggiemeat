conda create -n checkm python=3.9
conda activate checkm
conda install numpy matplotlib pysam
conda install hmmer prodigal pplacer
pip3 install checkm-genome # NOT working
conda install -c bioconda checkm-genome

cd /data/Unit_LMM/selberherr-group/roch/00_databases
mkdir checkm_db
cd checkm_db

wget  wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar -xzvf checkm_data_2015_01_16.tar.gz

export CHECKM_DATA_PATH=/data/Unit_LMM/selberherr-group/roch/00_databases/checkm_db

#get data for checkm
workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/
mkdir ${workdir}08_checkm
#copy genomes directories and rename it
for file in ${workdir}03_MinIon_Tormes/minion-tormes*/annotation/*/*.fna;do
   cp ${file} $(dirname "$file" )/..
 done
 for file in ${workdir}03_MinIon_Tormes/minion-tormes*/annotation/*.fna;do
    mv ${file} $(dirname "$file" )/..
  done
for file in ${workdir}03_MinIon_Tormes/minion-tormes*/*.fna; do
  fp=$(dirname "$file"); fl=$(basename "$file"); mv "$fp/$fl" "$fp"_"$fl"
done

mkdir ${workdir}08_checkm/fna_folder
#move genomes directories to new folder
for file in ${workdir}03_MinIon_Tormes/*.fna; do
 mv $file ${workdir}08_checkm/fna_folder
done
#shorten the filenames
for filename in ${workdir}08_checkm/fna_folder/*; do mv "$filename" "$(echo "$filename" | sed -e 's/minion-tormes-//g')";  done

mkdir ${workdir}08_checkm/results
workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/
checkm lineage_wf ${workdir}08_checkm/fna_folder ${workdir}08_checkm/results -t20

checkm qa ${workdir}08_checkm/results/lineage.ms /data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/08_checkm/results > ${workdir}08_checkm/results/quality_assessment_new.txt

checkm qa ${workdir}08_checkm/results/lineage.ms ${workdir}08_checkm/results > ${workdir}08_checkm/results/quality_assessment.txt
