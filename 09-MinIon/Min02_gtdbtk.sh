conda create -n gtdbtk-2.1.1 -c conda-forge -c bioconda gtdbtk=2.1.1

download-db.sh
conda env config vars set GTDBTK_DATA_PATH=/home/roch/miniconda3/envs/gtdbtk-2.1.1/share/gtdbtk-2.1.1/db

workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/
mkdir ${workdir}09_gtdbtk
gtdbtk classify_wf --genome_dir ${workdir}08_checkm/fna_folder --out_dir ${workdir}09_gtdbtk --extension fna --cpus 40

gtdbtk classify_wf --genome_dir ${workdir}08_checkm/fna_folder_test --out_dir ${workdir}09_gtdbtk --extension fna --cpus 40

gtdbtk classify_wf --genome_dir ${workdir}08_checkm/fna_folder --out_dir ${workdir}09_gtdbtk --extension fna --cpus 100

workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data


gtdbtk classify_wf --genome_dir ${workdir}/08_checkm/fna_folder --out_dir ${workdir}/09_gtdbtk --extension fna --cpus 64


for file in ${workdir}/
