conda create -n BTyper

conda install -c bioconda python=2.7
conda install -c bioconda biopython=1.73
conda install -c bioconda blast
conda install -c bioconda spades
conda install -c bioconda sra-tools

cd /home/roch/miniconda3/envs/BTyper/bin
ls
wget https://github.com/lmc297/BTyper/raw/master/archive/btyper-2.3.3.tar.gz
tar -xzvf btyper-2.3.3.tar.gz
rm btyper-2.3.3.tar.gz

cd btyper-2.3.3


workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data
mkdir ${workdir}/06_BTyper

mkdir ${workdir}/06_BTyper/20221220_90
python btyper -t seq -i ${workdir}/05_rMLST/fasta/20221220_90.fasta -o ${workdir}/06_BTyper/20221220_90 --draft_genome

mkdir ${workdir}/06_BTyper/20220811_42
python btyper -t seq -i ${workdir}/05_rMLST/fasta/20220811_42.fasta -o ${workdir}/06_BTyper/20220811_42 --draft_genome --anib True

mkdir ${workdir}/06_BTyper/20220818_7
python btyper -t seq -i ${workdir}/05_rMLST/fasta/20220818_7.fasta -o ${workdir}/06_BTyper/20220818_7 --draft_genome

mkdir ${workdir}/06_BTyper/20220818_31
python btyper -t seq -i ${workdir}/05_rMLST/fasta/20220818_31.fasta -o ${workdir}/06_BTyper/20220818_31 --draft_genome

mkdir ${workdir}/06_BTyper/20221220_7
python btyper -t seq -i ${workdir}/05_rMLST/fasta/20221220_7.fasta -o ${workdir}/06_BTyper/20221220_7 --draft_genome



python build_btyper_anib_db.py -db published
mkdir ${workdir}/06_BTyper/20220811_42_anib
python btyper -t seq -i ${workdir}/05_rMLST/fasta/20220811_42.fasta -o ${workdir}/06_BTyper/20220811_42_anib --draft_genome --anib True
mkdir ${workdir}/06_BTyper/20220811_42_anib_panClat
python btyper -t seq -i ${workdir}/05_rMLST/fasta/20220811_42.fasta -o ${workdir}/06_BTyper/20220811_42_anib_panClat --draft_genome --anib True --panC_database latest
mkdir ${workdir}/06_BTyper/20221220_90_anib_panClat
python btyper -t seq -i ${workdir}/05_rMLST/fasta/20221220_90.fasta -o ${workdir}/06_BTyper/20221220_90_anib_panClat --draft_genome --anib True --panC_database latest
