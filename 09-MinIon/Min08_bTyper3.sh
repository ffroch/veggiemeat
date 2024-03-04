conda create -n btyper3
conda activate btyper3
conda install -c bioconda btyper
btyper -h

workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data
mkdir ${workdir}/06_BTyper3

mkdir ${workdir}/06_BTyper3/20221220_90
btyper3 -i ${workdir}/05_rMLST/fasta/20221220_90.fasta -o ${workdir}/06_BTyper3/20221220_90 --ani_geneflow True --ani_typestrains True


mkdir ${workdir}/06_BTyper3/20220811_42
btyper3 -i ${workdir}/05_rMLST/fasta/20220811_42.fasta -o ${workdir}/06_BTyper3/20220811_42 --ani_geneflow True --ani_typestrains True

mkdir ${workdir}/06_BTyper3/20220818_7
btyper3 -i ${workdir}/05_rMLST/fasta/20220818_7.fasta -o ${workdir}/06_BTyper3/20220818_7 --ani_geneflow True --ani_typestrains True

mkdir ${workdir}/06_BTyper3/20220818_31
btyper3 -i ${workdir}/05_rMLST/fasta/20220818_31.fasta -o ${workdir}/06_BTyper3/20220818_31 --ani_geneflow True --ani_typestrains True

mkdir ${workdir}/06_BTyper3/20221220_7
btyper3 -i ${workdir}/05_rMLST/fasta/20221220_7.fasta -o ${workdir}/06_BTyper3/20221220_7 --ani_geneflow True --ani_typestrains True

