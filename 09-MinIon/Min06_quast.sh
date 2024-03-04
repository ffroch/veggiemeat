conda create -n quast
conda acticate quast
conda install -c bioconda quast

workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data
mkdir ${workdir}/11_quast

python ~/miniconda3/envs/quast/bin/quast.py -o ${workdir}/11_quast -r ${workdir}/00_ReferenceFasta/NC_007795_1_Staphylococcus_aureus_subsp_aureus_NCTC8325_chromosome.fasta ${workdir}/03_MinIon_Tormes/minion-tormes-20220522_24/genomes/86.fasta --conserved-genes-finding
python ~/miniconda3/envs/quast/bin/quast.py -o ${workdir}/11_quast ${workdir}/03_MinIon_Tormes/minion-tormes-20220522_24/genomes/86.fasta --conserved-genes-finding

python ~/miniconda3/envs/quast/bin/quast.py -o ${workdir}/11_quast -r ${workdir}/00_ReferenceFasta/NC_006270.3_Bacillus_licheniformis.fasta ${workdir}/03_MinIon_Tormes/minion-tormes-20221220/genomes/66.fasta --conserved-genes-finding
