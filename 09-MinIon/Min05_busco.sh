conda create -n busco
conda install -c conda-forge -c bioconda busco=5.4.4

workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data
mkdir ${workdir}/07_busco

workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data
mkdir ${workdir}/07_busco/results

cd ${workdir}/07_busco/results
for file in ${workdir}/08_checkm/fna_folder/*.fna; do
  fl=$(basename "$file")
  busco -m genome -i ${file} -o ${fl} --auto-lineage-prok --cpu 64
done

workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data
cd ${workdir}/07_busco/results
#run specific busco example
busco -m genome -i ${workdir}/08_checkm/fna_folder/20221220_42.fna -o specificLeuco -l lactobacillales_odb10 -f
busco -m genome -i ${workdir}/08_checkm/fna_folder/20221220_42.fna -o specificLeucbac -l bacteria_odb10  -f

workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data
cd ${workdir}/07_busco/results
busco -m genome -i /data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/03_MinIon_Tormes/minion-tormes-20221220/annotation/18_annotation/18.ffn -o specificLeuco -l lactobacillales_odb10 -f
busco -m prot -i /data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/03_MinIon_Tormes/minion-tormes-20221220/annotation/18_annotation/18.faa -o specificLeuco -l lactobacillales_odb10 -f
