#get data for roary
workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data
mkdir ${workdir}/13_roary
#copy genomes directories and rename it
for file in ${workdir}/03_MinIon_Tormes/minion-tormes*/annotation/*/*.gff;do
   cp ${file} $(dirname "$file" )/..
 done
 for file in ${workdir}/03_MinIon_Tormes/minion-tormes*/annotation/*.gff;do
    mv ${file} $(dirname "$file" )/..
  done
for file in ${workdir}/03_MinIon_Tormes/minion-tormes*/*.gff; do
  fp=$(dirname "$file"); fl=$(basename "$file"); mv "$fp/$fl" "$fp"_"$fl"
done

for file in ${workdir}/03_MinIon_Tormes/*.gff; do
 mv $file ${workdir}/13_roary
done

#shorten the filenames
for filename in ${workdir}/13_roary/*; do mv "$filename" "$(echo "$filename" | sed -e 's/minion-tormes-//g')";  done

workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data

tar -xvf ${workdir}/13_roary/Bacillus/genome_assemblies_genome_gff1.tar -C ${workdir}/13_roary/Bacillus/

for i in 2 3 4 5 6; do
tar -xvf ${workdir}/13_roary/Bacillus/genome_assemblies_genome_gff${i}.tar -C ${workdir}/13_roary/Bacillus/
done

for file in ${workdir}/13_roary/Bacillus/ncbi*/*.gz; do
  mv ${file} $(dirname "$file" )/..
done

for file in ${workdir}/13_roary/Bacillus/*.gz; do
gunzip ${file}
done

for file in ${workdir}/13_roary/Bacillus/*.tar; do
rm ${file}
done

rm -r ${workdir}/13_roary/Bacillus/ncbi*
rm -r ${workdir}/13_roary/Bacillus/report.txt

workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data
roary -e --mafft -p 20 ${workdir}/13_roary/test/*.gff -f ${workdir}/13_roary/test

roary  ${workdir}/13_roary/Bacillus/*.gff -p 64 -f ${workdir}/13_roary/Bacillus/ -e --mafft -r
