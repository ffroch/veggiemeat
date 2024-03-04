workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/


#copy genomes directories and rename it
for dir in ${workdir}03_MinIon_Tormes/minion-tormes*/genomes; do
  dp=$(dirname "$dir")
 cp -R $dir ${dp}_genomes
done

#move genomes directories to new folder
for dir in ${workdir}03_MinIon_Tormes/minion-tormes*_genomes/; do
 mv $dir ${workdir}05_rMLST
done

#rename the files in the directories
for files in ${workdir}05_rMLST/minion-tormes*_genomes/*.fasta; do
 fp=$(dirname "$files"); fl=$(basename "$files"); mv "$fp/$fl" "$fp"_"$fl"
done

#move them to new folder
mkdir ${workdir}05_rMLST/fasta
for files in ${workdir}05_rMLST/*.fasta; do
  mv ${files} ${workdir}05_rMLST/fasta/$(basename "$files")
done

#delete old unused folders
for dir in ${workdir}05_rMLST/*_genomes; do
  rm -d ${dir}
done

#shorten the filenames
for filename in ${workdir}05_rMLST/fasta/*; do mv "$filename" "$(echo "$filename" | sed -e 's/minion-tormes-//g')";  done
for filename in ${workdir}05_rMLST/fasta/*; do mv "$filename" "$(echo "$filename" | sed -e 's/_genomes//g')";  done



mkdir ${workdir}05_rMLST/output
for files in ${workdir}05_rMLST/fasta/*; do
(echo -n '{"base64":true,"details":true,"sequence": "'; base64 ${files}; echo '"}') | curl -s -H "Content-Type: application/json" -X POST "http://rest.pubmlst.org/db/pubmlst_rmlst_seqdef_kiosk/schemes/1/sequence" -d @- >> ${workdir}05_rMLST/output/$(basename "$files").json
done

for filename in ${workdir}05_rMLST/output/*; do mv "$filename" "$(echo "$filename" | sed -e 's/.fasta//g')";  done

for files in ${workdir}05_rMLST/output/*.json; do
if ! echo "${files}" | grep -q "taxon_prediction"; then
  part1=$(echo -n $(basename "${files}"))
  part2=$(grep -o "taxon_prediction.*" < ${workdir}05_rMLST/output/$(basename "${files}"))
  part3=$(echo $(expr substr "${part2}" 1 300))
  part5=$(echo "${part3}" | grep -o -P '(?<=taxon_prediction).*(?=\})')
  echo ${part1} >> ${workdir}05_rMLST/one
  echo ${part5} >> ${workdir}05_rMLST/two
fi
done

paste ${workdir}05_rMLST/one ${workdir}05_rMLST/two > ${workdir}05_rMLST/rMLST_ouput_summary.txt
rm ${workdir}05_rMLST/one
rm ${workdir}05_rMLST/two


#panC search
for files in ${workdir}03_MinIon_Tormes/*/annotation/*/*.tsv; do
  if ! echo "${files}" | grep -q "panC*"; then
    echo "${files}" | grep -o "panC*"
  fi
done

grep -E -w -R "panC.*" ${workdir}03_MinIon_Tormes/*/annotation/*/*.tsv >> ${workdir}panCpositives.txt
grep -E -w -R "ces.*" ${workdir}03_MinIon_Tormes/*/annotation/*/*.tsv >> ${workdir}ces_positives.txt
grep -E -w -R "Nhe.*" ${workdir}03_MinIon_Tormes/*/annotation/*/*.tsv >> ${workdir}Nhe_positives.txt
grep -E -w -R "Hbl.*" ${workdir}03_MinIon_Tormes/*/annotation/*/*.tsv >> ${workdir}Hbl_positives.txt
grep -E -w -R "cytK.*" ${workdir}03_MinIon_Tormes/*/annotation/*/*.tsv >> ${workdir}cytK_positives.txt
