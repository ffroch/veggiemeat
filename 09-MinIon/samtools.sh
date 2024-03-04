#to get the reads per contig use calls to draft file in medaka and run samtools

conda activate samtools
samtools idxstats calls_to_draft.bam




#for all samples
workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/02_MinIon_ProcessedData
mkdir ${workdir}/calls_to_draft_stats
for file in ${workdir}/flye/*/medaka/calls_to_draft.bam*; do
  cp ${file} $(dirname "$file")/..
done
for file in ${workdir}/flye/*/calls_to_draft.bam*; do
  fp=$(dirname "$file"); fl=$(basename "$file"); mv "$fp/$fl" "$fp"_"$fl"
done
for file in ${workdir}/flye/*calls_to_draft.bam; do
  fl=$(basename "$file")
 samtools idxstats ${file} | tee ${workdir}/calls_to_draft_stats/${fl}.txt
done
for file in ${workdir}/flye/*calls_to_draft.bam*; do
  rm ${file}
done


#get part of a sequence
workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data
mkdir ${workdir}/99_miscellaneous
samtools faidx ${workdir}/03_MinIon_Tormes/minion-tormes-20221220/genomes/66.fasta contig_1:1463358-1606264 >${workdir}/99_miscellaneous/extraseq_Blicheniformis.fasta



#coverage statistics
#do for all files in the folder and save the output in a file with the same prefix
workdir=/data/Unit_LMM/selberherr-group/roch/02_projects/roch_veggie_unknownjournal_2022/data/02_MinIon_ProcessedData
for file in ${workdir}/flye/*/medaka/calls_to_draft.bam*; do
  cp ${file} $(dirname "$file")/..
done
for file in ${workdir}/flye/*/calls_to_draft.bam*; do
  fp=$(dirname "$file"); fl=$(basename "$file"); mv "$fp/$fl" "$fp"_"$fl"
done
for file in ${workdir}/flye/*calls_to_draft.bam; do
  fl=$(basename "$file")
 samtools depth ${file} | tee ${workdir}/calls_to_draft_stats/depth_${fl}.txt
done
for file in ${workdir}/flye/*calls_to_draft.bam*; do
  rm ${file}
done




conda deactivate
