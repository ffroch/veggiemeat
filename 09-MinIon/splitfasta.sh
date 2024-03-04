#!/bin/bash

numseqs=$(grep -c ">" "$1");
numlines=$(wc -l < "$1");
if (( "$numlines" > $(( 2*$numseqs )) )); then
    echo "The fasta file needs to be linearised before this function will work.";
    return 1;
fi;

while read line; do
    if [ "${line:0:1}" == ">" ]; then
        header="$line";
        filename=$(echo "${line#>}" | sed 's/\ .*//g');
        touch "$filename".fasta
        echo "$header" >> "${filename}".fasta;
    else
        seq="$line";
        echo "$seq" >> "${filename}".fasta;
    fi;
done < $1
