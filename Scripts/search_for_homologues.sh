#!/bin/bash

cd ./diamond_searches/output

ls *.output.txt > reads-files.txt

while read FILE; do

awk '{print $2}' $FILE > "$FILE"_seqname.txt

done < reads-files.txt

cd ../../mollicutes_data/proteomes

mkdir seq_name_files

mv ../../diamond_searches/output/*seqname.txt seq_name_files

cat seq_name_files/* > total_seq_names.txt
cat *.fa > fasta_file.fa

while read FILE_1; do
echo '>'$FILE_1 >>MHJ0450_homologues.fa
grep -A 100000 -w $FILE_1 fasta_file.fa | sed -n -e '1,/>/ {/>/ !{'p''}} >>MHJ0450_homologues.fa
done <total_seq_names.txt
 
