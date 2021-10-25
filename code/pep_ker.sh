#!/bin/bash

rm pep_ker_scores.txt

cat Data/filtered_dataset.csv | while read line
do
    declare -a fields=($line)

    cat peptides/${fields[0]}.txt | ./pep2score_db_kernel -kmin 3 -kmax 12 -blf ./Matrices/BLOSUM50_new -blqij ./Matrices/blosum62.qij -t 2 -- peptides/${fields[1]}.txt | grep -v "#" >> pep_ker_scores.txt

done

rm Data/calculated_metrics_2.txt

cut -d " " -f 4 pep_ker_scores.txt | paste -d "," Data/calculated_metrics.txt - >> Data/calculated_metrics_2.txt
