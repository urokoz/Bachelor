#!/bin/bash

pep_ker(){

rm $2

cat $1 | while read line
do
    declare -a fields=($line)

    echo ${fields[0]} ${fields[1]} $(cat $3${fields[0]}.txt | ./pep2score_db_kernel -kmin 3 -kmax 12 -blf ./Matrices/BLOSUM50_new -blqij ./Matrices/blosum62.qij -t 2 -- $3${fields[1]}.txt | grep -v "#" | cut -d " " -f 4) >> $2

done
}

pep_ker $1 $2 $3
