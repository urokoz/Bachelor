#!/bin/bash

rm test1.txt
for i in {0..4}
do
    Data/pairlistscore_db_kernel Data/birch/prepped_data/test_$i Data/birch/prepped_data/train_$i -blf ./Matrices/BLOSUM50_new -blqij ./Matrices/blosum62.qij | grep -v "#" >> test1.txt
done

cut -d " " -f4,8 ./test1.txt | ./xycorr
