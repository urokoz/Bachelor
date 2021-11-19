#!/bin/bash

species="ragweed"

rm ${species}_inherits_from_self.txt
for i in {0..4}
do
    Data/pairlistscore_db_kernel Data/${species}/prepped_data/test_$i Data/${species}/prepped_data/train_$i -blf ./Matrices/BLOSUM50_new -blqij ./Matrices/blosum62.qij | grep -v "#" >> ${species}_inherits_from_self.txt
done

cut -d " " -f4,8 ${species}_inherits_from_self.txt | ./spear_xycorr
