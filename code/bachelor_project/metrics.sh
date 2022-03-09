#!/bin/bash

WD=$1

#./pep_ker.sh ${WD}prepped_data/unfiltered_dataset.csv ${WD}metrics/pep_ker_scores.txt ${WD}peptides/peptides/

./calculate_metrics.py -df ${WD}prepped_data/log_unfiltered_dataset.csv -pf ${WD}peptides/peptides.txt
./calculate_metrics.py -df ${WD}prepped_data/log_filtered_dataset.csv -pf ${WD}peptides/peptides.txt
./calculate_metrics.py -df ${WD}prepped_data/unfiltered_dataset.csv -pf ${WD}peptides/peptides.txt
./calculate_metrics.py -df ${WD}prepped_data/filtered_dataset.csv -pf ${WD}peptides/peptides.txt
