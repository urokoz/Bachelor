#!/bin/bash

# dataset
# transformation
# sig filtering
# HLA filtering

dataset=""

print_usage() {
  echo "Usage: ./data_prep.sh -d LOOKUP_FILE"
  echo "-d    LOOKUP_FILE   File containing information about the dataset"
}

while getopts d: flag
do
    case "${flag}" in
        d) dataset=${OPTARG};;
        *) print_usage
           exit 1 ;;
    esac
done

./dataload_and_filter.py -f $dataset
./dataload_and_filter.py -f $dataset -hs -ol 3
./dataload_and_filter.py -f $dataset -log
./dataload_and_filter.py -f $dataset -log -hs -ol 3
