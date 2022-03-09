#!/bin/bash

print_usage() {
    echo "Usage: ./make_peptide_files.sh -d PEP_FILE"
    echo "-p    PEP_FILE    File containing information about the peptides"
    echo "-d    dir         Directory to store the peptide files in"
}

PEP_FILE=""
dir=""

while getopts p:d: flag
do
    case "${flag}" in
        p) PEP_FILE=${OPTARG};;
        d) dir=${OPTARG};;
        *) print_usage
           exit 1 ;;
    esac
done

if [ -z $PEP_FILE ] || [ -z $dir ]
then
    print_usage
    exit 1
fi

cat $PEP_FILE | while read line
do
    declare -a fields=($line)

    echo ${fields[1]} > $dir${fields[0]}.txt

done
