#!/usr/bin/env python
# coding: utf-8

## PSSM
# Allel   -   binders   -   peptides    -   beta        -   SW        -   PCC

## SMM
# Allel   -   binders   -   peptides    -   lambda      -   gd/mc   -   PCC

## ANN
# Allel   -   bidners   -   peptides    -   n_hidden    -   eta     -   PCC


from argparse import ArgumentParser
import re

## Argument parsing:

parser = ArgumentParser(description="Extracts useful data from data files.")
parser.add_argument("-f", action="store", dest="data_file", type=str, help="File with data")

args = parser.parse_args()
data_file = args.data_file



infile = open(data_file, "r")

allel_stats_file = open("Allel_Statistics.txt", "r")
allel_stats = dict()

for line in allel_stats_file:
    if line.startswith("#"):
        continue

    stats = line.split()

    allel_stats[stats[0]] = [stats[1], stats[2]]


PSSM_list = []


for line in infile:
    PSSM = re.search(r"(\w+) PSSM beta: (\d+) SW: (\d) [\w ]+ N= (\d+) data: ([\d+-.]+) ([\d+-.e]+)", line)

    if PSSM:
        binders = allel_stats[PSSM.group(1)][1]
        PSSM_list.append([PSSM.group(1), binders, PSSM.group(4), PSSM.group(2), PSSM.group(3), PSSM.group(5), PSSM.group(6)])

for row in PSSM_list:
    print(row)
