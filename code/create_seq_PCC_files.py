#!/usr/bin/env python

import numpy as np
from sklearn import model_selection
from argparse import ArgumentParser


def load_pep_dict(filename):
    infile = open(filename, "r")
    pep_dict = dict()
    for line in infile:
        line = line.split("\t")
        name = line[0]
        seq = line[1]
        HLA = [[float(pair.split(",")[0]),pair.split(",")[1]] for pair in line[2].split()]

        pep_dict[name] = [seq, HLA]
    infile.close()

    return pep_dict


parser = ArgumentParser(description="Extracts useful data from data files.")
parser.add_argument("-data", action="store", dest="data_file", type=str, default="Data/ragweed/metrics/log_filtered_metrics.txt", help="File with data")
parser.add_argument("-pep", action="store", dest="pep_file", type=str, default="Data/ragweed/peptides/peptides.txt", help="File with pepides")

args = parser.parse_args()
data_file = args.data_file
pep_file = args.pep_file

pep_dict = load_pep_dict(pep_file)

data = np.loadtxt(data_file, delimiter=",", dtype = str)

outfile_name = "Data/ragweed/prepped_data/ragweed_seq_PCC"
outfile = open(outfile_name, "w")
for line in data:
    pep1_name = line[0].split()[0]
    pep2_name = line[0].split()[2]
    print(pep_dict[pep1_name][0], pep_dict[pep2_name][0], line[1], sep=" ", file=outfile)
outfile.close()
