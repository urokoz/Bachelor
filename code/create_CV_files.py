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
parser.add_argument("-data", action="store", dest="data_file", type=str, default="Data/birch/metrics/log_filtered_metrics.txt", help="File with data")
parser.add_argument("-pep", action="store", dest="pep_file", type=str, default="Data/birch/peptides/peptides.txt", help="File with pepides")

args = parser.parse_args()
data_file = args.data_file
pep_file = args.pep_file

pep_dict = load_pep_dict(pep_file)

data = np.loadtxt(data_file, delimiter=",", dtype = str)

# K-fold crossvalidation
K = 5
CV_split = model_selection.KFold(K, shuffle=True)
# CV_inner = model_selection.KFold(K, shuffle=True)
# K = len(y)
# CV_outer = model_selection.LeaveOneOut()
# CV_inner = model_selection.LeaveOneOut()


for (k, (train_index, test_index)) in enumerate(CV_split.split(data)):
    print("Outer fold: {0}/{1}".format(k + 1, K))
    # initialize outer CV fold
    data_train = data[train_index, :]
    data_test = data[test_index, :]


    train_file_name = "Data/birch/prepped_data/train_{}".format(k)
    train_file = open(train_file_name, "w")
    for line in data_train:
        pep1_name = line[0].split()[0]
        pep2_name = line[0].split()[2]
        print(pep_dict[pep1_name][0], pep_dict[pep2_name][0], line[1], sep=" ", file=train_file)
    train_file.close()

    test_file_name = "Data/birch/prepped_data/test_{}".format(k)
    test_file = open(test_file_name, "w")
    for line in data_test:
        pep1_name = line[0].split()[0]
        pep2_name = line[0].split()[2]
        print(pep_dict[pep1_name][0], pep_dict[pep2_name][0], line[1], sep="\t", file=test_file)
    test_file.close()
