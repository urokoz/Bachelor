#!/usr/bin/env python
import numpy as np
import glob
import pickle
from func_file import *

## SETUP BACKGROUND DICT
print("Setting up background dict...")
paths = glob.glob("../data/NetMHCIIpan/results/*")

bg_HLA_dict = dict()
n_paths = len(paths)
for i, path in enumerate(paths):
    print("\tReading file {}/{}\t '{}'".format(i+1, n_paths, path), end="\r")
    bg_HLA_dict = load_bg_HLA_data(path, bg_HLA_dict)
print()

with open('../data/dicts/bg_HLA_dict.pkl', 'wb') as f:
    pickle.dump(bg_HLA_dict, f)


## SETUP PEPTIDE HLA BINDING DICT
print("Setting up peptide HLA binding dict...")
pep_HLA_dict = load_pep_HLA_data("../data/NetMHCIIpan/ragweed_donor1_NetMHCIIpan.xls")
pep_HLA_dict = load_pep_HLA_data("../data/NetMHCIIpan/ragweed_donor2_NetMHCIIpan.xls", pep_HLA_dict)
pep_HLA_dict = load_pep_HLA_data("../data/NetMHCIIpan/tree_donor_NetMHCIIpan.xls", pep_HLA_dict)

with open('../data/dicts/pep_HLA_dict.pkl', 'wb') as f:
    pickle.dump(pep_HLA_dict, f)

## SETUP DONOR ALLELE DICT
print("Setting up donor allele dict dict...")
donor_allele_dict = load_donor_HLA_alleles("../data/ragweed_donor_HLA.txt")
donor_allele_dict = load_donor_HLA_alleles("../data/tree_donor_HLA.txt", donor_allele_dict)

with open('../data/dicts/donor_allele_dict.pkl', 'wb') as f:
    pickle.dump(donor_allele_dict, f)

## SETUP DONOR PEPTIDE DICT
print("Setting up donor peptide dict...")
donor_pep_dict = load_donor_pep_dict("../data/donor_react_ragweed.tab", donor_allele_dict, pep_HLA_dict, bg_HLA_dict)
donor_pep_dict = load_donor_pep_dict("../data/donor_react_tree.tab", donor_allele_dict, pep_HLA_dict, bg_HLA_dict, donor_pep_dict)

with open('../data/dicts/donor_pep_dict.pkl', 'wb') as f:
    pickle.dump(donor_pep_dict, f)

print("Finished setting up dicts!")
