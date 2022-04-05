#!/usr/bin/env python
import pickle
import sys
import numpy as np


def percent_v_bg(rank, HLA_allele, bg_dict):
    return 100*np.mean([int(bg_rank < rank) for bg_rank in bg_dict[HLA_allele]])



with open(sys.argv[1], 'rb') as f:
    loaded_dict = pickle.load(f)

rank = 0.5
allele = "DRB1_1401"

print(percent_v_bg(rank, allele, loaded_dict))
