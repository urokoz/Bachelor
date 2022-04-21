import numpy as np
from func_file import *
import pickle
import matplotlib.pyplot as plt

## donor_pep_pair_dict gives a list of peptide pairs that the donor has data for
#   donor --> [(pep_pair_1A, pep_pair_1B), ..., (pep_pair_nA, pep_pair_nB)]
with open("../data/dicts/donor_pep_pair_dict_ragweed.pkl", 'rb') as f:
    donor_pep_pair_dict = pickle.load(f)

## donor_pep_dict gives data on donors relation to peptide
#   donor, pep --> SI, best_rank, best_cor_rank, best_core, best_cor_core
with open("../data/dicts/donor_pep_dict.pkl", 'rb') as f:
    donor_pep_dict = pickle.load(f)

high = []
low = []
for donor, pep_pair_list in donor_pep_pair_dict.items():
    if donor not in donor_pep_dict:
        continue

    for pep_pair in pep_pair_list:
        pep_1 = pep_pair[0]
        pep_2 = pep_pair[1]

        SI_1 = sigmoid(float(donor_pep_dict[donor][pep_1][0]))
        SI_2 = sigmoid(float(donor_pep_dict[donor][pep_2][0]))

        if SI_1 > SI_2:
            SI_h = SI_1
            SI_l = SI_2
        else:
            SI_h = SI_2
            SI_l = SI_1

        high.append(SI_h)
        low.append(SI_l)

fig, ax = plt.subplots()
ax.scatter(high, low)

plt.show()


#set threshold: 0.4 (y-axis) and 0.5 (x-axis)
