import numpy as np
from func_file import *
import pickle
import matplotlib.pyplot as plt
import scipy.stats as st

## donor_pep_pair_dict gives a list of peptide pairs that the donor has data for
#   donor --> [(pep_pair_1A, pep_pair_1B), ..., (pep_pair_nA, pep_pair_nB)]
with open("../data/dicts/donor_pep_pair_dict_ragweed.pkl", 'rb') as f:
    donor_pep_pair_dict = pickle.load(f)

## donor_pep_dict gives data on donors relation to peptide
#   donor, pep --> SI, best_rank, best_cor_rank, best_core, best_cor_core
with open("../data/dicts/donor_pep_dict.pkl", 'rb') as f:
    donor_pep_dict = pickle.load(f)


#alphabet_file = alphabet_upload.values()
#alphabet_file = "https://raw.githubusercontent.com/brunoalvarez89/data/master/algorithms_in_bioinformatics/part_3/alphabet"
alphabet_file = "../data/Matrices/alphabet"
alphabet = np.loadtxt(alphabet_file, dtype=str)


#blosum_file = "https://raw.githubusercontent.com/brunoalvarez89/data/master/algorithms_in_bioinformatics/part_3/blosum50"
blosum_file = "../data/Matrices/BLOSUM50"
_blosum50 = np.loadtxt(blosum_file, dtype=float).reshape((24, -1)).T

blosum50 = {}

for i, letter_1 in enumerate(alphabet):

    blosum50[letter_1] = {}

    for j, letter_2 in enumerate(alphabet):

        blosum50[letter_1][letter_2] = _blosum50[i, j]


high = []
low = []
for donor, pep_pair_list in donor_pep_pair_dict.items():
    if donor not in donor_pep_dict:
        continue

    for pep_pair in pep_pair_list:
        pep_1 = pep_pair[0]
        pep_2 = pep_pair[1]

        SI_1 = sigmoid(float(donor_pep_dict[donor][pep_1][0]))
        core_1 = donor_pep_dict[donor][pep_1][4]
        SI_2 = sigmoid(float(donor_pep_dict[donor][pep_2][0]))
        core_2 = donor_pep_dict[donor][pep_2][4]

        if SI_1 > SI_2:
            SI_h = SI_1
            core_h = core_1
            SI_l = SI_2
            core_l = core_2
        else:
            SI_h = SI_2
            core_h = core_2
            SI_l = SI_1
            core_l = core_1

        if SI_h >= 0.5:
            core_id, core_bl = score_cores(core_h, core_l, blosum50)

            if SI_l >= 0.4:
                high.append(core_bl)
            else:
                low.append(core_bl)


fig, ax = plt.subplots()
p_val = st.ttest_ind(high,low, equal_var=False)[1]
ax.boxplot([high, low], vert = 0)
ax.set_yticklabels(["high", "low"])
ax.set_xlabel("core similarity")
ax.set_title("core similarity for high and low p-val = %.10f" % p_val)
plt.show()


#         high.append(SI_h)
#         low.append(SI_l)
#
# fig, ax = plt.subplots()
# ax.scatter(high, low)
#
# plt.show()


#set threshold: 0.4 (y-axis) and 0.5 (x-axis)
