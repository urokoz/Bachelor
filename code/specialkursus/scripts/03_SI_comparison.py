import sys
import numpy as np
from func_file import *
import pickle
import matplotlib.pyplot as plt
import scipy.stats as st

species = sys.argv[1]

## donor_pep_pair_dict gives a list of peptide pairs that the donor has data for
#   donor --> [(pep_pair_1A, pep_pair_1B), ..., (pep_pair_nA, pep_pair_nB)]
with open("../data/dicts/donor_pep_pair_dict_{}.pkl".format(species), 'rb') as f:
    donor_pep_pair_dict = pickle.load(f)

## donor_pep_dict gives data on donors relation to peptide
#   donor, pep --> SI, best_rank, best_cor_rank, best_core, best_cor_core
with open("../data/dicts/donor_pep_dict.pkl", 'rb') as f:
    donor_pep_dict = pickle.load(f)

blosum50 = load_blosum50()

high_core = []
low_core = []
high_SI = []
low_SI = []

both_binder_high = []
both_binder_low = []

rank_split = 5
for donor, pep_pair_list in donor_pep_pair_dict.items():
    if donor not in donor_pep_dict:
        continue

    for pep_pair in pep_pair_list:
        pep_1 = pep_pair[0]
        pep_2 = pep_pair[1]

        SI_1 = sigmoid(float(donor_pep_dict[donor][pep_1][0]))
        core_1 = donor_pep_dict[donor][pep_1][4]
        rank_1 = float(donor_pep_dict[donor][pep_1][2])
        SI_2 = sigmoid(float(donor_pep_dict[donor][pep_2][0]))
        core_2 = donor_pep_dict[donor][pep_2][4]
        rank_2 = float(donor_pep_dict[donor][pep_2][2])

        if SI_1 > SI_2:
            SI_h = SI_1
            SI_l = SI_2
            core_h = core_1
            core_l = core_2
            rank_h = rank_1
            rank_l = rank_2
        else:
            SI_h = SI_2
            SI_l = SI_1
            core_h = core_2
            core_l = core_1
            rank_h = rank_2
            rank_l = rank_1

        high_SI.append(SI_h)
        low_SI.append(SI_l)

        if SI_h >= 0.5:
            core_id, core_bl = score_cores(core_h, core_l, blosum50)# , weighting = [3,1,1,3,1,3,1,1,3])

            if SI_l >= 0.4:
                high_core.append(core_bl)
                if rank_h < rank_split and rank_l < rank_split:
                    both_binder_high.append(core_bl)
            else:
                low_core.append(core_bl)
                if rank_h < rank_split and rank_l < rank_split:
                    both_binder_low.append(core_bl)


print("Low both binder:", len(both_binder_low))
print("High both binder:", len(both_binder_high))


# Boxplot of core similarities for SI_h >= 0.5 and SI_l higher or lower than 0.4
fig, ax = plt.subplots()
p_val = st.ttest_ind(high_core,low_core, equal_var=False)[1]
ax.boxplot([high_core, low_core], vert = 0)
ax.set_yticklabels(["high SI", "low SI"])
ax.set_xlabel("core similarity")
ax.set_title("core similarity for high and low p-val = %.10f" % p_val)
plt.show(block = False)


# Boxplot of core similarities for both binder and not both binder
fig, ax = plt.subplots()
p_val = st.ttest_ind(both_binder_high, both_binder_low, equal_var=False)[1]
ax.boxplot([both_binder_high, both_binder_low], vert = 0)
ax.set_yticklabels(["both binder high SI", "both binder low SI"])
ax.set_xlabel("core similarity")
ax.set_title("core similarity for high and low p-val = %.10f" % p_val)
plt.show(block = False)


#         high.append(SI_h)
#         low.append(SI_l)
#
fig, ax = plt.subplots()
ax.scatter(high_SI, low_SI)
plt.show()


#set threshold: 0.4 (y-axis) and 0.5 (x-axis)
