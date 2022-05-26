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

with open("../data/dicts/pep_seq_dict_{}.pkl".format(species), 'rb') as f:
    pep_seq_dict = pickle.load(f)

blosum50 = load_blosum50()

high_core = []
low_core = []
none_core = []
high_SI = []
low_SI = []

both_binder_high = []
both_binder_low = []
both_binder_none = []

print("donor", "pep_h", "pep_l", "SI_h", "SI_l",
      "bvc_core_h", "bvc_core_l", "bvc_rank_h", "bvc_rank_l", "bvc_core_id", "bvc_core_bl",
      "bvb_core_h", "bvb_core_l", "bvb_rank_h", "bvb_rank_l", "bvb_core_id", "bvb_core_bl",
      "rand_core_id", "rand_core_bl", sep="\t")

rank_split = 5
for donor, pep_pair_list in donor_pep_pair_dict.items():
    if donor not in donor_pep_dict:
        continue

    for pep_pair in pep_pair_list:
        pep_1 = pep_pair[0]
        pep_2 = pep_pair[1]

        SI_1 = sigmoid(float(donor_pep_dict[donor][pep_1][0]))
        rank_allele_list_1 = donor_pep_dict[donor][pep_1][1]

        SI_2 = sigmoid(float(donor_pep_dict[donor][pep_2][0]))
        rank_allele_list_2 = donor_pep_dict[donor][pep_2][1]

        rand_core_1 = random_9mer(pep_seq_dict[pep_1])
        rand_core_2 = random_9mer(pep_seq_dict[pep_2])

        rand_core_id, rand_core_bl = score_cores(rand_core_1, rand_core_2, blosum50)

        ## Best vs best core
        bvb_rank_1 = 100
        for allele, rank, corr_rank, core in rank_allele_list_1:
            # if rank < bvb_rank_1:
            #     bvb_rank_1 = rank
            #     bvb_core_1 = core

            if corr_rank < bvb_rank_1:
                bvb_rank_1 = corr_rank
                bvb_core_1 = core

        bvb_rank_2 = 100
        for allele, rank, corr_rank, core in rank_allele_list_2:
            # if rank < bvb_rank_2:
            #     bvb_rank_2 = rank
            #     bvb_core_2 = core

            if corr_rank < bvb_rank_2:
                bvb_rank_2 = corr_rank
                bvb_core_2 = core

        bvb_core_id, bvb_core_bl = score_cores(bvb_core_1, bvb_core_2, blosum50)# , weighting = [3,1,1,3,1,3,1,1,3])


        ## Best vs corresponding
        bvc_rank = 100
        bvc_collective_rank = 200
        for [allele1, rank1, corr_rank1, core1], [allele2, rank2, corr_rank2, core2] in zip(rank_allele_list_1, rank_allele_list_2):
            if min(corr_rank1, corr_rank2) < bvc_rank:
                bvc_rank = min(corr_rank1, corr_rank2)
                bvc_rank_1 = corr_rank1
                bvc_rank_2 = corr_rank2
                bvc_core_1 = core1
                bvc_core_2 = core2

            # if (corr_rank1 + corr_rank2) < bvc_collective_rank:
            #     bvc_collective_rank = corr_rank1 + corr_rank2
            #     bvc_rank_1 = corr_rank1
            #     bvc_rank_2 = corr_rank2
            #     bvc_core_1 = core1
            #     bvc_core_2 = core2

        bvc_core_id, bvc_core_bl = score_cores(bvc_core_1, bvc_core_2, blosum50)# , weighting = [3,1,1,3,1,3,1,1,3])


        if SI_1 > SI_2:
            pep_h = pep_1
            pep_l = pep_2
            SI_h = SI_1
            SI_l = SI_2

            bvc_core_h = bvc_core_1
            bvc_core_l = bvc_core_2
            bvc_rank_h = bvc_rank_1
            bvc_rank_l = bvc_rank_2

            bvb_core_h = bvb_core_1
            bvb_core_l = bvb_core_2
            bvb_rank_h = bvb_rank_1
            bvb_rank_l = bvb_rank_2
        else:
            pep_h = pep_2
            pep_l = pep_1
            SI_h = SI_2
            SI_l = SI_1

            bvc_core_h = bvc_core_2
            bvc_core_l = bvc_core_1
            bvc_rank_h = bvc_rank_2
            bvc_rank_l = bvc_rank_1

            bvb_core_h = bvb_core_2
            bvb_core_l = bvb_core_1
            bvb_rank_h = bvb_rank_2
            bvb_rank_l = bvb_rank_1

        print(donor, pep_h, pep_l, SI_h, SI_l,
              bvc_core_h, bvc_core_l, bvc_rank_h, bvc_rank_l, bvc_core_id, bvc_core_bl,
              bvb_core_h, bvb_core_l, bvb_rank_h, bvb_rank_l, bvb_core_id, bvb_core_bl,
              rand_core_id, rand_core_bl, sep="\t")


        # high_SI.append(SI_h)
        # low_SI.append(SI_l)
        #
        if SI_h >= 0.5:

            if SI_l >= 0.5:
                high_core.append(rand_core_bl)
                # if rank_h < rank_split and rank_l < rank_split:
                #     both_binder_high.append(rand_core_bl)
            else:
                low_core.append(rand_core_bl)
                # if rank_h < rank_split and rank_l < rank_split:
                #     both_binder_low.append(rand_core_bl)
        else:
            none_core.append(rand_core_bl)
            # if rank_h < rank_split and rank_l < rank_split:
            #     both_binder_none.append(rand_core_bl)


p_val = st.ttest_ind(high_core,low_core, equal_var=False)[1]
print("Rand p-val: {}".format(p_val), file=sys.stderr)


# print("Low-Low both binder:", len(both_binder_none))
# print("High-Low both binder:", len(both_binder_low))
# print("High-High both binder:", len(both_binder_high))
# print("Low-Low all:", len(none_core))
# print("High-Low all:", len(low_core))
# print("High-High all:", len(high_core))

# # Boxplot of core similarities for SI_h >= 0.5 and SI_l higher or lower than 0.4
# fig, ax = plt.subplots()
# p_val = st.ttest_ind(high_core,low_core, equal_var=False)[1]
# ax.boxplot([high_core, low_core], vert = 0)
# ax.set_yticklabels(["high SI", "low SI"])
# ax.set_xlabel("core similarity")
# ax.set_title("core similarity for high and low p-val = %.10f" % p_val)
# plt.show(block = False)
#
#
# # Boxplot of core similarities for both binder and not both binder
# fig, ax = plt.subplots()
# p_val = st.ttest_ind(both_binder_high, both_binder_low, equal_var=False)[1]
# ax.boxplot([both_binder_high, both_binder_low], vert = 0)
# ax.set_yticklabels(["both binder high SI", "both binder low SI"])
# ax.set_xlabel("core similarity")
# ax.set_title("core similarity for high and low p-val = %.10f" % p_val)
# plt.show(block = False)
#
#
# fig, ax = plt.subplots()
# ax.scatter(high_SI, low_SI)
# plt.show()


#set threshold: 0.4 (y-axis) and 0.5 (x-axis)
