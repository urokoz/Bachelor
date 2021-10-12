#!/usr/bin/env python

import numpy as np
from scipy.stats import spearmanr, pearsonr
import matplotlib.pyplot as plt

def heatmap(pep_list, donor_list, donor_reaction_dict):
    # Heatmap generation
    donor_reaction_overview = np.zeros((len(donor_list), len(pep_list)))
    for i in range(len(pep_list)):
        for j in range(len(donor_list)):
            donor_reaction_overview[j,i] = donor_reaction_dict.get(donor_list[j]).get(pep_list[i][1],-1)

    fig, ax = plt.subplots()
    c = plt.imshow(donor_reaction_overview, interpolation='nearest', vmax=25, aspect = "auto")
    ax.set_title('Donor SI per peptide heatmap', fontsize=18)
    ax.set_xlabel("Peptides", fontsize=12)
    ax.set_ylabel("Donors", fontsize=12)
    plt.colorbar(c)
    plt.savefig("../../Figures/Heatmap.png", dpi=500, bbox_inches="tight")
    plt.show(block=False)


def load_pep_HLA_data(datafile="Data/2860_NetMHCIIpan.xls"):
    """ Reads HLA binding for 7 allels and outputs it as a dict.

        output format:
        pep_HLA_dict[pep_name] = [[rank1,core1],...,[rank7,core7]]
    """

    infile = open(datafile,"r")

    infile.readline()
    infile.readline()
    pep_HLA_dict = dict()
    old_pep = 0
    for i, line in enumerate(infile):
        line = line.split()

        cur_pep = line[2]
        if old_pep != cur_pep:
            # if old_pep:
            #     pep_HLA_dict[old_pep] = [2 if a<1 else 1 if a<5 else 0 for a in pep_HLA_dict[old_pep]]

            old_pep = cur_pep

        HLA_bind_rank = [[float(line[i]), line[i-2]] for i in range(6,25,3)]

        if cur_pep in pep_HLA_dict:
            pep_HLA_dict[cur_pep] = [old if old[0] < new[0] else new for old, new in zip(pep_HLA_dict[cur_pep], HLA_bind_rank)]
        else:
            pep_HLA_dict[cur_pep] = HLA_bind_rank

    return pep_HLA_dict


def load_peptide_pair_significance(filename):
    infile = open(filename, "r")
    sig_list = [float(line.split()[-1]) for line in infile]
    infile.close()
    return sig_list


## Main
lower_cutoff = False
bottom_sort_out = False
log_switch = True
outlier_sorting = 3

# Data format:
# ['5540', 'Amb', 'a', '1.0101', 'NSDKTIDGRGAKVEIINAGF', '3.74433',
#          'Amb', 'a', '1.0201', 'NSDKTIDGRGVKVNIVNAGL', '1.12407']

i = -1
n = 0
wanted_charts = 1000
old_ori_seq = ""
old_var_seq = ""
charts = []
seqs_for_FASTA = []
donor_list = []
pep_list = []

donor_reaction_dict = dict()
unique_seqs = set()
allergen_dict = dict()
pep_dict = dict()
pep_id_name = dict()


pep_HLA_dict = load_pep_HLA_data()

infile = open("Data/ragweed_Tcell_pairwise.MNi.tab", "r")
infile.readline()   # remove header
for line in infile:
    line = line.split()

    donor_id = line[0]
    if donor_id not in donor_list:
        donor_list.append(donor_id)

    ori_pepseq = line[4]
    var_pepseq = line[9]

    ori_SI = float(line[5])
    var_SI = float(line[10])

    if log_switch:
        ori_SI = np.log(ori_SI)
        var_SI = np.log(var_SI)

    if ori_pepseq != old_ori_seq or var_pepseq != old_var_seq:
        i += 1
        if i > wanted_charts-1:
            break
        n += 1
        old_ori_seq = ori_pepseq
        old_var_seq = var_pepseq

        ori_name = line[1] + "_" + line[2] + "_" + line[3]
        var_name = line[6] + "_" + line[7] + "_" + line[8]

        ori_id = ori_name + "_" + ori_pepseq
        var_id = var_name + "_" + var_pepseq

        if ori_id not in unique_seqs:
            unique_seqs.add(ori_id)

            if ori_name in allergen_dict:
                allergen_dict[ori_name] += 1
            else:
                allergen_dict[ori_name] = 1

            full_ori_name = ori_name + "_" + str(allergen_dict[ori_name])
            pep_list.append([full_ori_name, ori_pepseq, pep_HLA_dict[full_ori_name]])
            pep_id_name[ori_id] = full_ori_name

        if var_id not in unique_seqs:
            unique_seqs.add(var_id)

            if var_name in allergen_dict:
                allergen_dict[var_name] += 1
            else:
                allergen_dict[var_name] = 1

            full_var_name = var_name + "_" + str(allergen_dict[var_name])
            pep_list.append([full_var_name, var_pepseq, pep_HLA_dict[full_var_name]])
            pep_id_name[var_id] = full_var_name

        # Save info about the peptide pairing and prep lists for SI values
        charts.append([[],[], pep_id_name[ori_id], pep_id_name[var_id]])

    # add SI to the chart for the peptide pairing
    charts[i][0].append(ori_SI)
    charts[i][1].append(var_SI)

    if donor_id not in donor_reaction_dict:
        donor_reaction_dict[donor_id] = dict()

    donor_reaction_dict[donor_id][ori_pepseq] = ori_SI
    donor_reaction_dict[donor_id][var_pepseq] = var_SI
infile.close()

SRC_sig_list = load_peptide_pair_significance("Data/log_sampled_corr_SRC.txt")
PCC_sig_list = load_peptide_pair_significance("Data/log_sampled_corr_PCC.txt")

outfile = open("Data/filtered_dataset.csv","w")
for chart, PCC_sig, SCC_sig in zip(charts, PCC_sig_list, SRC_sig_list):
    ori_SI = chart[0]
    var_SI = chart[1]
    ori_name = chart[2]
    var_name = chart[3]

    PCC, PCC_p = pearsonr(ori_SI,var_SI)
    SCC, SCC_p = spearmanr(ori_SI,var_SI)

    # by SRC significance
    if outlier_sorting == 1 or outlier_sorting == 3:
        if SCC_sig > 0.05 and SCC > 0.5:
            continue

        elif SCC_sig < 0.05 and SCC < -0.25 and bottom_sort_out:
            continue
        elif SCC < -0.25 and lower_cutoff:
            SCC = -0.25
    # by PCC significance
    if outlier_sorting == 2 or outlier_sorting == 3:
        if PCC_sig > 0.05 and PCC > 0.5:
            continue

        elif PCC_sig < 0.05 and PCC < -0.25 and bottom_sort_out:
            continue
        elif PCC < -0.25 and lower_cutoff:
            PCC = -0.25

    ori_bind = bool(sum([rank <= 5 for [rank, core] in pep_HLA_dict[ori_name]]))
    var_bind = bool(sum([rank <= 5 for [rank, core] in pep_HLA_dict[var_name]]))

    if not (ori_bind and var_bind):
        continue

    print(ori_name, var_name, ",".join(str(e) for e in ori_SI), ",".join(str(e) for e in var_SI), sep="\t", file=outfile)
outfile.close()

outfile = open("Data/filtered_pep_list.csv", "w")
for pep in pep_list:
    pep_name = pep[0]
    pep_seq = pep[1]
    pep_HLA = pep[2]
    pep_bind = bool(sum([rank <= 5 for [rank, core] in pep_HLA]))

    if not pep_bind:
        continue

    print(pep_name, pep_seq, " ".join(",".join(str(f) for f in e) for e in pep_HLA), sep="\t", file = outfile)

    ### For file 2
    # # by SRC significance
    # if outlier_sorting == 1 or outlier_sorting == 3:
    #     if SRC_sig > 0.05 and SRC > 0.5:
    #         continue
    #
    #     elif SRC_sig < 0.05 and SRC < -0.25 and bottom_sort_out:
    #         continue
    #     elif SRC < -0.25 and lower_cutoff:
    #         SRC = -0.25
    # # by PCC significance
    # if outlier_sorting == 2 or outlier_sorting == 3:
    #     if PCC_sig > 0.05 and PCC > 0.5:
    #         continue
    #
    #     elif PCC_sig < 0.05 and PCC < -0.25 and bottom_sort_out:
    #         continue
    #     elif PCC < -0.25 and lower_cutoff:
    #         PCC = -0.25
