#!/usr/bin/env python

import numpy as np
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

# Data format:
# ['5540', 'Amb', 'a', '1.0101', 'NSDKTIDGRGAKVEIINAGF', '3.74433',
#          'Amb', 'a', '1.0201', 'NSDKTIDGRGVKVNIVNAGL', '1.12407']

i = -1
old_ori_seq = ""
old_var_seq = ""
charts = []
n = 0
seqs_for_FASTA = []
donor_reaction_dict = dict()
donor_list = []
wanted_charts = 1000


unique_seqs = set()
allergen_dict = dict()
pep_list = []
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
        charts.append([[],[], [pep_id_name[ori_id], pep_id_name[var_id]])

    # add SI to the chart for the peptide pairing
    charts[i][0].append(ori_SI)
    charts[i][1].append(var_SI)

    if donor_id not in donor_reaction_dict:
        donor_reaction_dict[donor_id] = dict()

    donor_reaction_dict[donor_id][ori_pepseq] = ori_SI
    donor_reaction_dict[donor_id][var_pepseq] = var_SI
infile.close()

outfile = open("filtered_dataset.csv")
for chart in charts
