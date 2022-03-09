#!/usr/bin/env python

from argparse import ArgumentParser
import numpy as np


parser = ArgumentParser(description="Preps and filters the data and peptides")
parser.add_argument("-f", action="store", dest="data_file", type=str, help="Raw datafile")


## Main

# Data format:
# ['5540', 'Amb', 'a', '1.0101', 'NSDKTIDGRGAKVEIINAGF', '3.74433',
#          'Amb', 'a', '1.0201', 'NSDKTIDGRGVKVNIVNAGL', '1.12407']

pep_list = []

pep_pair_dict = dict()
seen_peptide_pairs = set()
donor_reaction_dict = dict()
unique_seqs = set()
allergen_dict = dict()
pep_id_name = dict()


infile = open(data_file, "r")
infile.readline()   # remove header
for line in infile:
    line = line.replace("\n","").split("\t")


    donor_id = line[0]

    ori_pepseq = line[2]
    var_pepseq = line[5]

    ori_SI = float(line[3])
    var_SI = float(line[6])

    ori_name = line[1].strip().replace(" ", "_")
    var_name = line[4].strip().replace(" ", "_")

    ori_id = ori_name + "_" + ori_pepseq
    var_id = var_name + "_" + var_pepseq

    if log_switch:
        ori_SI = np.log(ori_SI)
        var_SI = np.log(var_SI)

    if frozenset([ori_id,var_id]) not in seen_peptide_pairs:
        seen_peptide_pairs.add(frozenset([ori_id,var_id]))

        if ori_id not in unique_seqs:
            unique_seqs.add(ori_id)

            if ori_name in allergen_dict:
                allergen_dict[ori_name] += 1
            else:
                allergen_dict[ori_name] = 1

            full_ori_name = ori_name + "_" + str(allergen_dict[ori_name])
            pep_list.append([full_ori_name, ori_pepseq])
            pep_id_name[ori_id] = full_ori_name

        if var_id not in unique_seqs:
            unique_seqs.add(var_id)

            if var_name in allergen_dict:
                allergen_dict[var_name] += 1
            else:
                allergen_dict[var_name] = 1

            full_var_name = var_name + "_" + str(allergen_dict[var_name])
            pep_list.append([full_var_name, var_pepseq])
            pep_id_name[var_id] = full_var_name

        # Save info about the peptide pairing and prep lists for SI values
        full_ori_name = pep_id_name[ori_id]
        full_var_name = pep_id_name[var_id]
        if full_ori_name < full_var_name:
            pep_pair = full_ori_name + "\t" + full_var_name
        else:
            pep_pair = full_var_name + "\t" + full_ori_name
        pep_pair_dict[pep_pair] = [[],[]]

    # add SI to the chart for the peptide pairing
    full_ori_name = pep_id_name[ori_id]
    full_var_name = pep_id_name[var_id]
    if pep_id_name[ori_id] < pep_id_name[var_id]:
        pep_pair = full_ori_name + "\t" + full_var_name
        pep_pair_dict[pep_pair][0].append(ori_SI)
        pep_pair_dict[pep_pair][1].append(var_SI)
    else:
        pep_pair = full_var_name + "\t" + full_ori_name
        pep_pair_dict[pep_pair][0].append(var_SI)
        pep_pair_dict[pep_pair][1].append(ori_SI)



    if donor_id not in donor_reaction_dict:
        donor_reaction_dict[donor_id] = dict()

    donor_reaction_dict[donor_id][ori_pepseq] = ori_SI
    donor_reaction_dict[donor_id][var_pepseq] = var_SI
infile.close()
