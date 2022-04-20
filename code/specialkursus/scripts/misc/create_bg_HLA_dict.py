#!/usr/bin/env python
import glob
import pickle

def load_bg_HLA_data(datafile, bg_HLA_dict = dict()):
    """ Reads predicted MHCII binding for of 10000 random uniprot peptides for
        different HLA allels and outputs it as a dict.

        output format:
        bg_HLA_dict[HLA_allele] = [10000 random ranks]
    """

    infile = open(datafile,"r")

    allele_name = infile.readline().strip()
    if not bg_HLA_dict.get(allele_name):
        bg_HLA_dict[allele_name] = []

    infile.readline()
    old_pep = ""
    first_flag = True
    for line in infile:
        line = line.split()
        cur_pep = line[2]
        rank = float(line[6])

        if cur_pep == old_pep:
            best_rank = min(rank, best_rank)
        else:
            if first_flag:
                best_rank = rank
                old_pep = cur_pep
                first_flag = False
            else:
                bg_HLA_dict[allele_name].append(best_rank)
                best_rank = rank
                old_pep = cur_pep

    bg_HLA_dict[allele_name].append(best_rank)

    print(len(bg_HLA_dict[allele_name]))

    infile.close()

    return bg_HLA_dict

paths = glob.glob("/home/mathias/Bachelor/Bachelor/code/specialkursus/data/NetMHCIIpan/results/*")

bg_HLA_dict = dict()
for path in paths:
    bg_HLA_dict = load_bg_HLA_data(path, bg_HLA_dict)

with open('bg_HLA_dict.pkl', 'wb') as f:
    pickle.dump(bg_HLA_dict, f)
