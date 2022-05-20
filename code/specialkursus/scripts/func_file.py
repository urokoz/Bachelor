#!/usr/bin/env python
import numpy as np
import math
import random


def load_pep_HLA_data(datafile, pep_HLA_dict = dict()):
    """ Reads predicted MHCII binding for different HLA allels and outputs
        it as a dict.

        output format:
        pep_HLA_dict[pep_name][HLA_allele] = [rank, core]
    """

    infile = open(datafile,"r")

    allele_names = infile.readline().split()
    infile.readline()

    for i, line in enumerate(infile):
        line = line.split()
        cur_pep = line[2]

        for i in range(6, len(line)-2, 3):
            rank = float(line[i])
            core = line[i-2]
            HLA_allele = allele_names[int((i-4)/3)]

            if pep_HLA_dict.get(cur_pep):
                if pep_HLA_dict[cur_pep].get(HLA_allele):
                    old_rank = pep_HLA_dict[cur_pep][HLA_allele][0]
                    if rank < old_rank:
                        pep_HLA_dict[cur_pep][HLA_allele] = [rank, core]
                else:
                    pep_HLA_dict[cur_pep][HLA_allele] = [rank, core]
            else:
                pep_HLA_dict[cur_pep] = dict()
                pep_HLA_dict[cur_pep][HLA_allele] = [rank, core]


    infile.close()

    return pep_HLA_dict


def load_donor_HLA_alleles(donor_file, HLA_dict = dict()):
    infile = open(donor_file, "r")


    for line in infile:
        donor, alleles = line.split()
        HLA_dict[donor] = alleles.split(",")

    return HLA_dict


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

    infile.close()

    return bg_HLA_dict


def load_donor_pep_dict(data_file, donor_allele_dict, pep_HLA_dict, bg_dict, donor_pep_dict = dict()):
    infile = open(data_file, "r")

    lines = infile.readlines()
    infile.close()

    n_lines = len(lines)
    for i, line in enumerate(lines):
        print("\tReading '{}' line {}/{}".format(data_file, i+1, n_lines), end="\r")
        donor, peptide, SI = line.split()

        if not donor_allele_dict.get(donor):
            continue

        rank_allele_list = []
        for allele in donor_allele_dict[donor]:
            rank, core = pep_HLA_dict[peptide][allele]
            corrected_rank = percent_v_bg(rank, allele, bg_dict)
            rank_allele_list.append([allele, rank, corrected_rank, core])

        if not donor_pep_dict.get(donor):
            donor_pep_dict[donor] = dict()

        donor_pep_dict[donor][peptide] = [SI, rank_allele_list]
    print()
    return donor_pep_dict


def score_cores(core1, core2, blosum50, weighting = None):
    if weighting:
        core_matches = 0
        core_blosum = 0
        for a, b, w in zip(core1, core2, weighting):
            core_matches += w*int(a==b)
            core_blosum += w*blosum50[a][b]

            core_ident = core_matches/len(core1)*100
    else:
        core_matches = 0
        core_blosum = 0
        for a,b in zip(core1, core2):
            core_matches += int(a==b)
            core_blosum += blosum50[a][b]

            core_ident = core_matches/len(core1)*100

    return core_ident, core_blosum


def load_blosum50():
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

    return blosum50


def percent_v_bg(rank, HLA_allele, bg_dict):
    return 100*np.mean([int(bg_rank < rank) for bg_rank in bg_dict[HLA_allele]])


def sigmoid(si):
    if si > 3.5:
        return 1 / (1 + np.exp(-(si-3.5)*0.5))
    else:
        return 1 / (1 + np.exp(-(si-3.5)*2))


def random_9mer(seq):
    idx = random.randint(0, len(seq)-9)
    return seq[idx:idx+9]
