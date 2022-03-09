#!/usr/bin/env python

def load_pep_HLA_data(datafile, pep_HLA_dict = dict()):
    """ Reads predicted MHCII binding for different HLA allels and outputs
        it as a dict.

        output format:
        pep_HLA_dict[pep_name][HLA_allele] = [rank, core]
    """

    infile = open(datafile,"r")

    allele_names = infile.readline().split()
    print(allele_names)
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


def load_donor_HLA_alleles(donor_file):
    infile = open(donor_file, "r")

    HLA_dict = dict()
    for line in infile:
        donor, alleles = line.split()
        HLA_dict[donor] = alleles.split(",")

    return HLA_dict
