#!/usr/bin/env python


def load_donor_HLA_alleles(donor_file):
    infile = open(donor_file, "r")

    HLA_dict = dict()
    for line in infile:
        donor, alleles = line.split()
        HLA_dict[donor] = alleles.split(",")

    return HLA_dict
