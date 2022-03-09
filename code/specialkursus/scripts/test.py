#!/usr/bin/env python

#!/usr/bin/env python

from argparse import ArgumentParser
import numpy as np
from func_file import *


parser = ArgumentParser(description="Preps and filters the data and peptides")
parser.add_argument("-f", action="store", dest="data_file", type=str, help="Raw datafile")

args = parser.parse_args()
data_file = args.data_file

HLA_dict = load_pep_HLA_data("../data/NetMHCIIpan/ragweed_7allele_NetMHCIIpan.xls")
HLA_dict = load_pep_HLA_data("../data/NetMHCIIpan/ragweed_donor1_NetMHCIIpan.xls", HLA_dict)
HLA_dict = load_pep_HLA_data("../data/NetMHCIIpan/ragweed_donor2_NetMHCIIpan.xls", HLA_dict)

donor_allele_dict = load_donor_HLA_alleles("../data/ragweed_donor_HLA.txt")

seven_alleles = ["DRB1_0301", "DRB1_0701", "DRB1_1501", "DRB3_0101", "DRB3_0202", "DRB4_0101", "DRB5_0101"]

infile = open(data_file)

for line in infile:
    donor, peptide, SI = line.split()

    if not donor_allele_dict.get(donor):
        continue

    best_7allele_rank = 100
    for allele in seven_alleles:
        rank = HLA_dict[peptide][allele][0]

        if rank < best_7allele_rank:
            best_7allele_rank = rank

    best_donor_rank = 100
    for allele in donor_allele_dict[donor]:
        rank = HLA_dict[peptide][allele][0]

        if rank < best_donor_rank:
            best_donor_rank = rank

    print(donor, peptide, SI, best_7allele_rank, best_donor_rank)
