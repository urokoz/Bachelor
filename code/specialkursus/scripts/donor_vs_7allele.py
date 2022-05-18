#!/usr/bin/env python

from argparse import ArgumentParser
import numpy as np
from func_file import *
import matplotlib.pyplot as plt
import pickle

#----- functions -------

#def pearsons_cc(y_est, y_true):
    #n = len(y_est)
    #sum_x = sum(y_est)
    #sum_y = sum(y_true)
    #sum_xy = sum([x*y for x,y in zip(y_est,y_true)])
    #sum_xx = sum([x*x for x in y_est])
    #sum_yy = sum([y*y for y in y_true])

    #r_xy = (n*sum_xy - sum_x*sum_y)/(math.sqrt(n*sum_xx-sum_x**2)*math.sqrt(n*sum_yy-sum_y**2))
    #return r_xy


def sim_scatterplot(x, y, plot_title, xlabel, ylabel,blocker=False):
    #PCC = pearsons_cc(x,y)
    #SCC, _ = st.spearmanr(x,y)
    fig, ax = plt.subplots()
    ax.scatter(x, y)
    ax.set_title(plot_title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.show()
    

SI_list = []
rank_list_7allele = {}
rank_list_donor = {}

parser = ArgumentParser(description="Preps and filters the data and peptides")
parser.add_argument("-f", action="store", dest="data_file", type=str, help="Raw datafile")

args = parser.parse_args()
data_file = args.data_file

if False:
    HLA_dict = load_pep_HLA_data("../data/NetMHCIIpan/ragweed_7allele_NetMHCIIpan.xls")
    HLA_dict = load_pep_HLA_data("../data/NetMHCIIpan/ragweed_donor1_NetMHCIIpan.xls", HLA_dict)
    HLA_dict = load_pep_HLA_data("../data/NetMHCIIpan/ragweed_donor2_NetMHCIIpan.xls", HLA_dict)
    donor_allele_dict = load_donor_HLA_alleles("../data/ragweed_donor_HLA.txt")
else:
    HLA_dict = load_pep_HLA_data("../data/NetMHCIIpan/tree_7allele_NetMHCIIpan.xls")
    HLA_dict = load_pep_HLA_data("../data/NetMHCIIpan/tree_donor_NetMHCIIpan.xls", HLA_dict)
    donor_allele_dict = load_donor_HLA_alleles("../data/tree_donor_HLA.txt")


seven_alleles = ["DRB1_0301", "DRB1_0701", "DRB1_1501", "DRB3_0101", "DRB3_0202", "DRB4_0101", "DRB5_0101"]


with open("../data/dicts/bg_HLA_dict.pkl", 'rb') as f:
    bg_dict = pickle.load(f)

infile = open(data_file)
for line in infile:
    donor, peptide, SI = line.split()

    if not donor_allele_dict.get(donor):
        continue

    best_7allele_rank = 100
    best_cor_7allele_rank = 100
    for allele in seven_alleles:
        rank = HLA_dict[peptide][allele][0]

        if rank < best_7allele_rank:
            best_7allele_rank = rank

        corrected_rank = percent_v_bg(rank, allele, bg_dict)

        if corrected_rank < best_cor_7allele_rank:
            best_cor_7allele_rank = corrected_rank


    best_donor_rank = 100
    best_cor_donor_rank = 100
    for allele in donor_allele_dict[donor]:
        rank = HLA_dict[peptide][allele][0]

        if rank < best_donor_rank:
            best_donor_rank = rank

        corrected_rank = percent_v_bg(rank, allele, bg_dict)

        if corrected_rank < best_cor_donor_rank:
            best_cor_donor_rank = corrected_rank

    print(donor, peptide, SI, best_7allele_rank, best_donor_rank, round(best_cor_7allele_rank, 3), round(best_cor_donor_rank, 3), sep="\t")


# xlabel = "Best donor rank"
# ylabel = "7 allele rank"
# plot_title = ylabel + " vs. " + xlabel
# x = rank_list_donor
# y = rank_list_7allele
##print(sim_scatterplot(SI, best_7allele_rank, "7 Allele vs. SI", "SI", "Rank"))
#print(sim_scatterplot(SI, best_donor_rank, "Donor rank vs. SI", "SI", "Donor Rank"))
