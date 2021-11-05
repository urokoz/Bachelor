#!/usr/bin/env python

import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.stats as st
from argparse import ArgumentParser


def pearsons_cc(y_est, y_true):
    n = len(y_est)
    sum_x = sum(y_est)
    sum_y = sum(y_true)
    sum_xy = sum([x*y for x,y in zip(y_est,y_true)])
    sum_xx = sum([x*x for x in y_est])
    sum_yy = sum([y*y for y in y_true])

    r_xy = (n*sum_xy - sum_x*sum_y)/(math.sqrt(n*sum_xx-sum_x**2)*math.sqrt(n*sum_yy-sum_y**2))
    return r_xy


def sim_scatterplot(x, y, plot_title, xlabel, ylabel,blocker=True):
    PCC = pearsons_cc(x,y)
    fig, ax = plt.subplots()
    ax.scatter(x,y)
    ax.set_title(plot_title + " PCC: %.3f" % PCC)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.show(block=blocker)


def hist(x, plot_title, xlabel, ylabel, bar_names = False):
    q25, q75 = np.percentile(x, [.25, .75])
    bin_width = 2 * (q75 - q25) * len(x) ** (-1 / 3)
    bins = round((max(x) - min(x)) / bin_width)

    fig, ax = plt.subplots()
    ax.hist(x, density = False, bins=bins)
    ax.set_title(plot_title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.show()

def boxplot(x, y, plot_title, xlabel, ylabel1, ylabel2):

    fig, ax = plt.subplots()
    p_val = st.ttest_ind(x,y, equal_var=False)[1]
    ax.boxplot([x,y], vert = 0)
    ax.set_yticklabels([ylabel1, ylabel2])
    ax.set_xlabel(xlabel)
    ax.set_title(plot_title + " P-val = %.10f" % p_val)
    plt.show()

def load_peptide_pair_significance(filename):
    infile = open(filename, "r")
    sig_list = [float(line.split()[-1]) for line in infile]
    infile.close()
    return sig_list


parser = ArgumentParser(description="Extracts useful data from data files.")
parser.add_argument("-df", action="store", dest="data_file", type=str, default="Data/calculated_metrics_2.txt", help="File with data")

args = parser.parse_args()
data_file = args.data_file

data = np.loadtxt(data_file, delimiter=",", dtype = str)

#Index overview:
#0. PCC
#1. SCC
#2. Needleman-Wunch naive sim
#3. Needleman-Wunch naive score
#4. Needleman-Wunch blosum sim
#5. Needleman-Wunch blosum score
#6. Smith-Waterman sim
#7. Smith-Waterman blosum
#8. K-mer identity
#9. K-mer blosum
#10. Best core vs. best core sim
#11. Best core vs. best core blosum
#12. Best ori core vs. corresponding var sim
#13. Best ori core vs. corresponding var blosum
#14. Best matching cores sim
#15. Best matching cores blosum
#16. Delta rank best core vs. best core
#17. nw_naive_sim x (100-delta_rank)
#18. Pep kernel score

## New index overview
# 0. PCC
# 1. SCC
# 2. Needleman-Wunch naive sim
# 3. Needleman-Wunch naive score
# 4. Needleman-Wunch blosum sim
# 5. Needleman-Wunch blosum score
# 6. Smith-Waterman sim
# 7. Smith-Waterman blosum
# 8. K-mer identity
# 9. K-mer blosum
# 10. Pep kernel score
# 11. Best core vs. best core sim
# 12. Best core vs. best core blosum
# 13. Best ori core vs. corresponding var sim
# 14. Best ori core vs. corresponding var blosum
# 15. Best matching cores sim
# 16. Best matching cores blosum
# 17. Delta rank best core vs. best core
# 18. Pep 1 best rank
# 19. Pep 2 best rank
# 20. Pep 1 promiscuity
# 21. Pep 2 promiscuity
# 22. Binders in common
# 23. nw_naive_sim x (100-delta_rank)
# 24. combined rank (1/rank1*1/rank2)


pep_pair_names = data[:,0]
pep_pair_sims = data[:,1:].astype(float)



# Visualize sorting #############################################################
# charts = np.loadtxt("Data/log_filtered_dataset.csv", delimiter="\t", dtype = str)
# SRC_sig_list = np.array(load_peptide_pair_significance("Data/log_sampled_corr_SRC.txt"))
# PCC_sig_list = np.array(load_peptide_pair_significance("Data/log_sampled_corr_PCC.txt"))
#
# SCC_sig_index = SRC_sig_list > 0.05
# PCC_sig_index = PCC_sig_list > 0.05
# insig_index = np.add(SCC_sig_index, PCC_sig_index)
#
# CR_index = pep_pair_sims[:,0].astype(float) > 0.5
#
# pos_insig_index = insig_index & CR_index
# filtered_out_charts = charts[pos_insig_index,:]
# for i, (pep1,pep2,x,y) in enumerate(filtered_out_charts):
#     x = [float(a) for a in x.split(",")]
#     y = [float(a) for a in y.split(",")]
#
#     xlabel = "Ori log(SI)"
#     ylabel = "Var log(SI)"
#     plot_title = pep1 +" vs. " + pep2
#     sim_scatterplot(x, y, plot_title, xlabel, ylabel)
#
#     if i > 10:
#         break
#
# non_binder_index = np.any(pep_pair_sims == "None",axis=1)
#
# sort_out_index = np.add(pos_insig_index,non_binder_index)
#
# xlabel = "Naive global similarity(%)"
# ylabel = "PCC"
# plot_title = "Filtered " + ylabel + " vs. " + xlabel
# x = pep_pair_sims[:,2].astype(float)[np.invert(sort_out_index)]
# y = pep_pair_sims[:,0].astype(float)[np.invert(sort_out_index)]
# sim_scatterplot(x, y, plot_title, xlabel, ylabel, False)
#
# plot_title = "Unfiltered " + ylabel + " vs. " + xlabel
# x = pep_pair_sims[:,2].astype(float)
# y = pep_pair_sims[:,0].astype(float)
# sim_scatterplot(x, y, plot_title, xlabel, ylabel)
#################################################################################

# xlabel = "Naive global similarity(%)"
# ylabel = "SRC"
# plot_title = ylabel +" vs. " + xlabel
# x = pep_pair_sims[:,2]
# y = pep_pair_sims[:,1]
# sim_scatterplot(x, y, plot_title, xlabel, ylabel)

# xlabel = "Naive global similarity(%)"
# ylabel = "PCC"
# plot_title = ylabel + " vs. " + xlabel
# x = pep_pair_sims[:,2]
# y = pep_pair_sims[:,0]
# sim_scatterplot(x, y, plot_title, xlabel, ylabel)

xlabel = "% Similarity - Pairwise Global Alignment (BLOSUM50)"
ylabel = "SRC"
plot_title = ""
x = pep_pair_sims[:,4]
y = pep_pair_sims[:,0]
sim_scatterplot(x, y, plot_title, xlabel, ylabel)

xlabel = "% Similarity - Pairwise Local Alignment (BLOSUM50)"
ylabel = "SRC"
plot_title = ""
x = pep_pair_sims[:,7]
y = pep_pair_sims[:,0]
sim_scatterplot(x, y, plot_title, xlabel, ylabel)


xlabel = "9-mere % Identity"
ylabel = "SRC"
plot_title = ylabel +" vs. " + xlabel
x = pep_pair_sims[:,8]
y = pep_pair_sims[:,0]
sim_scatterplot(x, y, plot_title, xlabel, ylabel)

xlabel = "9-mere BLOSUM50 score"
ylabel = "SRC"
plot_title = ylabel +" vs. " + xlabel
x = pep_pair_sims[:,9]
y = pep_pair_sims[:,0]
sim_scatterplot(x, y, plot_title, xlabel, ylabel)


# xlabel = "Histogram"
# ylabel = "Count"
# plot_title = "Histogram"
# x = pep_pair_sims[:,8]
# hist(x, plot_title, xlabel, ylabel)
#
# xlabel = "Pep kernel score"
# ylabel = "SCC"
# plot_title = ylabel + " vs. " + xlabel
# x = pep_pair_sims[:,18]
# y = pep_pair_sims[:,1]
# sim_scatterplot(x, y, plot_title, xlabel, ylabel)

CR_delta_rank = pep_pair_sims[pep_pair_sims[:,0] > 0.5, 17]
NCR_delta_rank = pep_pair_sims[np.invert(pep_pair_sims[:,0] > 0.5), 17]

CR_prom = np.prod(pep_pair_sims[pep_pair_sims[:,0] > 0.5, 20:22], axis = 1)
NCR_prom = np.prod(pep_pair_sims[np.invert(pep_pair_sims[:,0] > 0.5), 20:22], axis = 1)
#print(CR_prom)
#print(NCR_prom)

xlabel = "Pep kernel score"
ylabel = "SRC"
plot_title = ylabel + " vs. " + xlabel
x = pep_pair_sims[:,10]
y = pep_pair_sims[:,0]
sim_scatterplot(x, y, plot_title, xlabel, ylabel)

xlabel = "Best core vs. best core (% Identity)"
ylabel = "SRC"
plot_title = ylabel + " vs. " + "(" + xlabel + ")"
x = pep_pair_sims[:,11]
y = pep_pair_sims[:,0]
sim_scatterplot(x, y, plot_title, xlabel, ylabel)

xlabel = "Best core vs. best core (BLOSUM50)"
ylabel = "SRC"
plot_title = ylabel + " vs. " + "(" + xlabel + ")"
x = pep_pair_sims[:,12]
y = pep_pair_sims[:,0]
sim_scatterplot(x, y, plot_title, xlabel, ylabel)

xlabel = "Best ori core vs. corresponding var sim"
ylabel = "SRC"
plot_title = ylabel + " vs. " + xlabel
x = pep_pair_sims[:,13]
y = pep_pair_sims[:,0]
sim_scatterplot(x, y, plot_title, xlabel, ylabel)

xlabel = "Best ori core vs. corresponding var (BLOSUM50)"
ylabel = "SRC"
plot_title = ylabel + " vs. " + xlabel
x = pep_pair_sims[:,14]
y = pep_pair_sims[:,0]
sim_scatterplot(x, y, plot_title, xlabel, ylabel)

xlabel = "Best matching cores (BLOSUM50)"
ylabel = "SRC"
plot_title = ylabel + " vs. " + xlabel
x = pep_pair_sims[:,15]
y = pep_pair_sims[:,0]
sim_scatterplot(x, y, plot_title, xlabel, ylabel)

xlabel = "Delta rank best core vs. best core"
ylabel = "SRC"
plot_title = ylabel + " vs. " + xlabel
x = pep_pair_sims[:,17]
y = pep_pair_sims[:,0]
sim_scatterplot(x, y, plot_title, xlabel, ylabel)

xlabel = "nw_naive_sim x (100-delta_rank)"
ylabel = "SRC"
plot_title = ylabel + " vs. " + xlabel
x = pep_pair_sims[:,22]
y = pep_pair_sims[:,0]
sim_scatterplot(x, y, plot_title, xlabel, ylabel)


fig, ax = plt.subplots()
p_val = st.ttest_ind(NCR_delta_rank,CR_delta_rank, equal_var=False)[1]
ax.boxplot([NCR_delta_rank,CR_delta_rank], vert = 0)
ax.set_yticklabels(["Non-CR", "CR"])
ax.set_xlabel("Delta rank")
ax.set_title("Delta rank for CR and non CR. p-val = %.10f" % p_val)
plt.show()



xlabel = "Promescuity"
ylabel1 = "Cross reactive"
ylabel2 = "None Cross reactive"
plot_title = "Promescuity score"
x = CR_prom
y = NCR_prom
boxplot(x, y, plot_title, xlabel, ylabel1, ylabel2)
mean1, sigma1 = np.mean(CR_prom), st.sem(CR_prom)
mean2, sigma2 = np.mean(NCR_prom), st.sem(NCR_prom)
conf_int_CR = st.norm.interval(alpha = 0.95, loc = mean1, scale = sigma1)
conf_int_NCR = st.norm.interval(alpha = 0.95, loc = mean2, scale = sigma2)
#print(np.mean(CR_prom))
#print(np.mean(NCR_prom))
#print(conf_int_CR)
#print(conf_int_NCR)
#print(len(CR_prom))
#print(len(NCR_prom))

reci_log = np.log(pep_pair_sims[:,23])

xlabel = "Reciprok rank"
ylabel = "SRC"
plot_title = ylabel + " vs. " + xlabel
x = reci_log
y = pep_pair_sims[:,0]
sim_scatterplot(x, y, plot_title, xlabel, ylabel)
