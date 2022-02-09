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


def sim_scatterplot(x, y, plot_title, xlabel, ylabel,blocker=False):
    PCC = pearsons_cc(x,y)
    SCC, _ = st.spearmanr(x,y)
    fig, ax = plt.subplots()
    ax.scatter(x, y)
    ax.set_title(plot_title + " PCC: {}".format(round(PCC,3)))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.show()


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


def boxplot(x, y, plot_title, xlabel, ylabel1, ylabel2, blocker=False):

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
parser.add_argument("-df", action="store", dest="data_file", type=str, default="Data/birch/metrics/log_filtered_metrics.txt", help="File with data")

args = parser.parse_args()
data_file = args.data_file

data = np.loadtxt(data_file, delimiter=",", dtype = str)

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
# 25. Promiscuity product (Pep 1 promiscuity * Pep 2 promiscuity)
# 26. nw_blosum + promiscuity product
# 27. bvc_and_prod_promis

pep_pair_names = data[:,0]
pep_pair_sims = data[:,1:].astype(float)

# Visualize sorting #############################################################
# charts = np.loadtxt("Data/ragweed/prepped_data/log_unfiltered_dataset.csv", delimiter="\t", dtype = str)
# SRC_sig_list = np.array(load_peptide_pair_significance("Data/ragweed/significances/ragweed_log_SCC_sig.txt"))
# PCC_sig_list = np.array(load_peptide_pair_significance("Data/ragweed/significances/ragweed_log_PCC_sig.txt"))

# SCC_sig_index = SRC_sig_list > 0.05
# PCC_sig_index = PCC_sig_list > 0.05
#
# insig_index = np.add(SCC_sig_index, PCC_sig_index)
#
# CR_index = pep_pair_sims[:,0].astype(float) > 0.5
#
# pos_insig_index = insig_index & CR_index
#
# for (pep1,pep2,x,y), SRC_p, PCC_p in zip(charts, SRC_sig_list, PCC_sig_list):
#     x = [float(a) for a in x.split(",")]
#     y = [float(a) for a in y.split(",")]
#     if (SRC_p > 0.05 or PCC_p > 0.05) and pearsons_cc(x,y) > 0.5:
#         print("PCC: ", PCC_p)
#         print("SRC: ", SRC_p)
#         xlabel = "Ori log(SI)"
#         ylabel = "Var log(SI)"
#         plot_title = pep1 +" vs. " + pep2
#         sim_scatterplot(x, y, plot_title, xlabel, ylabel, blocker=True)

# #
# non_binder_index = np.any(pep_pair_sims == "None",axis=1)
#
# sort_out_index = np.add(pos_insig_index,non_binder_index)
#
# # xlabel = "Naive global similarity(%)"
# # ylabel = "PCC"
# # plot_title = "Filtered " + ylabel + " vs. " + xlabel
# # x = pep_pair_sims[:,2].astype(float)[np.invert(sort_out_index)]
# # y = pep_pair_sims[:,0].astype(float)[np.invert(sort_out_index)]
# # sim_scatterplot(x, y, plot_title, xlabel, ylabel)
#
# # plot_title = "Unfiltered " + ylabel + " vs. " + xlabel
# x = pep_pair_sims[:,2].astype(float)
# y = pep_pair_sims[:,0].astype(float)
#
# x1 = x[np.invert(sort_out_index)]
# y1 = y[np.invert(sort_out_index)]
# x2 = x[sort_out_index]
# y2 = y[sort_out_index]
# # sim_scatterplot(x, y, plot_title, xlabel, ylabel, blocker=True)
#
# print(len(x))
# print(len(x1))
#
# PCC_before = pearsons_cc(x,y)
# PCC_after = pearsons_cc(x1,y1)
# fig, ax1 = plt.subplots()
# ax1.scatter(x1,y1)
# ax1.scatter(x2,y2)
# ax1.set_title("PCC vs. naive global similarity(%). " + "PCC = %.3f. After = %.3f" % (PCC_before, PCC_after))
# ax1.set_xlabel("Global similarity(%)")
# ax1.set_ylabel("PCC")
# plt.show()
#################################################################################

# xlabel = "Naive global alignment (% Identity)"
# ylabel = "PCC"
# plot_title = ylabel +" vs. " + xlabel
# x = pep_pair_sims[:,2]
# y = pep_pair_sims[:,0]
# sim_scatterplot(x, y, plot_title, xlabel, ylabel)
#
# xlabel = "Naive global alignment score"
# ylabel = "PCC"
# plot_title = ylabel + " vs. " + xlabel
# x = pep_pair_sims[:,3]
# y = pep_pair_sims[:,0]
# sim_scatterplot(x, y, plot_title, xlabel, ylabel)
#
# xlabel = "Global Alignment BLOSUM50 (% Identity)"
# ylabel = "PCC"
# plot_title = "Global BLOSUM50"
# x = pep_pair_sims[:,4]
# y = pep_pair_sims[:,0]
# sim_scatterplot(x, y, plot_title, xlabel, ylabel)
#
# xlabel = "Global Alignment BLOSUM50 score"
# ylabel = "PCC"
# plot_title = ylabel +" vs. " + xlabel
# x = pep_pair_sims[:,5]
# y = pep_pair_sims[:,0]
# sim_scatterplot(x, y, plot_title, xlabel, ylabel)
#
#
# xlabel = "Local Alignment BLOSUM50 (% Identity)"
# ylabel = "PCC"
# plot_title = ylabel +" vs. " + xlabel
# x = pep_pair_sims[:,6]
# y = pep_pair_sims[:,0]
# sim_scatterplot(x, y, plot_title, xlabel, ylabel)
#
# xlabel = "Local Alignment BLOSUM50 score"
# ylabel = "PCC"
# plot_title = ylabel +" vs. " + xlabel
# x = pep_pair_sims[:,7]
# y = pep_pair_sims[:,0]
# sim_scatterplot(x, y, plot_title, xlabel, ylabel, blocker = True)
#
#
# xlabel = "9-mere (% Identity)"
# ylabel = "PCC"
# plot_title = ylabel +" vs. " + xlabel
# x = pep_pair_sims[:,8]
# y = pep_pair_sims[:,0]
# sim_scatterplot(x, y, plot_title, xlabel, ylabel)
#
# xlabel = "9-mere BLOSUM50 score"
# ylabel = "PCC"
# plot_title = ylabel +" vs. " + xlabel
# x = pep_pair_sims[:,9]
# y = pep_pair_sims[:,0]
# sim_scatterplot(x, y, plot_title, xlabel, ylabel, blocker=True)
#
# # xlabel = "Histogram"
# # ylabel = "Count"
# # plot_title = "Histogram"
# # x = pep_pair_sims[:,8]
# # hist(x, plot_title, xlabel, ylabel)
#
# xlabel = "Pep kernel score"
# ylabel = "PCC"
# plot_title = ylabel + " vs. " + xlabel
# x = pep_pair_sims[:,10]
# y = pep_pair_sims[:,0]
# sim_scatterplot(x, y, plot_title, xlabel, ylabel)
#
# xlabel = "Best core vs. best core (% Identity)"
# ylabel = "PCC"
# plot_title = ylabel + " vs. " + "(" + xlabel + ")"
# x = pep_pair_sims[:,11]
# y = pep_pair_sims[:,0]
# sim_scatterplot(x, y, plot_title, xlabel, ylabel)
#
# xlabel = "Best core vs. best core (BLOSUM50)"
# ylabel = "PCC"
# plot_title = ylabel + " vs. " + "(" + xlabel + ")"
# x = pep_pair_sims[:,12]
# y = pep_pair_sims[:,0]
# sim_scatterplot(x, y, plot_title, xlabel, ylabel)
#
# xlabel = "Best core vs. corresponding (% Identity)"
# ylabel = "PCC"
# plot_title = ylabel + " vs. " + xlabel
# x = pep_pair_sims[:,13]
# y = pep_pair_sims[:,0]
# sim_scatterplot(x, y, plot_title, xlabel, ylabel)
#
# xlabel = "Best core vs. corresponding (BLOSUM50)"
# ylabel = "PCC"
# plot_title = ylabel + " vs. " + xlabel
# # x = pep_pair_sims[:,14]
# # y = pep_pair_sims[:,0]
# # Removing left side clustering:
# x = pep_pair_sims[pep_pair_sims[:,14] > 15,14]
# y = pep_pair_sims[pep_pair_sims[:,14] > 15,0]
# sim_scatterplot(x, y, plot_title, xlabel, ylabel, blocker=True)
#
# xlabel = "Best matching cores (% Identity)"
# ylabel = "PCC"
# plot_title = ylabel + " vs. " + xlabel
# x = pep_pair_sims[:,15]
# y = pep_pair_sims[:,0]
# sim_scatterplot(x, y, plot_title, xlabel, ylabel)
#
# xlabel = "Best matching cores (BLOSUM50)"
# ylabel = "PCC"
# plot_title = ylabel + " vs. " + xlabel
# x = pep_pair_sims[:,16]
# y = pep_pair_sims[:,0]
# sim_scatterplot(x, y, plot_title, xlabel, ylabel, blocker=True)
#
# xlabel = "Delta rank best core vs. best core"
# ylabel = "PCC"
# plot_title = ylabel + " vs. " + xlabel
# x = pep_pair_sims[:,17]
# y = pep_pair_sims[:,0]
# sim_scatterplot(x, y, plot_title, xlabel, ylabel)
#
# xlabel = "Naive global (% Identity) x (100-delta_rank)"
# ylabel = "PCC"
# plot_title = ylabel + " vs. " + xlabel
# x = pep_pair_sims[:,23]
# y = pep_pair_sims[:,0]
# sim_scatterplot(x, y, plot_title, xlabel, ylabel, blocker=False)
#
# xlabel = "Global blosum score + promiscuity score"
# ylabel = "PCC"
# plot_title = ylabel + " vs. " + xlabel
# x = pep_pair_sims[:,26]
# y = pep_pair_sims[:,0]
# sim_scatterplot(x, y, plot_title, xlabel, ylabel, blocker=False)
#
# xlabel = "Best vs. corresponding + promiscuity score"
# ylabel = "PCC"
# plot_title = ylabel + " vs. " + xlabel
# x = pep_pair_sims[:,27]
# y = pep_pair_sims[:,0]
# sim_scatterplot(x, y, plot_title, xlabel, ylabel, blocker=False)
#
#
#
#
CR_delta_rank = pep_pair_sims[pep_pair_sims[:,0] > 0.5, 17]
NCR_delta_rank = pep_pair_sims[np.invert(pep_pair_sims[:,0] > 0.5), 17]

CR_prom = np.prod(pep_pair_sims[pep_pair_sims[:,0] > 0.5, 20:22], axis = 1)
NCR_prom = np.prod(pep_pair_sims[np.invert(pep_pair_sims[:,0] > 0.5), 20:22], axis = 1)

fig, ax = plt.subplots()
p_val = st.ttest_ind(CR_delta_rank,NCR_delta_rank, equal_var=False)[1]
ax.boxplot([NCR_delta_rank,CR_delta_rank], vert = 0)
ax.set_yticklabels(["Non CR", "CR"])
ax.set_xlabel("Delta rank")
ax.set_title("Delta rank for CR and non CR. p-val = %.10f" % p_val)
#plt.show()

CR_binder = pep_pair_sims[pep_pair_sims[:,0] > 0.5, 22]
NCR_binder = pep_pair_sims[np.invert(pep_pair_sims[:,0] > 0.5), 22]

fig, ax = plt.subplots()
p_val = st.ttest_ind(NCR_binder,CR_binder, equal_var=False)[1]
ax.boxplot([CR_binder,NCR_binder], vert = 0)
ax.set_yticklabels(["Non CR", "CR"])
ax.set_xlabel("Binders in common")
ax.set_title("Binders in common for CR and non CR. p-val = %.10f" % p_val)
mean1, sigma1 = np.mean(CR_binder), st.sem(CR_binder)
mean2, sigma2 = np.mean(NCR_binder), st.sem(NCR_binder)
conf_int_CR = st.norm.interval(alpha = 0.95, loc = mean1, scale = sigma1)
conf_int_NCR = st.norm.interval(alpha = 0.95, loc = mean2, scale = sigma2)
print(np.mean(CR_binder))
print(np.mean(NCR_binder))
print(conf_int_CR)
print(conf_int_NCR)
print(len(CR_binder))
print(len(NCR_binder))
#plt.show()





xlabel = "Promescuity"
ylabel1 = "CR"
ylabel2 = "Non CR"
plot_title = "Promescuity score"
x = CR_prom
y = NCR_prom



fig, ax = plt.subplots()
p_val = st.ttest_ind(NCR_prom,CR_prom, equal_var=False)[1]
ax.boxplot([CR_prom,NCR_prom], vert = 0)
ax.set_yticklabels(["Non CR", "CR"])
ax.set_xlabel("Promescuity score")
ax.set_title("Promescuity score for CR and non CR. p-val = %.10f" % p_val)
plt.show()

mean1, sigma1 = np.mean(CR_prom), st.sem(CR_prom)
mean2, sigma2 = np.mean(NCR_prom), st.sem(NCR_prom)
conf_int_CR = st.norm.interval(alpha = 0.95, loc = mean1, scale = sigma1)
conf_int_NCR = st.norm.interval(alpha = 0.95, loc = mean2, scale = sigma2)
print(np.mean(CR_prom))
print(np.mean(NCR_prom))
#print(conf_int_CR)
#print(conf_int_NCR)
#print(len(CR_prom))
#print(len(NCR_prom))
#
# reci_log = np.log(pep_pair_sims[:,23])
#
# xlabel = "Reciprok rank"
# ylabel = "PCC"
# plot_title = ylabel + " vs. " + xlabel
# x = reci_log
# y = pep_pair_sims[:,0]
# sim_scatterplot(x, y, plot_title, xlabel, ylabel)
#
# #printing all performances (BLOSUM50 + PCC)
#
#Global
print("Global naive ident:", pearsons_cc(pep_pair_sims[:,2], pep_pair_sims[:,0]))
print("Global naive score:", pearsons_cc(pep_pair_sims[:,3], pep_pair_sims[:,0]))
print("Global blosum ident:", pearsons_cc(pep_pair_sims[:,4], pep_pair_sims[:,0]))
print("Global blosum score:", pearsons_cc(pep_pair_sims[:,5], pep_pair_sims[:,0]))
# #Local
print("Local blosum score:", pearsons_cc(pep_pair_sims[:,7], pep_pair_sims[:,0]))
# #9-mere
print("9-mer blosum score:", pearsons_cc(pep_pair_sims[:,9], pep_pair_sims[:,0]))
# #pep_kernal
print("Peptide kernel score:", pearsons_cc(pep_pair_sims[:,10], pep_pair_sims[:,0]))
# #best vs. best
print("Best v. best core ident:", pearsons_cc(pep_pair_sims[:,11], pep_pair_sims[:,0]))
print("Best v. best core blosum:", pearsons_cc(pep_pair_sims[:,12], pep_pair_sims[:,0]))
# #bestvs.corresponding
print("Best v. corresponding core ident:", pearsons_cc(pep_pair_sims[:,13], pep_pair_sims[:,0]))
print("Best v. corresponding core blosum:", pearsons_cc(pep_pair_sims[:,14], pep_pair_sims[:,0]))
# #best matching cores
print("Best matching binding cores ident:", pearsons_cc(pep_pair_sims[:,15], pep_pair_sims[:,0]))
print("Best matching binding cores blosum:", pearsons_cc(pep_pair_sims[:,16], pep_pair_sims[:,0]))
# #combined rank
print("Combined rank:", pearsons_cc(pep_pair_sims[:,24], pep_pair_sims[:,0]))
#naive_sim * (100-delta_rank)
print("Naive sim x delta_rank:", pearsons_cc(pep_pair_sims[:,23], pep_pair_sims[:,0]))
# #nw_blosum + prom
print("NW_blosum + prom:", pearsons_cc(pep_pair_sims[:,26], pep_pair_sims[:,0]))
