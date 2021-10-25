#!/usr/bin/env python

import numpy as np
import math
import matplotlib.pyplot as plt


def pearsons_cc(y_est, y_true):
    n = len(y_est)
    sum_x = sum(y_est)
    sum_y = sum(y_true)
    sum_xy = sum([x*y for x,y in zip(y_est,y_true)])
    sum_xx = sum([x*x for x in y_est])
    sum_yy = sum([y*y for y in y_true])

    r_xy = (n*sum_xy - sum_x*sum_y)/(math.sqrt(n*sum_xx-sum_x**2)*math.sqrt(n*sum_yy-sum_y**2))
    return r_xy


def sim_scatterplot(x, y, plot_title, xlabel, ylabel):
    PCC = pearsons_cc(x,y)
    fig, ax = plt.subplots()
    ax.scatter(x,y)
    ax.set_title(plot_title + " PCC: %.3f" % PCC)
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

data = np.loadtxt("Data/calculated_metrics_2.txt", delimiter=",", dtype = str)

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


pep_pair_names = data[:,0]
pep_pair_sims = data[:,1:].astype(float)

# xlabel = "Naive global similarity(%)"
# ylabel = "SRC"
# plot_title = ylabel +" vs. " + xlabel
# x = pep_pair_sims[:,2]
# y = pep_pair_sims[:,1]
# sim_scatterplot(x, y, plot_title, xlabel, ylabel)
#
# xlabel = "Naive global similarity(%)"
# ylabel = "PCC"
# plot_title = ylabel + " vs. " + xlabel
# x = pep_pair_sims[:,2]
# y = pep_pair_sims[:,0]
# sim_scatterplot(x, y, plot_title, xlabel, ylabel)
#
# xlabel = "Local similarity(%)"
# ylabel = "SRC"
# plot_title = ylabel +" vs. " + xlabel
# x = pep_pair_sims[:,6]
# y = pep_pair_sims[:,1]
# sim_scatterplot(x, y, plot_title, xlabel, ylabel)
#
xlabel = "9-mere (% Identity)"
ylabel = "SRC"
plot_title = ylabel +" vs. " + xlabel
x = pep_pair_sims[:,8]
y = pep_pair_sims[:,0]
sim_scatterplot(x, y, plot_title, xlabel, ylabel)

xlabel = "9-mere (Blosum50)"
ylabel = "SRC"
plot_title = ylabel +" vs. " + xlabel
x = pep_pair_sims[:,9]
y = pep_pair_sims[:,0]
sim_scatterplot(x, y, plot_title, xlabel, ylabel)
#
# xlabel = "Histogram"
# ylabel = "Count"
# plot_title = "Histogram"
# x = pep_pair_sims[:,8]
# hist(x, plot_title, xlabel, ylabel)

xlabel = "Pep kernel score"
ylabel = "SCC"
plot_title = ylabel + " vs. " + xlabel
x = pep_pair_sims[:,18]
y = pep_pair_sims[:,1]
sim_scatterplot(x, y, plot_title, xlabel, ylabel)
