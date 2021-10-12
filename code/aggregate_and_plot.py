#!/usr/bin/env python

import numpy as np
import math
import matplotlib.pyplot as plt


def pearsons_cc(y_est, y_true):
    """ Calucates the Pearson's correlation coefficient.
    Author: Mathias Rahbek-Borre, Fall 2020.
    """
    n = len(y_est)
    sum_x = sum(y_est)
    sum_y = sum(y_true)
    sum_xy = sum([x*y for x,y in zip(y_est,y_true)])
    sum_xx = sum([x*x for x in y_est])
    sum_yy = sum([y*y for y in y_true])

    r_xy = (n*sum_xy - sum_x*sum_y)/(math.sqrt(n*sum_xx-sum_x**2)*math.sqrt(n*sum_yy-sum_y**2))
    return r_xy


def sim_plot(x,y,plot_title, xlabel, ylabel):
    PCC = pearsons_cc(x,y)
    fig, ax = plt.subplots()
    ax.scatter(x,y)
    ax.set_title(plot_title + " PCC: %.3f" % PCC)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.show()


a = np.loadtxt("Data/calculated_metrics.txt", delimiter=",", dtype = str)

pep_pair_names = a[:,0]
pep_pair_sims = a[:,1:].astype(float)

xlabel = "Naive global similarity(%)"
ylabel = "SCC"
plot_title = ylabel +" vs. " + xlabel
x = pep_pair_sims[:,2]
y = pep_pair_sims[:,1]

sim_plot(x, y, plot_title, xlabel, ylabel)
