from argparse import ArgumentParser
import numpy as np
from func_file import *
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


def sim_scatterplot(x, y, plot_title, xlabel, ylabel,blocker=False):
    #PCC = pearsons_cc(x,y)
    #SCC, _ = st.spearmanr(x,y)
    fig, ax = plt.subplots()
    ax.scatter(x, y)
    ax.set_title(plot_title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.show()



data = np.loadtxt("donor_data.txt", delimiter=" ", dtype = str)

#indicies:
#0 = donor nr.
#1 = allergen
#2 = SI
#3 = 7 allele rank
#4 = donor rank

xlabel = "SI"
ylabel = "Donor rank"
plot_title = ylabel + " vs. " + xlabel
x = data[:,2]
y = data[:,4]
sim_scatterplot(x, y, plot_title, xlabel, ylabel)
