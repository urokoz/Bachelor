from argparse import ArgumentParser
import numpy as np
from func_file import *
import matplotlib.pyplot as plt
import scipy.stats as st

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
    plt.show(block=blocker)





data = np.loadtxt("../data/tree_donor_data.txt", delimiter=" ", dtype = str)

data = data[:,2:].astype(float)

# data[:,1:] = 100 - data[:,1:]

if False:
    data[data == 0] = 0.01
    data = np.log(data)
    # data[:,0] = np.log(data[:,0])

# indicies:
# 0 = SI
# 1 = 7 allele rank
# 2 = donor rank

# xlabel = "SI"
# ylabel = "Donor rank"
# plot_title = ylabel + " vs. " + xlabel
# x = data[:,0]
# y = data[:,2]
# sim_scatterplot(x, y, plot_title, xlabel, ylabel)
#
# xlabel = "SI"
# ylabel = "7 allele rank"
# plot_title = ylabel + " vs. " + xlabel
# x = data[:,0]
# y = data[:,1]
# sim_scatterplot(x, y, plot_title, xlabel, ylabel)
#
# xlabel = "Donor rank"
# ylabel = "7 allele rank"
# plot_title = ylabel + " vs. " + xlabel
# x = data[:,2]
# y = data[:,1]
# sim_scatterplot(x, y, plot_title, xlabel, ylabel,blocker=True)


react = data[data[:,0] >= 2,:]
non_react = data[data[:,0] < 2,:]


fig, ax = plt.subplots()
p_val = st.ttest_ind(react[:,2],non_react[:,2], equal_var=False)[1]
ax.boxplot([non_react[:,2],react[:,2]], vert = 0)
ax.set_yticklabels(["Non react", "React"])
ax.set_xlabel("Donor rank")
ax.set_title("Donor rank for CR(SI>=2) and non CR. p-val = %.10f" % p_val)
plt.show(block=False)


fig, ax = plt.subplots()
p_val = st.ttest_ind(react[:,1],non_react[:,1], equal_var=False)[1]
ax.boxplot([non_react[:,1],react[:,1]], vert = 0)
ax.set_yticklabels(["Non react", "React"])
ax.set_xlabel("7allele rank")
ax.set_title("7allele rank for CR(SI>=2) and non CR. p-val = %.10f" % p_val)
plt.show()
