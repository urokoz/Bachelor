from argparse import ArgumentParser
import numpy as np
from func_file import *
import matplotlib.pyplot as plt
import scipy.stats as st
import pandas as pd

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
    SCC, _ = st.spearmanr(x,y)
    fig, ax = plt.subplots()
    ax.scatter(x, y)
    ax.set_title(plot_title + " SCC: {}".format(round(SCC,3)))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.show(block=blocker)


def barplot(x, y, plot_title, xlabel):
    fig, ax = plt.subplots()
    ax.bar(x, y)
    ax.set_title(plot_title)
    ax.set_xlabel(xlabel)
    plt.show()


#data = np.loadtxt("donor_data.txt", delimiter=" ", dtype = str)

data = pd.read_csv("donor_data.txt", sep='\s+', header = None)
data = data[data[1].str.match('Amb_a')]
data = data[data[2] >= 2]
data = data[data[3] <= 6]
data = data[data[4] <= 6]

#if False:
#    data[data == 0] = 0.01
#    data = np.log(data)
    # data[:,0] = np.log(data[:,0])

# indicies:
# 2 = SI
# 3 = 7 allele rank
# 4 = donor rank

xlabel = "Donor rank"
ylabel = "SI"
plot_title = ylabel + " vs. " + xlabel
y = np.asarray(data[2])
x = np.asarray(data[4])
sim_scatterplot(x, y, plot_title, xlabel, ylabel)
#
xlabel = "7 allele rank"
ylabel = "SI"
plot_title = ylabel + " vs. " + xlabel
y = np.asarray(data[2])
x = np.asarray(data[3])
sim_scatterplot(x, y, plot_title, xlabel, ylabel)
plt.show()
#
# xlabel = "Donor rank"
# ylabel = "7 allele rank"
# plot_title = ylabel + " vs. " + xlabel
# x = data[:,2]
# y = data[:,1]
# sim_scatterplot(x, y, plot_title, xlabel, ylabel,blocker=True)
#

#react = data[data[:,0] >= 2,:]
#non_react = data[data[:,0] < 2,:]

# print(np.sort(react[:,1]))
# print(np.sort(non_react[:,1]))

#fig, ax = plt.subplots()
#p_val = st.ttest_ind(react[:,2],non_react[:,2], equal_var=False)[1]
#ax.boxplot([non_react[:,2],react[:,2]], vert = 0)
#ax.set_yticklabels(["Non react", "React"])
#ax.set_xlabel("Donor rank")
#ax.set_title("Donor rank for CR(SI>=2) and non CR. p-val = %.10f" % p_val)
#fig, ax = plt.subplots()
#p_val = st.ttest_ind(react[:,1],non_react[:,1], equal_var=False)[1]
#ax.boxplot([non_react[:,1],react[:,1]], vert = 0)
#ax.set_yticklabels(["Non react", "React"])
#ax.set_xlabel("7allele rank")
#ax.set_title("7allele rank for CR(SI>=2) and non CR. p-val = %.10f" % p_val)
#plt.show()


#HLA_data = np.loadtxt("../results/donor/tree_HLA_count.txt", delimiter = "\t", dtype = str)

#fig, ax = plt.subplots()
#ax.bar(HLA_data[:,0], HLA_data[:,1].astype(float))
#ax.set_xlabel("HLA class II allele")
#ax.set_ylabel("Count")
#ax.set_title("HLA class II allele occurence in ragweed dataset")
#ax.set_xticklabels(HLA_data[:,0], rotation=45)
#ax.legend(fontsize=5, title_fontsize=15)
#plt.show()
#
