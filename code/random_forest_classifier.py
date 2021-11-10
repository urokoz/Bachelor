#!/usr/bin/env python

import numpy as np
import math
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn import model_selection
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt
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


def sim_scatterplot(x, y, plot_title, xlabel, ylabel):
    PCC = pearsons_cc(x,y)
    fig, ax = plt.subplots()
    ax.scatter(x,y)
    ax.set_title(plot_title + " PCC: %.3f" % PCC)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.show()

parser = ArgumentParser(description="Extracts useful data from data files.")
parser.add_argument("-train", action="store", dest="train_file", type=str, default="Data/ragweed/metrics/log_filtered_metrics.txt", help="File with data")
parser.add_argument("-test", action="store", dest="test_file", type=str, default="Data/birch/metrics/log_filtered_metrics.txt", help="File with data")

args = parser.parse_args()
train_file = args.train_file
test_file = args.test_file

train_data = np.loadtxt(train_file, delimiter=",", dtype = str)
test_data = np.loadtxt(test_file, delimiter=",", dtype = str)

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


train_names = train_data[:,0]
train_sims = train_data[:,1:].astype(float)

X_train = train_sims[:, [5,9,10,12,14,16,18,19,20,21,22]]
y_train = train_sims[:,0]

test_names = test_data[:,0]
test_sims = test_data[:,1:].astype(float)

X_test = test_sims[:, [5,9,10,12,14,16,18,19,20,21,22]]
y_test = test_sims[:,0]

# K-fold crossvalidation
# K = 8
# CV_outer = model_selection.KFold(K, shuffle=True)
# # CV_inner = model_selection.KFold(K, shuffle=True)
# K = len(y)
# CV_outer = model_selection.LeaveOneOut()
# CV_inner = model_selection.LeaveOneOut()

forest = RandomForestRegressor()

forest.fit(X_train, y_train)

y_est = forest.predict(X_test)


# mse = []
# ft = None
# y_true = []
# y_est = []

# for (kout, (train_index_out, test_index_out)) in enumerate(CV_outer.split(X, y)):
#     print("Outer fold: {0}/{1}".format(kout + 1, K))
#     # initialize outer CV fold
#     X_out_train = X[train_index_out, :]
#     X_out_test = X[test_index_out, :]
#     y_out_train = y[train_index_out]
#     y_out_test = y[test_index_out]
#
#     forest.fit(X_out_train, y_out_train)
#     if ft:
#         ft_im = np.vstack((ft_im, forest.feature_importances_))
#     else:
#         ft_im = forest.feature_importances_
#         ft = True
#
#     y_pred = forest.predict(X_out_test)
#     y_est.append(y_pred)
#     y_true.append(y_out_test)
#
#     mse.append(mean_squared_error(y_out_test, y_pred, squared=False))
#
# print(np.mean(ft_im,axis=0))
print(forest.feature_importances_)

xlabel = "Y est"
ylabel = "Y true"
plot_title = ylabel +" vs. " + xlabel

sim_scatterplot(y_est, y_test, plot_title, xlabel, ylabel)
# print(np.mean(mse))
