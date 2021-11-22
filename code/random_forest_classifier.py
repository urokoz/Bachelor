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


def ft_im_heatmap(ft_im, training_features_names, K, train_species):
    fig, ax = plt.subplots()
    c = plt.imshow(ft_im, interpolation='nearest', aspect = 1, vmin = 0)
    ax.set_title('Feature importance across crossvalidation folds - ' + train_species, fontsize=18)
    ax.set_ylabel("CV-fold", fontsize=12)
    ax.set_xticks(np.arange(len(training_features_names)))
    ax.set_xticklabels(training_features_names, fontsize=14)
    plt.xticks(rotation=90)
    ax.set_yticks(np.arange(0,K))
    ax.set_yticklabels(np.arange(1,K+1))
    plt.colorbar(c)
    fig.tight_layout()
    plt.show()


## New index overview
# 0. PCC
# 1. SCC
# 2. Needleman-Wunch naive sim
# 3. Needleman-Wunch naive score
# 4. Needleman-Wunch blosum sim
# 5. Needleman-Wunch blosum score   *
# 6. Smith-Waterman sim
# 7. Smith-Waterman blosum
# 8. K-mer identity
# 9. K-mer blosum   *
# 10. Pep kernel score  *
# 11. Best core vs. best core sim
# 12. Best core vs. best core blosum    *
# 13. Best ori core vs. corresponding var sim
# 14. Best ori core vs. corresponding var blosum    *
# 15. Best matching cores sim
# 16. Best matching cores blosum
# 17. Delta rank best core vs. best core    *
# 18. Pep 1 best rank
# 19. Pep 2 best rank
# 20. Pep 1 promiscuity *
# 21. Pep 2 promiscuity *
# 22. Binders in common *
# 23. nw_naive_sim x (100-delta_rank)
# 24. combined rank (1/rank1*1/rank2)
# 25. Promiscuity product (Pep 1 promiscuity * Pep 2 promiscuity)
# 26. nw_blosum + promiscuity product
# 27. bvc_and_prod_promis

all_features = ["PCC", "SCC", "Needleman-Wunch naive sim", "Needleman-Wunch naive score",
                "Needleman-Wunch blosum sim", "Needleman-Wunch blosum score",
                "Smith-Waterman sim", "Smith-Waterman blosum", "K-mer identity",
                "K-mer blosum", "Pep kernel score", "Best core vs. best core sim",
                "Best core vs. best core blosum", "Best core vs. corresponding sim",
                "Best core vs. corresponding blosum", "Best matching cores sim",
                "Best matching cores blosum", "Delta rank best core vs. best core",
                "Pep 1 best rank", "Pep 2 best rank", "Pep 1 promiscuity",
                "Pep 2 promiscuity", "Binders in common", "nw_naive_sim x (100-delta_rank)",
                "combined rank (1/rank1*1/rank2)", "Promiscuity product", "nw_blosum + promiscuity product",
                "bvc_and_prod_promis"]

parser = ArgumentParser(description="Extracts useful data from data files.")
parser.add_argument("-train", action="store", dest="train_species", type=str, default="birch", help="File with data")
parser.add_argument("-test", action="store", dest="test_species", type=str, default="ragweed", help="File with data")
parser.add_argument("-sf", action="store_true", default=False)

args = parser.parse_args()
train_species = args.train_species
test_species = args.test_species
single_file = args.sf

train_file = "Data/" + train_species + "/metrics/log_filtered_metrics.txt"
test_file = "Data/" + test_species + "/metrics/log_filtered_metrics.txt"

training_features = [10,14,17,25]
training_features_names = np.array(all_features)[training_features]

# K-fold crossvalidation
K = 5
CV = model_selection.KFold(K, shuffle=True)
# K = len(y)
# CV = model_selection.LeaveOneOut()

train_data = np.loadtxt(train_file, delimiter=",", dtype = str)

train_names = train_data[:,0]
train_sims = train_data[:,1:].astype(float)

X = train_sims[:, training_features]
y = train_sims[:, 0]

if not single_file:
    test_data = np.loadtxt(test_file, delimiter=",", dtype = str)

    test_names = test_data[:,0]
    test_sims = test_data[:,1:].astype(float)

    X_val = test_sims[:, training_features]
    y_val = test_sims[:, 0]

forest = RandomForestRegressor()

ft = None
y_true = []
y_est = []
y_val_est = []
y_val_true = []
for _ in range(1):
    for (k, (train_index, test_index)) in enumerate(CV.split(X, y)):
        print("Outer fold: {0}/{1}".format(k + 1, K))
        # initialize outer CV fold
        X_train = X[train_index, :]
        X_test = X[test_index, :]
        y_train = y[train_index]
        y_test = y[test_index]

        forest.fit(X_train, y_train)
        if ft:
            ft_im = np.vstack((ft_im, forest.feature_importances_))
        else:
            ft_im = forest.feature_importances_
            ft = True

        y_pred = forest.predict(X_test)
        y_est.extend(y_pred)
        y_true.extend(y_test)

        if not single_file:
            y_val_pred = forest.predict(X_val)
            y_val_est.extend(y_val_pred)
            y_val_true.extend(y_val)


pearson_training = pearsons_cc(y_true, y_est)
print("Pearson for training set(" + train_species + ") across CV folds:", pearson_training)
if not single_file:
    pearson_val = pearsons_cc(y_val_true, y_val_est)
    print("Pearson for validation set(" + test_species + ") across CV folds:", pearson_val)

ft_im_heatmap(ft_im, training_features_names, K, train_species)
