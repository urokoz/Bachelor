# -*- coding: utf-8 -*-
"""Correlation.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1EI9LpecrRCWa46R_lCSZyxMmxAEOu6nI
"""

#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.stats as st
import pandas as pd
from scipy.stats import spearmanr


#alphabet_file = alphabet_upload.values()
#alphabet_file = "https://raw.githubusercontent.com/brunoalvarez89/data/master/algorithms_in_bioinformatics/part_3/alphabet"
alphabet_file = "Matrices/alphabet"
alphabet = np.loadtxt(alphabet_file, dtype=str)


#blosum_file = "https://raw.githubusercontent.com/brunoalvarez89/data/master/algorithms_in_bioinformatics/part_3/blosum50"
blosum_file = "Matrices/BLOSUM50"
_blosum50 = np.loadtxt(blosum_file, dtype=float).reshape((24, -1)).T

blosum50 = {}

for i, letter_1 in enumerate(alphabet):

    blosum50[letter_1] = {}

    for j, letter_2 in enumerate(alphabet):

        blosum50[letter_1][letter_2] = _blosum50[i, j]


def smith_waterman_alignment(query="VLLP", database="VLILP", scoring_scheme={}, gap_open=-5, gap_extension=-1):

    # Matrix imensions
    M = len(query)
    N = len(database)

    # D matrix change to float
    D_matrix = np.zeros((M+1, N+1), int)

    # P matrix
    P_matrix = np.zeros((M+1, N+1), int)

    # Q matrix
    Q_matrix = np.zeros((M+1, N+1), int)

    # E matrix
    E_matrix = np.zeros((M+1, N+1), dtype=object)

    # Initialize matrices
    for i in range(M, 0, -1):
        # Here you might include  penalties for end gaps, i.e
        # alignment_matrix[i-1, N] = alignment_matrix[i, N] + gap_open
        D_matrix[i-1, N] = 0
        P_matrix[i-1, N] = 0
        Q_matrix[i-1, N] = 0
        E_matrix[i-1, N] = 0

    for j in range(N, 0, -1):
        # Here you might include  penalties for end gaps, i.e
        #alignment_matrix[M, j-1] = alignment_matrix[M, j] + gap_open
        D_matrix[M, j-1] = 0
        P_matrix[M, j-1] = 0
        Q_matrix[M, j-1] = 0
        E_matrix[M, j-1] = 0


    # Main loop
    D_matrix_max_score, D_matrix_i_max, D_matrix_j_max = -9, -9, -9
    for i in range(M-1, -1, -1):
        for j in range(N-1, -1, -1):

            # Q_matrix[i,j] entry
            gap_open_database = D_matrix[i+1,j] + gap_open
            gap_extension_database = Q_matrix[i+1,j] + gap_extension
            max_gap_database = max(gap_open_database, gap_extension_database)

            Q_matrix[i,j] = max_gap_database

            # P_matrix[i,j] entry
            gap_open_query = D_matrix[i,j+1] + gap_open
            gap_extension_query = P_matrix[i,j+1] + gap_extension
            max_gap_query = max(gap_open_query, gap_extension_query)

            P_matrix[i,j] = max_gap_query

            # D_matrix[i,j] entry
            diagonal_score = D_matrix[i+1,j+1] + scoring_scheme[query[i]][database[j]]

            # E_matrix[i,j] entry
            candidates = [(1, diagonal_score),
                          (2, gap_open_database),
                          (4, gap_open_query),
                          (3, gap_extension_database),
                          (5, gap_extension_query)]

            direction, max_score = max(candidates, key=lambda x: x[1])

            # check entry sign
            if max_score > 0:
                E_matrix[i,j] = direction
                D_matrix[i, j] = max_score
            else:
                E_matrix[i,j] = 0
                D_matrix[i, j] = 0

            # fetch global max score
            if max_score > D_matrix_max_score:
                D_matrix_max_score = max_score
                D_matrix_i_max = i
                D_matrix_j_max = j

    return P_matrix, Q_matrix, D_matrix, E_matrix, D_matrix_i_max, D_matrix_j_max, D_matrix_max_score


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


# ## Alignment Matrix Traceback

# ### Traceback

# In[12]:


def smith_waterman_traceback(E_matrix, D_matrix, i_max, j_max, query="VLLP", database="VLILP", gap_open=-5, gap_extension=-1):

    M = len(query)
    N = len(database)

    aligned_query = []
    aligned_database = []

    # start from max_i, max_j
    i, j = i_max, j_max
    matches = 0
    while i < M and j < N:

        # E[i,j] = 0, stop back tracking
        if E_matrix[i, j] == 0:
            break


        # E[i,j] = 1, match
        if E_matrix[i, j] == 1:
            aligned_query.append(query[i])
            aligned_database.append(database[j])
            if ( query[i] == database[j]):
                matches += 1
            i += 1
            j += 1


        # E[i,j] = 2, gap opening in database
        if E_matrix[i, j] == 2:
            aligned_database.append("-")
            aligned_query.append(query[i])
            i += 1


        # E[i,j] = 3, gap extension in database
        if E_matrix[i, j] == 3:

            count = i + 2
            score = D_matrix[count, j] + gap_open + gap_extension

            # Find length of gap
            while((score - D_matrix[i, j])*(score - D_matrix[i, j]) >= 0.00001):
                count += 1
                score = D_matrix[count, j] + gap_open + (count-i-1)*gap_extension

            for k in range(i, count):
                aligned_database.append("-")
                aligned_query.append(query[i])
                i += 1


        # E[i,j] = 4, gap opening in query
        if E_matrix[i, j] == 4:
            aligned_query.append("-")
            aligned_database.append(database[j])
            j += 1


        # E[i,j] = 5, gap extension in query
        if E_matrix[i, j] == 5:

            count = j + 2
            score = D_matrix[i, count] + gap_open + gap_extension

            # Find length of gap
            while((score - D_matrix[i, j])*(score - D_matrix[i, j]) >= 0.0001):
                count += 1
                score = D_matrix[i, count] + gap_open + (count-j-1)*gap_extension

            for k in range(j, count):
                aligned_query.append("-")
                aligned_database.append(database[j])
                j += 1


    return aligned_query, aligned_database, matches

def needleman_wunsch(ori, var, match = 1, mismatch = -1, gap = -1):
  """ Aligns two sequences using Needleman-Wunsch alignment.
      Outputs the two aligned sequences, number of matches and a percent
      similarity.

      Author: Mathias Rahbek-Borre, 6-9-2021
      Inspired by: https://gist.github.com/slowkow/06c6dba9180d013dfd82bec217d22eb5
  """

    n_ori = len(ori)
    n_var = len(var)

    D = np.zeros((n_ori+1, n_var+1))
    P = np.zeros((n_ori+1, n_var+1))

    D[:,0] = np.arange(0,-(n_ori+1),-1)
    D[0,:] = np.arange(0,-(n_var+1),-1)
    P[:,0] = 4
    P[0,:] = 3
    P[0,0] = 0

    for i in range(n_ori):
        for j in range(n_var):
            if ori[i] == var[j]:
                diag = D[i,j] + match
            else:
                diag = D[i,j] + mismatch

            up = D[i+1,j] + gap
            left = D[i,j+1] + gap

            max_score = max(diag,up,left)

            D[i+1,j+1] = max_score

            if diag == max_score:
                P[i+1,j+1] = 2
            elif up == max_score:
                P[i+1,j+1] = 3
            elif left == max_score:
                P[i+1,j+1] = 4

    # print(D)
    # print(P)

    # traceback
    ori_align = ""
    var_align = ""

    i = n_ori
    j = n_var

    while i > 0 or j > 0:
        # print(i,j,P[i,j])

        if P[i,j] == 2:
            i -= 1
            j -= 1
            ori_align = ori[i] + ori_align
            var_align = var[j] + var_align
        elif P[i,j] == 3:
            j -= 1
            ori_align = "-" + ori_align
            var_align = var[j] + var_align
        elif P[i,j] == 4:
            i -= 1
            ori_align = ori[i] + ori_align
            var_align = "-" + var_align

      aligned_similarity = 0
      matches = 0
      for j in range(len(ori_align)):
          if ori_align[j] == var_align[j]:
              aligned_similarity += 100/len(ori_align)
              matches += 1

      return ori_align, var_align, aligned_similarity, matches


## Main
infile = open("Data/ragweed_Tcell_pairwise.MNi.tab", "r")

infile.readline()


# ['5540', 'Amb', 'a', '1.0101', 'NSDKTIDGRGAKVEIINAGF', '3.74433',
#          'Amb', 'a', '1.0201', 'NSDKTIDGRGVKVNIVNAGL', '1.12407']

i = -1
old_ori_seq = ""
old_var_seq = ""
charts = []
wanted_charts = 100000
n = 0
seqs_for_FASTA = []

for line in infile:
    line = line.split()

    ori_pepseq = line[4]
    var_pepseq = line[9]

    ori_SI = float(line[5])
    var_SI = float(line[10])

    if ori_pepseq != old_ori_seq or var_pepseq != old_var_seq:
        old_ori_seq = ori_pepseq
        old_var_seq = var_pepseq

        ori_name = line[1] + " " + line[2] + " " + line[3]
        var_name = line[6] + " " + line[7] + " " + line[8]

        title = ori_pepseq + "(" + ori_name + ")" + " vs. " + var_pepseq + "(" + var_name + ")"

        seqs_for_FASTA.append([ori_name, ori_pepseq])
        seqs_for_FASTA.append([var_name, var_pepseq])

        ## Similarity measurements

        # Global alignment with Needleman-Wunsch

        ori_align, var_align, aligned_similarity, nw_matches = needleman_wunsch(ori_pepseq, var_pepseq)

        # print(ori_align)
        # print(var_align)
        # print()

        # Dumb version of percent identity, where early gaps ruin similarity
        # while len(ori_pepseq) > len(var_pepseq):
        #     var_pepseq = var_pepseq + "x"


        # while len(ori_pepseq) < len(var_pepseq):
        #     ori_pepseq = ori_pepseq + "x"

        # naive_similarity = 0
        # for j in range(len(ori_pepseq)):
        #     if ori_pepseq[j] == var_pepseq[j]:
        #         naive_similarity += 100/len(ori_pepseq)

        # local alignment? with Smith-Waterman (O2)
        i += 1
        if i > wanted_charts-1:
            break

        scoring_scheme = blosum50
        gap_open = -11
        gap_extension = -1

        P_matrix, Q_matrix, D_matrix, E_matrix, i_max, j_max, max_score = smith_waterman_alignment(ori_pepseq, var_pepseq, scoring_scheme, gap_open, gap_extension)
        aligned_query, aligned_database, matches = smith_waterman_traceback(E_matrix, D_matrix, i_max, j_max, ori_pepseq, var_pepseq, gap_open, gap_extension)

        # print("ALN", "Origin", len(ori_pepseq), "Variant", len(var_pepseq), len(aligned_query), matches, max_score)
        # print("QAL", i_max, ''.join(aligned_query))
        # print("DAL", j_max,''.join(aligned_database))
        # print("")


        charts.append([[],[], title, aligned_similarity, 0])

    charts[i][0].append(ori_SI)
    charts[i][1].append(var_SI)
    charts[i][4] += 1

# output sequences in fasta format for HLA profiling.



coef_sim_matrix = [[],[],[]]
cross_react_count = [[],[],[]]
for chart in charts:
    # fig = plt.figure()
    # ax=fig.add_axes([0,0,1,1])
    # ax.scatter(chart[0], chart[1])
    # ax.set_xlabel("Ori SI")
    # ax.set_ylabel("Var SI")
    # ax.set_title(chart[2])
    # plt.show()

    n_sim = chart[3]

    PCC = pearsons_cc(chart[0], chart[1])
    SRC, p = spearmanr(chart[0], chart[1])
    # print("{:<8} {:<12} {:<12} {:<10}".format("n = %.d" % chart[4], "PCC: %.3f" % PCC, "SRC: %.3f" % SRC, "N_sim: %.d " % chart[3]))

    coef_sim_matrix[0].append(PCC)
    coef_sim_matrix[1].append(SRC)
    coef_sim_matrix[2].append(n_sim)

    CR = 1 if PCC > 0.5 else 0

    if n_sim < 50:
        cross_react_count[0].append(CR)
    elif n_sim >= 80:
        cross_react_count[2].append(CR)
    else:
        cross_react_count[1].append(CR)


mean_CR = [np.mean(cross_react_count[0]), np.mean(cross_react_count[1]), np.mean(cross_react_count[2])]

fig = plt.figure()
ax = fig.add_axes([0,0,1,1])
ax.scatter(coef_sim_matrix[2],coef_sim_matrix[0])
plt.show()
sim_v_PCC_PCC = pearsons_cc(coef_sim_matrix[2],coef_sim_matrix[0])
print(sim_v_PCC_PCC)

fig = plt.figure()
ax = fig.add_axes([0,0,1,1])
ax.bar(["<50%", "50%-80%", ">=80%"], mean_CR)
plt.show()

#hey
