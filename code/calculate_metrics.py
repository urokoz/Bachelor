#!/usr/bin/env python

import numpy as np
from scipy.stats import spearmanr
import math

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


def needleman_wunsch(ori, var, match = 1, mismatch = -1, gap = -1, scoring_scheme = "naive"):
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
            if scoring_scheme == "blosum":
                diag = D[i,j] + blosum50[ori[i]][var[j]]
            else:
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


    matches = 0
    for j in range(len(ori_align)):
        if ori_align[j] == var_align[j]:
            matches += 1

    aligned_similarity = matches/min(n_ori,n_var)*100
    if scoring_scheme == "blosum":
        matches = max_score
    return ori_align, var_align, aligned_similarity, matches, max_score


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


def load_pep_dict(filename):
    infile = open(filename, "r")
    pep_dict = dict()
    for line in infile:
        line = line.split("\t")
        name = line[0]
        seq = line[1]
        HLA = [[float(pair.split(",")[0]),pair.split(",")[1]] for pair in line[2].split()]

        pep_dict[name] = [seq, HLA]
    infile.close()

    return pep_dict


def score_cores(core1, core2):
    core_matches = 0
    core_blosum = 0
    for a,b in zip(core1, core2):
        core_matches += int(a==b)
        core_blosum += blosum50[a][b]

        core_ident = core_matches/len(core1)*100

    return core_ident, core_blosum


def load_peptide_pair_significance(filename):
    infile = open(filename, "r")
    sig_list = [float(line.split()[-1]) for line in infile]
    infile.close()
    return sig_list


pep_dict = load_pep_dict("Data/filtered_pep_list.csv")

infile = open("Data/filtered_dataset.csv", "r")
outfile = open("Data/calculated_metrics.txt","w")

SRC_sig_list = load_peptide_pair_significance("Data/log_sampled_corr_SRC.txt")
PCC_sig_list = load_peptide_pair_significance("Data/log_sampled_corr_PCC.txt")

chart_sim_list = []

for line in infile:
    line = line.split()
    ori_name = line[0]
    var_name = line[1]
    ori_SI = [float(e) for e in line[2].split(",")]
    var_SI = [float(e) for e in line[3].split(",")]
    ori_pepseq = pep_dict[ori_name][0]
    var_pepseq = pep_dict[var_name][0]
    ori_HLA = pep_dict[ori_name][1]
    var_HLA = pep_dict[var_name][1]

    sim_list = [ori_name + " vs. " + var_name]

    # calculate corr for CR
    PCC = pearsons_cc(ori_SI, var_SI)
    SCC = spearmanr(ori_SI, var_SI)[0]
    sim_list.extend([PCC,SCC])

    # Global similarity with Needleman-Wunsch
    ori_align, var_align, nw_naive_sim, nw_naive_matches, nw_naive_score = needleman_wunsch(ori_pepseq, var_pepseq)
    ori_align, var_align, nw_blosum_sim, nw_blosum_matches, nw_blosum_score = needleman_wunsch(ori_pepseq, var_pepseq, scoring_scheme = "blosum")
    sim_list.extend([nw_naive_sim, nw_naive_score, nw_blosum_sim, nw_blosum_score])

    # local alignment with Smith-Waterman (O2)
    scoring_scheme = blosum50
    gap_open = -11
    gap_extension = -1
    P_matrix, Q_matrix, D_matrix, E_matrix, i_max, j_max, max_score = smith_waterman_alignment(ori_pepseq, var_pepseq, scoring_scheme, gap_open, gap_extension)
    aligned_query, aligned_database, sw_matches = smith_waterman_traceback(E_matrix, D_matrix, i_max, j_max, ori_pepseq, var_pepseq, gap_open, gap_extension)
    sw_sim = sw_matches/min(len(ori_pepseq), len(var_pepseq))*100
    sim_list.extend([sw_sim, max_score])

    ## Core similarity vs PCC
    # best core vs best core
    rank1 = np.inf
    rank2 = np.inf
    for [rank, core] in ori_HLA:
        if rank < rank1:
            rank1 = rank
            best_core1 = core

    for [rank, core] in var_HLA:
        if rank < rank2:
            rank2 = rank
            best_core2 = core

    bvb_ident, bvb_blosum = score_cores(best_core1, best_core2)

    # difference between best ranks
    delta_rank = np.abs(rank1 - rank2)

    # best core from ori and corresponding core for var
    best_rank = np.inf
    for [rank_ori, core_ori], [rank_var, core_var]  in zip(ori_HLA, var_HLA):
        if rank_ori < best_rank:
            best_rank = rank_ori
            var_rank = rank_var
            best_core3 = core_ori
            best_core4 = core_var

    bvc_ident, bvc_blosum = score_cores(best_core3, best_core4)

    bm_ident = None
    bm_blosum = None

    for [rank_ori, core_ori] in ori_HLA:
        if rank_ori > 5:
            continue
        for [rank_var, core_var] in var_HLA:
            if rank_var > 5:
                continue

            ident_match, blosum_match = score_cores(core_ori, core_var)

            if bm_blosum == None:
                bm_blosum = blosum_match
            elif blosum_match > bm_blosum:
                bm_blosum = blosum_match

            if bm_ident == None:
                bm_ident = ident_match
            elif ident_match > bm_ident:
                bm_ident = ident_match

    sim_list.extend([bvb_ident, bvb_blosum, bvc_ident, bvc_blosum, bm_ident, bm_blosum, delta_rank])

    # combination metrics
    glob_sim_and_d_rank = nw_naive_sim*(100-delta_rank)

    sim_list.extend([glob_sim_and_d_rank])

    chart_sim_list.append(sim_list)

    print(*sim_list, sep=",", file=outfile)
