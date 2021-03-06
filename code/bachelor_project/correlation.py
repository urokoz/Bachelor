#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.stats as st
from scipy.stats import spearmanr
import pandas as pd
from sklearn.neighbors import LocalOutlierFactor
import seaborn as sns

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

    return ori_align, var_align, aligned_similarity, matches


# def name_peptides(seqs):
#     unique_seqs = set()
#     allergen_dict = dict()
#     pep_list = []
#     for [name, seq] in seqs_for_FASTA:
#         if seq not in unique_seqs:
#             unique_seqs.add(seq)
#
#             if name in allergen_dict:
#                 allergen_dict[name] += 1
#             else:
#                 allergen_dict[name] = 1
#
#             pep_list.append([name + "_" + str(allergen_dict[name]), seq])
#     return pep_list


def heatmap(pep_list, donor_list, donor_reaction_dict):
    # Heatmap generation
    donor_reaction_overview = np.zeros((len(donor_list), len(pep_list)))
    for i in range(len(pep_list)):
        for j in range(len(donor_list)):
            donor_reaction_overview[j,i] = donor_reaction_dict.get(donor_list[j]).get(pep_list[i][1],-1)

    fig, ax = plt.subplots()
    c = plt.imshow(donor_reaction_overview, interpolation='nearest', vmax=25, aspect = "auto")
    ax.set_title('Donor SI per peptide heatmap', fontsize=18)
    ax.set_xlabel("Peptides", fontsize=12)
    ax.set_ylabel("Donors", fontsize=12)
    plt.colorbar(c)
    plt.savefig("../../Figures/Heatmap.png", dpi=500, bbox_inches="tight")
    plt.show(block=False)


def print_corr_plot(chart, non_outliers_list=None, dest = "../../Figures/{}.png"):
    PCC = pearsons_cc(chart[0], chart[1])
    SRC, p = spearmanr(chart[0], chart[1])
    fig, ax = plt.subplots()
    ax.scatter(chart[0], chart[1], label="PCC: {} \nSRC: {}".format(round(PCC,3), round(SRC,3)))
    ax.legend()
    ax.set_xlabel("Ori SI")
    ax.set_ylabel("Var SI")
    ax.set_title("log_{}_vs._{}".format(chart[2][0], chart[2][1]).replace("_"," "))

    file_name = "log_{}_v_{}".format(chart[2][0], chart[2][1])
    file_path = dest.format(file_name)

    fig.savefig(file_path, dpi=500)
    plt.close()

    # fig = plt.figure()
    # ax1 = fig.add_subplot(111)
    # ax1.scatter(chart[0], chart[1], c = "b", label="PCC Before Outlier Selection: %.3f" % PCC)
    # ax1.scatter(non_outliers_list[0],non_outliers_list[1], c ="r", label="PCC After Outlier Selection: %.3f" % PCC_POS)
    # ax1.legend()
    # ax1.set_xlabel("Ori SI")
    # ax1.set_ylabel("Var SI")
    # ax1.set_title(chart[2])
    #
    # plt.show()

def print_stats(bins):
        print("ttest <50% vs. 50%-80%:")
        print(st.ttest_ind(bins[0], bins[1], equal_var=False))
        print("ttest <50% vs. >80%:")
        print(st.ttest_ind(bins[0], bins[2], equal_var=False))
        print("ttest 50%-80% vs. >80%:")
        print(st.ttest_ind(bins[1], bins[2], equal_var=False))


def corr_v_sim_func(cross_react_count, coef_sim_matrix):
    mean_CR = [np.mean(cross_react_count[0]), np.mean(cross_react_count[1]), np.mean(cross_react_count[2])]

    x1, x2, y1, y2 = [],[],[],[]

    for i in range(len(coef_sim_matrix[0])):
        if coef_sim_matrix[0][i] > 0.5:
            y1.append(coef_sim_matrix[0][i])
            x1.append(coef_sim_matrix[2][i])
        else:
            y2.append(coef_sim_matrix[0][i])
            x2.append(coef_sim_matrix[2][i])

    fig, (ax1, ax2) = plt.subplots(2, 1)
    ax1.scatter(x1,y1, c="blue")
    ax1.scatter(x2,y2, c="black")
    ax1.set_ylabel("Pearson corr. coeff.")
    ax2.bar(["<50%", "50%-80%", ">=80%"], mean_CR)
    ax2.set_xlabel("% Sequence identity", labelpad=5)
    ax2.set_ylabel("Fraction significant")
    #plt.savefig("../../Figures/PCC_v_sim.png")
    sim_v_PCC_PCC = pearsons_cc(coef_sim_matrix[2],coef_sim_matrix[0])
    print("PCC for scatterplot.",sim_v_PCC_PCC)
    plt.show()


def pcc_src_comparison(coef_sim_matrix):
    #PCC histogram
    #calculation of optimal number of bins
    q25, q75 = np.percentile(coef_sim_matrix[0],[.25,.75])
    bin_width = 2*(q75 - q25)*len(coef_sim_matrix)**(-1/3)
    bins = round((max(coef_sim_matrix[0]) - min(coef_sim_matrix[0]))/bin_width)

    #Histograms (PCC & SRC)
    fig, (ax1, ax2) = plt.subplots(2, 1)
    ax1.hist(coef_sim_matrix[0], density = False, bins = bins)
    ax1.set_xlabel("PCC value")
    ax1.set_ylabel("Count")
    ax2.hist(coef_sim_matrix[1], density = False, bins = bins)
    ax2.set_xlabel("SRC value")
    ax2.set_ylabel("Count")

    #SCR/PCC scatterplot
    fig, ax1 = plt.subplots()
    ax1.scatter(coef_sim_matrix[0],coef_sim_matrix[1])
    ax1.set_xlabel("PCC")
    ax1.set_ylabel("SRC")
    plt.show()

    #cross-reaction table
    total_peptide_pairs = len(coef_sim_matrix[0])
    n_PCC_CR = sum([pcc > 0.5 for pcc in coef_sim_matrix[0]])
    n_SRC_CR = sum([src > 0.5 for src in coef_sim_matrix[1]])
    print("Peptide pair crossreativity overview:")
    print("Total: {:<8} PCC:{:<12} SRC:{:<12}".format(total_peptide_pairs, n_PCC_CR, n_SRC_CR))


def load_pep_HLA_data(datafile="Data/ragweed/peptides/2860_NetMHCIIpan.xls"):
    infile = open(datafile,"r")

    infile.readline()
    infile.readline()
    pep_HLA_dict = dict()
    old_pep = 0
    for i, line in enumerate(infile):
        line = line.split()

        cur_pep = line[2]
        if old_pep != cur_pep:
            # if old_pep:
            #     pep_HLA_dict[old_pep] = [2 if a<1 else 1 if a<5 else 0 for a in pep_HLA_dict[old_pep]]

            old_pep = cur_pep

        HLA_bind_rank = [[float(line[i]), line[i-2]] for i in range(6,25,3)]

        if cur_pep in pep_HLA_dict:
            pep_HLA_dict[cur_pep] = [old if old[0] < new[0] else new for old, new in zip(pep_HLA_dict[cur_pep], HLA_bind_rank)]
        else:
            pep_HLA_dict[cur_pep] = HLA_bind_rank

    return pep_HLA_dict


def load_peptide_pair_significance(filename):
    infile = open(filename, "r")
    sig_list = [float(line.split()[-1]) for line in infile]
    infile.close()
    return sig_list


def blosum_score(seq1, seq2):
    assert len(seq1) == len(seq2), "blosum_score: sequences are not the same length"
    return sum(blosum50[a][b] for a,b in zip(seq1,seq2))


def outlier_using_IQR(x_values, y_values):
    """ Find outliers using 1.5 x IQR method
    """
    # convert data to np array
    x_values = np.array(chart[0])
    y_values = np.array(chart[1])

    corr_data = np.stack((x_values, y_values))

    #Interquartile range for first array in each pair
    q3_1, q1_1 = np.percentile(corr_data[0], [75, 25])
    IQR1 = q3_1 - q1_1

    #Interquartile range for second array in each pair
    q3_2, q1_2 = np.percentile(corr_data[1], [75, 25])
    IQR2 = q3_2 - q1_2

    non_outliers_list = corr_data[:, np.invert(np.add(corr_data[0, :] > (q3_1 + 1.5 * IQR1), corr_data[1, :] > (q3_2 + 1.5 * IQR2)))]
    non_outliers_list = non_outliers_list[:, np.invert(np.add(non_outliers_list[0, :] < (q1_1 - 1.5 * IQR1), non_outliers_list[1, :] < (q1_2 - 1.5 * IQR2)))]

    return non_outliers_list

#def LOF(array, KNN_n = 2):
    #data = [[a,b] for a,b in zip(array[0], array[1])]
    # for i, j in enumerate(array[0]):
    #     data_point = [array[0][i], array[1][i]]
    #     data.append(data_point)

    # Lav LOF fit
    #clf = LocalOutlierFactor(n_neighbors=KNN_n)
    #clf.fit_predict(data)  # Returns -1 for anomalies/outliers and 1 for inliers.
    #output = clf.negative_outlier_factor_

    #return output

PCC_SRC_switch = 0  # 0 for PCC and 1 for SRC
log_switch = True  # Log transform the correlation data
outlier_sorting = 3    # 0 is nothing, 1 is SRC sig, 2 is PCC sig, 3 is both PCC and SRC sig, 4 is abs(PCC-SRC)
lower_cutoff = 0
bottom_sort_out = 0
HLA_sort = 0

## Main
infile = open("Data/ragweed/ragweed_Tcell_pairwise.MNi.tab", "r")
infile.readline()   # remove header

# Data format:
# ['5540', 'Amb', 'a', '1.0101', 'NSDKTIDGRGAKVEIINAGF', '3.74433',
#          'Amb', 'a', '1.0201', 'NSDKTIDGRGVKVNIVNAGL', '1.12407']

i = -1
old_ori_seq = ""
old_var_seq = ""
charts = []
n = 0
seqs_for_FASTA = []
donor_reaction_dict = dict()
donor_list = []
wanted_charts = 1000

pep_HLA_dict = load_pep_HLA_data()

unique_seqs = set()
allergen_dict = dict()
pep_list = []
pep_dict = dict()
pep_id_name = dict()

for line in infile:
    line = line.split()

    donor_id = line[0]
    if donor_id not in donor_list:
        donor_list.append(donor_id)

    ori_pepseq = line[4]
    var_pepseq = line[9]

    ori_SI = float(line[5])
    var_SI = float(line[10])

    if log_switch:
        ori_SI = np.log(ori_SI)
        var_SI = np.log(var_SI)

    if ori_pepseq != old_ori_seq or var_pepseq != old_var_seq:
        i += 1
        if i > wanted_charts-1:
            break
        n += 1
        old_ori_seq = ori_pepseq
        old_var_seq = var_pepseq

        ori_name = line[1] + "_" + line[2] + "_" + line[3]
        var_name = line[6] + "_" + line[7] + "_" + line[8]

        ori_id = ori_name + "_" + ori_pepseq
        var_id = var_name + "_" + var_pepseq

        if ori_id not in unique_seqs:
            unique_seqs.add(ori_id)

            if ori_name in allergen_dict:
                allergen_dict[ori_name] += 1
            else:
                allergen_dict[ori_name] = 1

            full_ori_name = ori_name + "_" + str(allergen_dict[ori_name])
            pep_list.append([full_ori_name, ori_pepseq, pep_HLA_dict[full_ori_name]])
            pep_dict[full_ori_name] = [ori_pepseq, pep_HLA_dict[full_ori_name]]
            pep_id_name[ori_id] = full_ori_name

        if var_id not in unique_seqs:
            unique_seqs.add(var_id)

            if var_name in allergen_dict:
                allergen_dict[var_name] += 1
            else:
                allergen_dict[var_name] = 1

            full_var_name = var_name + "_" + str(allergen_dict[var_name])
            pep_list.append([full_var_name, var_pepseq, pep_HLA_dict[full_var_name]])
            pep_dict[full_var_name] = [var_pepseq, pep_HLA_dict[full_var_name]]
            pep_id_name[var_id] = full_var_name

        ## Similarity measurements
        # Global alignment with Needleman-Wunsch

        ori_align, var_align, nw_ident, nw_matches = needleman_wunsch(ori_pepseq, var_pepseq)
        ori_align1, var_align1, nw_ident, nw_blosum = needleman_wunsch(ori_pepseq, var_pepseq, scoring_scheme = "blosum")

        # local alignment? with Smith-Waterman (O2)
        scoring_scheme = blosum50
        gap_open = -11
        gap_extension = -1

        P_matrix, Q_matrix, D_matrix, E_matrix, i_max, j_max, max_score = smith_waterman_alignment(ori_pepseq, var_pepseq, scoring_scheme, gap_open, gap_extension)
        aligned_query, aligned_database, sw_matches = smith_waterman_traceback(E_matrix, D_matrix, i_max, j_max, ori_pepseq, var_pepseq, gap_open, gap_extension)

        # print("ALN", full_ori_name, len(ori_pepseq), full_var_name, len(var_pepseq), len(aligned_query), matches, max_score)
        # print("QAL", i_max, ''.join(aligned_query))
        # print("DAL", j_max,''.join(aligned_database))
        # print("")

        # Save info about the peptide pairing and prep lists for SI values
        charts.append([[],[], [pep_id_name[ori_id], pep_id_name[var_id]], [nw_ident,sw_matches,max_score,nw_blosum], 0])

    # add SI to the chart for the peptide pairing
    charts[i][0].append(ori_SI)
    charts[i][1].append(var_SI)
    charts[i][4] += 1       # n donor

    if donor_id not in donor_reaction_dict:
        donor_reaction_dict[donor_id] = dict()

    donor_reaction_dict[donor_id][ori_pepseq] = ori_SI
    donor_reaction_dict[donor_id][var_pepseq] = var_SI

infile.close()
no_bind = 0
weak_bind = 0
strong_bind = 0
for pep in pep_list:
    bind_str = min([a[0] for a in pep[2]])
    if bind_str > 5:
        pep_dict[pep[0]].append(0)
        no_bind += 1
    elif bind_str >= 1:
        pep_dict[pep[0]].append(1)
        weak_bind += 1
    else:
        pep_dict[pep[0]].append(2)
        strong_bind += 1


# print("No binders:", no_bind, "Weak binders:", weak_bind, "Strong binders:", strong_bind)

# # Print seqs for fasta format
# outfile = open("seqs_for_HLA_profiling.fsa", "w")
# for [name, seq] in pep_list:
#     print(">" + name, file=outfile)
#     print(seq, file=outfile)
# outfile.close()

# heatmap(pep_list, donor_list, donor_reaction_dict)

coef_sim_matrix = [[],[],[],[],[],[],[],[],[],[],[],[],[]]
cross_react_count = [[],[],[]]
PCC_bins = [[],[],[]]
sensitive_plots = []
HLA_binder_table = [np.zeros((2,3)),np.zeros((2,3)),np.zeros((2,3)),np.zeros((2,3)),np.zeros((2,3)),np.zeros((2,3))]
HLA_binder_table_2 = [np.zeros((2,3)),np.zeros((2,3)),np.zeros((2,3)),np.zeros((2,3)),np.zeros((2,3)),np.zeros((2,3))]
table_count = [0,0,0,0,0,0]
CR_delta_rank = []
NCR_delta_rank = []
threshold_values = [0.17 if log_switch else 0.3]
PCC_t = []

# load the significances for each chart
if log_switch:
    SRC_sig_list = load_peptide_pair_significance("Data/ragweed/significances/ragweed_log_PCC_sig.txt")
    PCC_sig_list = load_peptide_pair_significance("Data/ragweed/significances/ragweed_log_SCC_sig.txt")
else:
    SRC_sig_list = load_peptide_pair_significance("Data/ragweed/significances/ragweed_SCC_sig.txt")
    PCC_sig_list = load_peptide_pair_significance("Data/ragweed/significances/ragweed_PCC_sig.txt")

for t in threshold_values:

    for i, (chart, SRC_sig, PCC_sig) in enumerate(zip(charts, SRC_sig_list, PCC_sig_list)):
        # outfile = open("charts/log_chart{}.txt".format(i), "w")
        # for x,y in zip(chart[0], chart[1]):
        #     print(x,y,sep="\t",file=outfile)
        # outfile.close()

        PCC = pearsons_cc(chart[0], chart[1])
        SRC, p = spearmanr(chart[0], chart[1])
        global_ident = chart[3][0]
        global_blosum = chart[3][3]
        local_matches = chart[3][1]
        local_blosum = chart[3][2]

        ## Outlier chart sorting
        # by SRC significance
        sorted_out = 0
        if outlier_sorting == 1 or outlier_sorting == 3:
            if SRC_sig > 0.05 and SRC > 0.5:
                sorted_out = 1

            elif SRC_sig < 0.05 and SRC < -0.25:
                sorted_out = bottom_sort_out
            elif SRC < -0.25 and lower_cutoff:
                SRC = -0.25
        # by PCC significance
        if outlier_sorting == 2 or outlier_sorting == 3:
            if PCC_sig > 0.05 and PCC > 0.5:
                sorted_out = 1

            elif PCC_sig < 0.05 and PCC < -0.25:
                sorted_out = bottom_sort_out
            elif PCC < -0.25 and lower_cutoff:
                PCC = -0.25

        # difference between PCC and SRC
        if outlier_sorting == 4:
            if np.abs(PCC-SRC) > t:
                sorted_out = 1


        #point = LOF(corr_data)
        #outliers.append(point)
        #print(outliers)
        #print(outliers)
        #print(corr_data)

        #for i, dp in enumerate(point):
            #if dp <= -10:
            #    point[i] = int(1)
            #else:
            #    point[i] = int(0)

        #point = point.astype(bool) #returns binary array for every dataset (1 = outlier, 0 = inlier)
        #outlier_list = corr_data[:, point]
        #non_outliers_list = corr_data[:, np.invert(point)]
        #print(corr_data)
        #print(non_outliers_list)
        # PCC_POS = pearsons_cc(*non_outliers_list) #POS = "post outlier selection"
        # Delta_PCC = abs(PCC_POS - PCC)


        # print_corr_plot(chart, non_outliers_list)

        #if Delta_PCC >= 0.3:
            #sensitive_plots.append(Delta_PCC)
        #print(sensitive_plots)

        # print("{:<8} {:<12} {:<12} {:<10}".format("n = %.d" % chart[4], "PCC: %.3f" % PCC, "SRC: %.3f" % SRC, "N_sim: %.d " % chart[3]))

        pep1_info = pep_dict[chart[2][0]]
        pep2_info = pep_dict[chart[2][1]]

        pep1_bind = pep1_info[1]
        pep2_bind = pep2_info[1]

        ori_bind = bool(sum([rank <= 5 for [rank, core] in pep_HLA_dict[chart[2][0]]]))
        var_bind = bool(sum([rank <= 5 for [rank, core] in pep_HLA_dict[chart[2][1]]]))

        if HLA_sort and not (ori_bind and var_bind):
            sorted_out = 1

        ## Core similarity vs PCC
        # best core vs best core
        rank1 = np.inf
        rank2 = np.inf
        for [rank, core] in pep1_bind:
            if rank < rank1:
                rank1 = rank
                best_core1 = core

        for [rank, core] in pep2_bind:
            if rank < rank2:
                rank2 = rank
                best_core2 = core

        # best core from ori and corresponding core for var
        best_rank = np.inf
        for [rank_ori, core_ori], [rank_var, core_var]  in zip(pep1_bind, pep2_bind):
            if rank_ori < best_rank:
                best_rank = rank_ori
                var_rank = rank_var
                best_core3 = core_ori
                best_core4 = core_var

        core_matches = 0
        core_blosum = 0
        for a,b in zip(best_core1, best_core2):
            core_matches += int(a==b)
            core_blosum += blosum50[a][b]

            core_ident = core_matches/len(best_core1)*100

        best_blosum_cores = "NA"
        best_ident_cores = "NA"
        best_ident = None
        best_blosum = None

        for [rank_ori, core_ori] in pep1_bind:
            if rank_ori > 5:
                continue
            for [rank_var, core_var] in pep2_bind:
                if rank_var > 5:
                    continue

                blosum_match = blosum_score(core_ori, core_var)
                ident = sum([int(a == b) for a,b in zip(core_ori, core_var)])/len(core_ori)*100

                if best_blosum == None:
                    best_blosum = blosum_match
                elif blosum_match > best_blosum:
                    best_blosum = blosum_match
                    best_blosum_cores = [core_ori, core_var]

                if best_ident == None:
                    best_ident = ident
                elif ident > best_ident:
                    best_ident = ident
                    best_ident_cores = [core_ori, core_var]


        delta_rank = np.abs(rank1 - rank2)

        glob_sim_and_d_rank = global_ident*(100-delta_rank)

        coef_sim_matrix[0].append(PCC)
        coef_sim_matrix[1].append(SRC)
        coef_sim_matrix[2].append(global_ident)
        coef_sim_matrix[3].append(core_ident)
        coef_sim_matrix[4].append(core_blosum)
        coef_sim_matrix[5].append(delta_rank)
        coef_sim_matrix[6].append(best_ident)
        coef_sim_matrix[7].append(best_blosum)
        coef_sim_matrix[8].append(glob_sim_and_d_rank)
        coef_sim_matrix[9].append(local_matches)
        coef_sim_matrix[10].append(local_blosum)
        coef_sim_matrix[11].append(sorted_out)
        coef_sim_matrix[12].append(global_blosum)


        str_bind = int(pep_dict[chart[2][0]][2] == 2) + int(pep_dict[chart[2][1]][2] == 2)
        weak_bind = int(pep_dict[chart[2][0]][2] > 0) + int(pep_dict[chart[2][1]][2] > 0)

        str_bind_2 = min(sum([int(a[0] < 1 and b[0] < 1) for a,b in zip(pep_dict[chart[2][0]][1],pep_dict[chart[2][1]][1])]),2)
        weak_bind_2 = min(sum([int(a[0] < 5 and b[0] < 5) for a,b in zip(pep_dict[chart[2][0]][1],pep_dict[chart[2][1]][1])]),2)

        CR = 1 if SRC > 0.5 else 0

        # for delta rank boxplot
        if CR:
            CR_delta_rank.append(delta_rank)
        else:
            NCR_delta_rank.append(delta_rank)

        # Print HLA binding tables and PCC for <50%, 50-80% and >80% bins.
        if global_ident < 50:
            cross_react_count[0].append(CR)
            PCC_bins[0].append(PCC)
            if not CR:
                HLA_binder_table[0][0, str_bind] += 1
                HLA_binder_table[0][1, weak_bind] += 1
                HLA_binder_table_2[0][0, str_bind_2] += 1
                HLA_binder_table_2[0][1, weak_bind_2] += 1
                table_count[0] += 1
            else:
                HLA_binder_table[3][0, str_bind] += 1
                HLA_binder_table[3][1, weak_bind] += 1
                HLA_binder_table_2[3][0, str_bind_2] += 1
                HLA_binder_table_2[3][1, weak_bind_2] += 1
                table_count[3] += 1
        elif global_ident >= 80:
            cross_react_count[2].append(CR)
            PCC_bins[2].append(PCC)
            if not CR:
                HLA_binder_table[2][0, str_bind] += 1
                HLA_binder_table[2][1, weak_bind] += 1
                HLA_binder_table_2[2][0, str_bind_2] += 1
                HLA_binder_table_2[2][1, weak_bind_2] += 1
                table_count[2] += 1
            else:
                HLA_binder_table[5][0, str_bind] += 1
                HLA_binder_table[5][1, weak_bind] += 1
                HLA_binder_table_2[5][0, str_bind_2] += 1
                HLA_binder_table_2[5][1, weak_bind_2] += 1
                table_count[5] += 1
        else:
            cross_react_count[1].append(CR)
            PCC_bins[1].append(PCC)
            if not CR:
                HLA_binder_table[1][0, str_bind] += 1
                HLA_binder_table[1][1, weak_bind] += 1
                HLA_binder_table_2[1][0, str_bind_2] += 1
                HLA_binder_table_2[1][1, weak_bind_2] += 1
                table_count[1] += 1
            else:
                HLA_binder_table[4][0, str_bind] += 1
                HLA_binder_table[4][1, weak_bind] += 1
                HLA_binder_table_2[4][0, str_bind_2] += 1
                HLA_binder_table_2[4][1, weak_bind_2] += 1
                table_count[4] += 1
#     PCC_0 = pearsons_cc(coef_sim_matrix[2], coef_sim_matrix[1])
#     PCC_t.append(PCC_0)
#
# PCC_t = np.array(PCC_t)
# print(threshold_values[PCC_t == max(PCC_t)])

#label = ["Low sim, low SRC","Mid sim, low SRC","High sim, low SRC","Low sim, high SRC","Mid sim, high SRC","high sim, high SRC"]
#print("Pool size")
#print(table_count[:3])
#print(table_count[3:])
#for i, (n, table) in enumerate(zip(table_count,HLA_binder_table_2)):
    #print(table)
    #print()
    #print(label[i])
    #print(np.round(table/n,2))

PCC_SRC_str = "PCC " if PCC_SRC_switch == 0 else "SRC "
log_norm_str = "log transformed " if log_switch else "normal "
sorts = ["none", "SRC sig", "PCC sig", "SRC and PCC sig", "abs(PCC-SRC)<t"]
sorting_str = ". Sorted: " + sorts[outlier_sorting]
overall_title = PCC_SRC_str + log_norm_str + "data" + sorting_str

sorted_out_list = np.array(coef_sim_matrix[11]).astype(bool)

fig, ax1 = plt.subplots()

# Global similarity vs PCC/SRC
array_1 = np.array((coef_sim_matrix[2], coef_sim_matrix[PCC_SRC_switch]))
array_1_in = array_1[:, np.invert(sorted_out_list)]
array_1_out = array_1[:, sorted_out_list]

print(len(array_1_in[0]))
print(len(array_1_out[0]))

PCC_1_before = pearsons_cc(*array_1)
PCC_1_after = pearsons_cc(*array_1_in)
ax1.scatter(*array_1_in)
ax1.scatter(*array_1_out)
ax1.set_title(PCC_SRC_str + "as a function of global similarity(percent). PCC = %.3f. After = %.3f" % (PCC_1_before, PCC_1_after), fontsize=14)
ax1.set_xlabel("Global similarity(percent)")
ax1.set_ylabel(PCC_SRC_str)
plt.show()

#
# # best v best core identity vs PCC/SRC
# array_2 = np.array((coef_sim_matrix[3], coef_sim_matrix[PCC_SRC_switch]))
# array_2_in = array_2[:, np.invert(sorted_out_list)]
# array_2_out = array_2[:, sorted_out_list]
#
# PCC_2_before = pearsons_cc(*array_2)
# PCC_2_after = pearsons_cc(*array_2_in)
# ax2.scatter(*array_2_in)
# ax2.scatter(*array_2_out)
# ax2.set_title(PCC_SRC_str + "as a function of core identity(percent). PCC = %.3f. After = %.3f" % (PCC_2_before, PCC_2_after), fontsize=14)
# ax2.set_xlabel("Core identity(%)")
# ax2.set_ylabel(PCC_SRC_str)
# # plt.show()
#
# # best v best core blosum vs PCC/SRC
# array_3 = np.array((coef_sim_matrix[4], coef_sim_matrix[PCC_SRC_switch]))
# array_3_in = array_3[:, np.invert(sorted_out_list)]
# array_3_out = array_3[:, sorted_out_list]
#
# PCC_3_before = pearsons_cc(*array_3)
# PCC_3_after = pearsons_cc(*array_3_in)
# ax3.scatter(*array_3_in)
# ax3.scatter(*array_3_out)
# ax3.set_title(PCC_SRC_str + "as a function of core identity(BLOSUM score). PCC = %.3f. After = %.3f" % (PCC_3_before, PCC_3_after), fontsize=14)
# ax3.set_xlabel("Core identity(BLOSUM score)")
# ax3.set_ylabel(PCC_SRC_str)
# # plt.show()
#
#
# # global sim x delta_rank vs PCC/SRC
# array_4 = np.array((coef_sim_matrix[8], coef_sim_matrix[PCC_SRC_switch]))
# array_4_in = array_4[:, np.invert(sorted_out_list)]
# array_4_out = array_4[:, sorted_out_list]
#
# PCC_4_before = pearsons_cc(*array_4)
# PCC_4_after = pearsons_cc(*array_4_in)
# ax4.scatter(*array_4_in)
# ax4.scatter(*array_4_out)
# ax4.set_title(PCC_SRC_str + "as a function of global similarity and delta_rank. PCC = %.3f. After = %.3f" % (PCC_4_before, PCC_4_after), fontsize=14)
# ax4.set_xlabel("Global similarity(percent)*(100-delta_rank)")
# ax4.set_ylabel(PCC_SRC_str)
#
# # best matching binding cores identity vs PCC/SRC
# array_5 = np.array((coef_sim_matrix[6], coef_sim_matrix[PCC_SRC_switch]))
# sorted_out_list1 = sorted_out_list[array_5[0,:] != None]
# array_5 = array_5[:,array_5[0,:] != None]
# array_5_in = array_5[:, np.invert(sorted_out_list1)]
# array_5_out = array_5[:, sorted_out_list1]
#
# PCC_5_before = pearsons_cc(*array_5)
# PCC_5_after = pearsons_cc(*array_5_in)
# ax5.scatter(*array_5_in)
# ax5.scatter(*array_5_out)
# ax5.set_title(PCC_SRC_str + "as a function of best matching cores identity(percent). PCC = %.3f. After = %.3f" % (PCC_5_before, PCC_5_after), fontsize=14)
# ax5.set_xlabel("Core identity(%)")
# ax5.set_ylabel(PCC_SRC_str)
# # plt.show()
#
# # best matching binding cores blosum vs PCC/SRC
# array_6 = np.array((coef_sim_matrix[7], coef_sim_matrix[PCC_SRC_switch]))
# sorted_out_list2 = sorted_out_list[array_6[0,:] != None]
# array_6 = array_6[:,array_6[0,:] != None]
# array_6_in = array_6[:, np.invert(sorted_out_list2)]
# array_6_out = array_6[:, sorted_out_list2]
#
# PCC_6_before = pearsons_cc(*array_6)
# PCC_6_after = pearsons_cc(*array_6_in)
# ax6.scatter(*array_6_in)
# ax6.scatter(*array_6_out)
# ax6.set_title(PCC_SRC_str + "as a function of best matching cores BLOSUM score. PCC = %.3f. After = %.3f" % (PCC_6_before, PCC_6_after), fontsize=14)
# ax6.set_xlabel("Core identity(BLOSUM score)")
# ax6.set_ylabel(PCC_SRC_str)
#
# plt.show()
#
# fig, ((ax9,ax7),(ax10,ax8)) = plt.subplots(2,2)
# # Global similarity vs PCC/SRC
# array_9 = np.array((coef_sim_matrix[2], coef_sim_matrix[PCC_SRC_switch]))
# array_9_in = array_9[:, np.invert(sorted_out_list)]
# array_9_out = array_9[:, sorted_out_list]
#
# PCC_9_before = pearsons_cc(*array_9)
# PCC_9_after = pearsons_cc(*array_9_in)
# ax9.scatter(*array_9_in)
# ax9.scatter(*array_9_out)
# ax9.set_title(PCC_SRC_str + "as a function of global similarity(percent). PCC = %.3f. After = %.3f" % (PCC_9_before, PCC_9_after), fontsize=14)
# ax9.set_xlabel("Global similarity(percent)")
# ax9.set_ylabel(PCC_SRC_str)
# # plt.show()
# # Global similarity vs PCC/SRC
# array_10 = np.array((coef_sim_matrix[12], coef_sim_matrix[PCC_SRC_switch]))
# array_10_in = array_10[:, np.invert(sorted_out_list)]
# array_10_out = array_10[:, sorted_out_list]
#
# PCC_10_before = pearsons_cc(*array_10)
# PCC_10_after = pearsons_cc(*array_10_in)
# ax10.scatter(*array_10_in)
# ax10.scatter(*array_10_out)
# ax10.set_title(PCC_SRC_str + "as a function of global similarity(percent). PCC = %.3f. After = %.3f" % (PCC_10_before, PCC_10_after), fontsize=14)
# ax10.set_xlabel("Global similarity(percent)")
# ax10.set_ylabel(PCC_SRC_str)
# # plt.show()
# # best matching binding cores identity vs PCC/SRC
# array_7 = np.array((coef_sim_matrix[9], coef_sim_matrix[PCC_SRC_switch]))
# array_7 = array_7[:,array_7[0,:] != None]
# array_7_in = array_7[:, np.invert(sorted_out_list)]
# array_7_out = array_7[:, sorted_out_list]
#
# PCC_7_before = pearsons_cc(*array_7)
# PCC_7_after = pearsons_cc(*array_7_in)
# ax7.scatter(*array_7_in)
# ax7.scatter(*array_7_out)
# ax7.set_title(PCC_SRC_str + "as a function of local alignment matches. PCC = %.3f. After = %.3f" % (PCC_7_before, PCC_7_after), fontsize=14)
# ax7.set_xlabel("Number of matches")
# ax7.set_ylabel(PCC_SRC_str)
# # plt.show()
#
# # best matching binding cores blosum vs PCC/SRC
# array_8 = np.array((coef_sim_matrix[10], coef_sim_matrix[PCC_SRC_switch]))
# array_8 = array_8[:,array_8[0,:] != None]
# array_8_in = array_8[:, np.invert(sorted_out_list)]
# array_8_out = array_8[:, sorted_out_list]
#
# PCC_8_before = pearsons_cc(*array_8)
# PCC_8_after = pearsons_cc(*array_8_in)
# ax8.scatter(*array_8_in)
# ax8.scatter(*array_8_out)
# ax8.set_title(PCC_SRC_str + "as a function of local alignment BLOSUM score. PCC = %.3f. After = %.3f" % (PCC_8_before, PCC_8_after), fontsize=14)
# ax8.set_xlabel("BLOSUM score)")
# ax8.set_ylabel(PCC_SRC_str)
#
# plt.show()
# fig, ax7 = plt.subplots(1,1)
# p_val = st.ttest_ind(NCR_delta_rank,CR_delta_rank, equal_var=False)[1]
# ax7.boxplot([NCR_delta_rank,CR_delta_rank], vert = 0)
# ax7.set_yticklabels(["Non-CR", "CR"])
# ax7.set_xlabel("Delta rank")
# ax7.set_title("Delta rank for CR and non CR. p-val = %.10f" % p_val)
# plt.show()

# fig, ax8 = plt.subplots(1,1)
# ax8.scatter(threshold_values, PCC_t)
# plt.show()
# print(max(PCC_t))

# print("Crossreaction frequency t-test")
# print_stats(cross_react_count)
# ---
# print("PCC t-test")
# print_stats(PCC_bins)

# print("Correlation vs. similarity + bar plot")
# corr_v_sim_func(cross_react_count, coef_sim_matrix)
#
# print("PCC and SRC histogram")
# pcc_src_comparison(coef_sim_matrix)
