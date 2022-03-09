#!/usr/bin/env python

import sys
import numpy as np
import math
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from scipy.stats import spearmanr, pearsonr

def heatmap(pep_list, donor_list, donor_reaction_dict):
    # Heatmap generation
    print("n peptides:", len(pep_list))
    print("n donors:", len(donor_list))

    donor_reaction_overview = np.zeros((len(donor_list), len(pep_list)))
    for i in range(len(pep_list)):
        for j in range(len(donor_list)):
            donor_reaction_overview[j,i] = donor_reaction_dict.get(donor_list[j]).get(pep_list[i][1],np.nan)

    fig, ax = plt.subplots()
    current_cmap = plt.cm.get_cmap()
    current_cmap.set_bad(color='grey')
    current_cmap.set_over(color='red')
    c = plt.imshow(donor_reaction_overview,vmax=25, interpolation='nearest', aspect = "auto")
    ax.set_xlabel("Peptides", fontsize=12)
    ax.set_ylabel("Donors", fontsize=12)
    plt.colorbar(c)
    plt.savefig("../../Figures/Heatmap.png", dpi=500, bbox_inches="tight")
    plt.show()


def load_pep_HLA_data(datafile="Data/2860_NetMHCIIpan.xls"):
    """ Reads HLA binding for 7 allels and outputs it as a dict.

        output format:
        pep_HLA_dict[pep_name] = [[rank1,core1],...,[rank7,core7]]
    """

    infile = open(datafile,"r")

    infile.readline()
    infile.readline()
    pep_HLA_dict = dict()
    old_pep = 0
    for i, line in enumerate(infile):
        line = line.split()

        cur_pep = line[2]
        if old_pep != cur_pep:
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
    PCC = pearsons_cc(x,y)
    #SCC, _ = st.spearmanr(x,y)
    fig, ax = plt.subplots()
    ax.scatter(x, y)
    ax.set_title(plot_title + " PCC: {}".format(round(PCC,3)))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.show()




## Argument parsing:

parser = ArgumentParser(description="Preps and filters the data and peptides")
parser.add_argument("-f", action="store", dest="data_file", type=str, default = "ragweed_lookup_file.txt", help="Lookup file with paths of datafiles")
parser.add_argument("-hs", action="store_true", default=False, help="Sort out nonbinders")
parser.add_argument("-log", action="store_true", default=False, help="Log-transform SI values")
parser.add_argument("-lc", action="store_true", default=False, help="Raise PCC/SCC values under -0.25 to -0.25")
parser.add_argument("-bs", action="store_true", default=False, help="Sort out significant datapoints under -0.25 PCC/SCC")
parser.add_argument("-ch", action="store_true", default=False, help="Print charts for permutation test")
parser.add_argument("-ol", action="store", type=int, default=0, help="Significance sorting: 0: none, 1: SRC sig, 2: PCC sig, 3: both PCC and SRC sig")
parser.add_argument("-fas_f", action="store", type=str, default=None, help="Print fasta of peptides")
parser.add_argument("-hm", action="store_true", default=False, help="Changes output to heatmap only")


args = parser.parse_args()
lookup_file = args.data_file
HLA_sort = args.hs
log_switch = args.log
lower_cutoff = args.lc
bottom_sort_out = args.bs
print_charts = args.ch
outlier_sorting = args.ol
fasta_name = args.fas_f
heatmap_switch = args.hm

infile = open(lookup_file, "r")
lookup_dict = dict()
for line in infile:
    line = line.split()
    lookup_dict[line[0]] = line[1]
infile.close()

data_dir = lookup_dict["data_dir"]
data_file = data_dir + lookup_dict["dataset"]
pep_file = data_dir + lookup_dict["pep_file"]
hla_file = data_dir + lookup_dict["HLA_file"]

## Main

# Data format:
# ['5540', 'Amb', 'a', '1.0101', 'NSDKTIDGRGAKVEIINAGF', '3.74433',
#          'Amb', 'a', '1.0201', 'NSDKTIDGRGVKVNIVNAGL', '1.12407']

seqs_for_FASTA = []
donor_list = []
pep_list = []
Ori_SI_means = []

pep_pair_dict = dict()
seen_peptide_pairs = set()
donor_reaction_dict = dict()
unique_seqs = set()
allergen_dict = dict()
pep_dict = dict()
pep_id_name = dict()


infile = open(data_file, "r")
infile.readline()   # remove header
for line in infile:
    line = line.replace("\n","").split("\t")

    donor_id = line[0]
    if donor_id not in donor_list:
        donor_list.append(donor_id)

    ori_pepseq = line[2]
    var_pepseq = line[5]

    ori_SI = float(line[3])
    var_SI = float(line[6])

    ori_name = line[1].strip().replace(" ", "_")
    var_name = line[4].strip().replace(" ", "_")

    ori_id = ori_name + "_" + ori_pepseq
    var_id = var_name + "_" + var_pepseq

    if log_switch:
        ori_SI = np.log(ori_SI)
        var_SI = np.log(var_SI)

    if frozenset([ori_id,var_id]) not in seen_peptide_pairs:
        seen_peptide_pairs.add(frozenset([ori_id,var_id]))

        if ori_id not in unique_seqs:
            unique_seqs.add(ori_id)

            if ori_name in allergen_dict:
                allergen_dict[ori_name] += 1
            else:
                allergen_dict[ori_name] = 1

            full_ori_name = ori_name + "_" + str(allergen_dict[ori_name])
            pep_list.append([full_ori_name, ori_pepseq])
            pep_id_name[ori_id] = full_ori_name

        if var_id not in unique_seqs:
            unique_seqs.add(var_id)

            if var_name in allergen_dict:
                allergen_dict[var_name] += 1
            else:
                allergen_dict[var_name] = 1

            full_var_name = var_name + "_" + str(allergen_dict[var_name])
            pep_list.append([full_var_name, var_pepseq])
            pep_id_name[var_id] = full_var_name

        # Save info about the peptide pairing and prep lists for SI values
        full_ori_name = pep_id_name[ori_id]
        full_var_name = pep_id_name[var_id]
        if full_ori_name < full_var_name:
            pep_pair = full_ori_name + "\t" + full_var_name
        else:
            pep_pair = full_var_name + "\t" + full_ori_name
        pep_pair_dict[pep_pair] = [[],[]]

    # add SI to the chart for the peptide pairing
    full_ori_name = pep_id_name[ori_id]
    full_var_name = pep_id_name[var_id]
    if pep_id_name[ori_id] < pep_id_name[var_id]:
        pep_pair = full_ori_name + "\t" + full_var_name
        pep_pair_dict[pep_pair][0].append(ori_SI)
        pep_pair_dict[pep_pair][1].append(var_SI)
    else:
        pep_pair = full_var_name + "\t" + full_ori_name
        pep_pair_dict[pep_pair][0].append(var_SI)
        pep_pair_dict[pep_pair][1].append(ori_SI)



    if donor_id not in donor_reaction_dict:
        donor_reaction_dict[donor_id] = dict()

    donor_reaction_dict[donor_id][ori_pepseq] = ori_SI
    donor_reaction_dict[donor_id][var_pepseq] = var_SI
infile.close()

if log_switch:
    SRC_sig_list = load_peptide_pair_significance(data_dir + lookup_dict["log_SCC_sig"])
    PCC_sig_list = load_peptide_pair_significance(data_dir + lookup_dict["log_PCC_sig"])
else:
    SRC_sig_list = load_peptide_pair_significance(data_dir + lookup_dict["SCC_sig"])
    PCC_sig_list = load_peptide_pair_significance(data_dir + lookup_dict["PCC_sig"])

pep_HLA_dict = load_pep_HLA_data(hla_file)

sig_counter = 0
HLA_counter = 0

if heatmap_switch:
    heatmap(pep_list, donor_list, donor_reaction_dict)
else:
    outfile_name = "prepped_data/"
    if log_switch:
        outfile_name += "log_"
    if HLA_sort:
        outfile_name += "filtered_"
    else:
        outfile_name += "unfiltered_"

    outfile = open(data_dir+outfile_name+"dataset.csv","w")
    for i, ((pep_pair, SI_vals), SCC_sig, PCC_sig) in enumerate(zip(pep_pair_dict.items(), SRC_sig_list, PCC_sig_list)):
        ori_SI = SI_vals[0]
        var_SI = SI_vals[1]
        ori_name = pep_pair.split()[0]
        var_name = pep_pair.split()[1]

        # Print SI values as chart files for permutation test
        if print_charts:
            chartname = "chart{}.txt".format(i)
            if log_switch:
                chartname = "log_" + chartname

            chartfile = open(data_dir + "charts/" + chartname,"w")

            #prints scatterplots args: (x,y,plot_title, xlabel, ylabel)
            print(sim_scatterplot(ori_SI, var_SI, ori_name + "vs." + var_name, ori_name, var_name))


            for x,y in zip(ori_SI, var_SI):
                print(x,y,sep="\t",file=chartfile)
            chartfile.close()


        ori_bind = bool(sum([rank <= 5 for [rank, core] in pep_HLA_dict[ori_name]]))
        var_bind = bool(sum([rank <= 5 for [rank, core] in pep_HLA_dict[var_name]]))

        if HLA_sort and not (ori_bind and var_bind):
            HLA_counter += 1
            continue

        PCC, PCC_p = pearsonr(ori_SI,var_SI)
        SCC, SCC_p = spearmanr(ori_SI,var_SI)

        # by SRC significance
        if outlier_sorting == 1 or outlier_sorting == 3:
            if SCC_sig > 0.05 and SCC > 0.5:
                sig_counter += 1
                continue

            elif SCC_sig < 0.05 and SCC < -0.25 and bottom_sort_out:
                sig_counter += 1
                continue
            elif SCC < -0.25 and lower_cutoff:
                SCC = -0.25
        # by PCC significance
        if outlier_sorting == 2 or outlier_sorting == 3:
            if PCC_sig > 0.05 and PCC > 0.5:
                sig_counter += 1
                continue

            elif PCC_sig < 0.05 and PCC < -0.25 and bottom_sort_out:
                sig_counter += 1
                continue
            elif PCC < -0.25 and lower_cutoff:
                PCC = -0.25


        print(ori_name, var_name, ",".join(str(e) for e in ori_SI), ",".join(str(e) for e in var_SI), sep="\t", file=outfile)
    outfile.close()

    outfile = open(pep_file, "w")

    if fasta_name != None:
        fasta_file = open(fasta_name, "w")

    for pep in pep_list:
        pep_name = pep[0]
        pep_seq = pep[1]
        pep_HLA = pep_HLA_dict[pep_name]
        pep_bind = bool(sum([rank <= 5 for [rank, core] in pep_HLA]))

        if fasta_name != None:
            print(">" + pep_name, file=fasta_file)
            print(pep_seq, file=fasta_file)


        print(pep_name, pep_seq, " ".join(",".join(str(f) for f in e) for e in pep_HLA), sep="\t", file = outfile)

print("Sig: ", sig_counter)
print("HLA: ", HLA_counter)
