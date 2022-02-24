#!/usr/bin/env python

infile = open("HLA_tree.csv", "r")

header = infile.readline().strip().split(";")

donor_HLA_dict = dict()
HLA_count_dict = dict()
for line in infile:
    line = line.strip().split(";")

    for i, type in enumerate(line):
        if i == 0:
            donor = type
            donor_HLA_dict[donor] = []
            continue

        if bool(int(type)):
            donor_HLA_dict[donor].append(header[i])
            if HLA_count_dict.get(header[i]):
                HLA_count_dict[header[i]] += 1
            else:
                HLA_count_dict[header[i]] = 1

outfile = open("tree_HLA_count.txt", "w")

for allele, count in sorted(HLA_count_dict.items()):
    print(allele, count, sep="\t", file = outfile)

outfile.close()

outfile = open("tree_donor_HLA.txt", "w")

for donor, alleles in donor_HLA_dict.items():
    print(donor, ",".join(sorted(alleles)), sep="\t", file = outfile)

outfile.close()
