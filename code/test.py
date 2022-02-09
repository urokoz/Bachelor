#!/usr/bin/env python

infile = open("HLA_ragweed.csv", "r")

header = infile.readline().strip().split(";")

HLA_set = set()

donor_HLA_dict = dict()
HLA_count_dict = dict()
for line in infile:
    line = line.strip().split(";")

    for i, type in enumerate(line):
        if i == 0:
            donor = type
            donor_HLA_dict[donor] = []
            continue

        HLA = header[i]

        word = ""

        if type == "-" or type == "NA":
            continue

        for char in type:
            if char == "*":
                HLA += "*"

            elif char == ":":
                HLA_type = HLA + word + ":"
                word = ""

            elif char == "/":
                if "-" in word:
                    start, finish = word.split("-")
                    for i in range(int(start),int(finish)+1):
                        final_HLA = HLA_type + str(i)
                        HLA_set.add(final_HLA)
                        donor_HLA_dict[donor].append(final_HLA)
                        if HLA_count_dict.get(final_HLA):
                            HLA_count_dict[final_HLA] += 1
                        else:
                            HLA_count_dict[final_HLA] = 1

                else:
                    final_HLA = HLA_type + word
                    HLA_set.add(final_HLA)
                    donor_HLA_dict[donor].append(final_HLA)
                word = ""
            else:
                word += char

HLA_list = list(HLA_set)

HLA_list = sorted(HLA_list)
i = 0
for allele, count in HLA_count_dict.items():
    if count > 1:
        print(allele, count, sep="\t")
        i += 1

print(i)
#
# for donor, alleles in donor_HLA_dict.items():
#     print(donor, alleles, sep="\t")
