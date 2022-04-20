#!/usr/bin/env python

infile = open("HLA_ragweed.csv", "r")

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

        HLA = header[i].split("-")[1]

        if HLA[:3] != "DRB":
            continue

        word = ""

        if type == "-" or type == "NA" or type == "Not Tested":
            continue

        for char in type:
            if char == "*":
                HLA += "_"

            elif char == ":":
                HLA_type = HLA + word #+ ":"
                word = ""

            elif char == "/":
                if "-" in word:
                    start, finish = word.split("-")

                    final_HLA = HLA_type + start

                    donor_HLA_dict[donor].append(final_HLA)

                    if HLA_count_dict.get(final_HLA):
                        HLA_count_dict[final_HLA] += 1
                    else:
                        HLA_count_dict[final_HLA] = 1

                else:
                    final_HLA = HLA_type + word

                    donor_HLA_dict[donor].append(final_HLA)
                    if HLA_count_dict.get(final_HLA):
                        HLA_count_dict[final_HLA] += 1
                    else:
                        HLA_count_dict[final_HLA] = 1
                word = ""
                break

            else:
                word += char
        if not (type == "-" or type == "NA" or type == "Not Tested" or word == ""):
            if "-" in word:
                start, finish = word.split("-")

                final_HLA = HLA_type + start

                donor_HLA_dict[donor].append(final_HLA)

                if HLA_count_dict.get(final_HLA):
                    HLA_count_dict[final_HLA] += 1
                else:
                    HLA_count_dict[final_HLA] = 1

            else:
                final_HLA = HLA_type + word

                donor_HLA_dict[donor].append(final_HLA)
                if HLA_count_dict.get(final_HLA):
                    HLA_count_dict[final_HLA] += 1
                else:
                    HLA_count_dict[final_HLA] = 1

outfile = open("ragweed_HLA_count.txt", "w")

for allele, count in sorted(HLA_count_dict.items()):
    print(allele, count, sep="\t", file = outfile)

outfile.close()

outfile = open("ragweed_donor_HLA.txt", "w")

for donor, alleles in donor_HLA_dict.items():
    print(donor, ",".join(sorted(alleles)), sep="\t", file = outfile)

outfile.close()
