#!/usr/bin/env python
import sys


infile = open(sys.argv[1], "r")

for line in infile:
    donor, alleles = line.split()
    for allele in alleles.split(","):
        print(donor, allele, sep="\t")
