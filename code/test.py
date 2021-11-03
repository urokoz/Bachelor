#!/usr/bin/env python

import numpy as np

load_pep_ker_dict(filename):

    pep_ker_dict = dict()
    pep_ker_file = open(filename, "r")
    for line in pep_ker_file:
        line = line.split()
        index = frozenset([line[0], line[1]])

        pep_ker_dict[index] = line[2]

    return pep_ker_dict
