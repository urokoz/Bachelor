#!/usr/bin/env python

import numpy as np

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

def score_cores(core1, core2):
    core_matches = 0
    core_blosum = 0
    for a,b in zip(core1, core2):
        core_matches += int(a==b)
        core_blosum += blosum50[a][b]

        core_ident = core_matches/len(core1)*100

    return core_ident, core_blosum


def k_mer(seq1,seq2, k = 9):
    best_ident = 0
    best_blosum = -np.inf
    for i in range(len(seq1)-(k-1)):
        for j in range(len(seq2)-(k-1)):
            core1 = seq1[i:i+k]
            core2 = seq2[j:j+k]
            ident, blosum = score_cores(core1,core2)

            best_ident = max(ident, best_ident)
            best_blosum = max(blosum, best_blosum)

    return best_ident, best_blosum


seq1 = "AKVEIINAGFTLNGVKNVII"
seq2 = "GDAIGISGGSQIWIDHSSLS"

best_ident, best_blosum = k_mer(seq1, seq2)

print(best_ident)
print(best_blosum)
