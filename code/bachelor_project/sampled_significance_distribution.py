#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

infile = open("Data/sampled_corr_SRC.txt","r")
list1 = [0]*20
list2 = [0]*20
for line in infile:
    line = line.split()
    SRC = float(line[6])
    pval = float(line[-1])
    # print(pval)
    lower = np.min((float(line[11]), float(line[12])))
    upper = np.max((float(line[11]), float(line[12])))

    i = int(np.floor((SRC+1)*10)) if SRC != 1.0 else 19

    list1[i] += int(pval <= 0.05)
    list2[i] += 1

# print(sum(list2))
bins = []
for x, n in zip(list1, list2):
    print(x,n)
    bins.append(0 if n == 0 else x/n)


fig, ax = plt.subplots()

ax.bar(np.linspace(-0.975, 0.975, 20), bins, width = 0.05)
ax.set_title("SRC significance freq")
plt.show()
