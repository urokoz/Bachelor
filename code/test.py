#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

infile = open("Data/sampled_corr.txt","r")
list1 = [0]*11
list2 = [0]*11
for line in infile:
    line = line.split()
    SRC = float(line[6])
    pval = float(line[-1])
    # print(pval)
    lower = np.min((float(line[11]), float(line[12])))
    upper = np.max((float(line[11]), float(line[12])))

    i = int(np.floor((SRC+1)*5))

    list1[i] += int(pval <= 0.05)
    list2[i] += 1

# print(sum(list2))
bins = []
for x, n in zip(list1, list2):
    print(x,n)
    bins.append(0 if n == 0 else x/n)


fig, ax = plt.subplots()

ax.bar(np.linspace(-1, 1, 11), bins, width = 0.1)
ax.set_title("SRC significance freq")
plt.show()
