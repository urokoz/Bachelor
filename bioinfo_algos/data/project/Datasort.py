#!/usr/bin/env python
# coding: utf-8

# Python Imports
import numpy as np
import math
import os

# Load file

# Accessing all subdirectories
os.walk("/home/kamille/22125/Exercises/data/Project/PSSM")
allel = next(os.walk('.'))[1]

# Going through all allels (subdirectories)
for i in range(0,len(allel)):
    data_dir = "/home/kamille/22125/Exercises/data/Project/PSSM" + str(allel[i]) + "/"

    # Going through all training files
    for i in range(0,5):
        training_file = data_dir + "f00" + str(i)
        pssm_training_file = data_dir + "pt00" + str(i)

        ## Reading the data
        t_data = open(training_file, 'r')
        t_outfile = open(pssm_training_file, "w")

        # Specifying the binding threshold
        threshold = 1-math.log(500)/math.log(50000)

        for line in t_data:
            IC50 = float(line.split()[1])
            if IC50 > threshold:
                t_sequence = line.split()[0]
                print(t_sequence, file = t_outfile)


        # Closes file
        t_data.close()
        t_outfile.close()

    # Going through all evaluation files
    for i in range(0,5):
        evaluation_file = data_dir + "c00" + str(i)
        pssm_evaluation_file = data_dir + "pe00" + str(i)

        # Reading the data
        e_data = open(evaluation_file, 'r')
        e_outfile = open(pssm_evaluation_file, "w")

        for line in e_data:
            e_sequence = line.split()[:2]
            print(str(e_sequence[0]) + " " + str(e_sequence[1]), file = e_outfile)

        # Closes file
        e_data.close()
        e_outfile.close()
