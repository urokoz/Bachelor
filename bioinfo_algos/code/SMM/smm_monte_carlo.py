#!/usr/bin/env python
# coding: utf-8

# # SMM with Monte Carlo

# ## Python Imports

# In[1]:


import numpy as np
import random
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import argparse

## Argument parsing:

parser = argparse.ArgumentParser(description="Hobohm1 program")
parser.add_argument("-l", action="store", dest="LAMB", type=float, default = 0.01, help="Lambda (default: 0.01)")
parser.add_argument("-t", action="store", dest="TRAINING_FILE", type=str, help="File with training data")
parser.add_argument("-e", action="store", dest="EVALUATION_FILE", type=str, help="File with evaluation data")
parser.add_argument("-s", action="store", dest="SEED", type=int, default=1, help="Seed for random numbers (default 1)")
parser.add_argument("-i", action="store", dest="ITERS", type=int, default=1000, help="Number of epochs to train (default 1000)")
parser.add_argument("-Ts", action="store", dest="T_I", type=float, default=0.01, help="Start Temp (0.01)")
parser.add_argument("-Te", action="store", dest="T_F", type=float, default=0.000001, help="End Temp (0.000001)")
parser.add_argument("-nT", action="store", dest="T_STEPS", type=int, default=10, help="Number of T steps (10)")


args = parser.parse_args()
lamb = args.LAMB
training_file = args.TRAINING_FILE
evaluation_file = args.EVALUATION_FILE
seed = args.SEED


# ## Data Imports

# ## DEFINE THE PATH TO YOUR COURSE DIRECTORY

# In[2]:


data_dir = "/home/mathias/bioinfo_algos/data/"


# ## Define run time parameters
# temperature vector
T_i = args.T_I
T_f = args.T_F
T_steps = args.T_STEPS
T_delta = (T_f - T_i) / T_steps
T = np.linspace(T_i,T_f,T_steps )

# iterations
iters = args.ITERS

# regularization lambda
#lamb = 1
#lamb = 0.00001
#lamb = 10


# ### Training Data

# In[4]:


#training_file = data_dir + "SMM/A0201_training"
#training_file = data_dir + "SMM/A2403_training"
training = np.loadtxt(training_file, dtype=str)


# ### Evaluation Data

# In[5]:


#evaluation_file = data_dir + "SMM/A0201_evaluation"
#evaluation_file = data_dir + "SMM/A2403_evaluation"
evaluation = np.loadtxt(evaluation_file, dtype=str)


# ### Alphabet

# In[6]:


alphabet_file = data_dir + "Matrices/alphabet"
alphabet = np.loadtxt(alphabet_file, dtype=str)


# ### Sparse Encoding Scheme

# In[7]:


sparse_file = data_dir + "Matrices/sparse"
_sparse = np.loadtxt(sparse_file, dtype=float)
sparse = {}

for i, letter_1 in enumerate(alphabet):

    sparse[letter_1] = {}

    for j, letter_2 in enumerate(alphabet):

        sparse[letter_1][letter_2] = _sparse[i, j]


# ## Peptide Encoding

# In[8]:


def encode(peptides, encoding_scheme, alphabet):

    encoded_peptides = []

    for peptide in peptides:

        encoded_peptide = []

        for peptide_letter in peptide:

            for alphabet_letter in alphabet:

                encoded_peptide.append(encoding_scheme[peptide_letter][alphabet_letter])

        encoded_peptides.append(encoded_peptide)

    return np.array(encoded_peptides)


# ## Error Function

# In[9]:


def cumulative_error(peptides, y, lamb, weights):

    error = 0

    for i in range(0, len(peptides)):

        # get peptide
        peptide = peptides[i]

        # get target prediction value
        y_target = y[i]

        # get prediction
        y_pred = np.dot(peptide, weights)

        # calculate error
        error += 1.0/2 * (y_pred - y_target)**2

    gerror = error + lamb*np.dot(weights, weights)
    error /= len(peptides)

    return gerror, error


# ## Predict value for a peptide list

# In[10]:


def predict(peptides, weights):

    pred = []

    for i in range(0, len(peptides)):

        # get peptide
        peptide = peptides[i]

        # get prediction
        y_pred = np.dot(peptide, weights)

        pred.append(y_pred)

    return pred


# ## Calculate MSE between two vectors

# In[11]:


def cal_mse(vec1, vec2):

    mse = 0

    for i in range(0, len(vec1)):
        mse += (vec1[i] - vec2[i])**2

    mse /= len(vec1)

    return( mse)



# ## Main Loop

# In[13]:


# Random seed
np.random.seed( 1 )

# peptides
peptides = training[:, 0]
peptides = encode(peptides, sparse, alphabet)

# target values
y = np.array(training[:, 1], dtype=float)

#evaluation peptides
evaluation_peptides = evaluation[:, 0]
evaluation_peptides = encode(evaluation_peptides, sparse, alphabet)

#evaluation targets
evaluation_targets = np.array(evaluation[:, 1], dtype=float)

# weights
input_dim  = len(peptides[0])
output_dim = 1
w_bound = 0.1
weights = np.random.uniform(-w_bound, w_bound, size=input_dim)

# error plots
gerror_plot = []
mse_plot = []
train_mse_plot = []
eval_mse_plot = []
train_pcc_plot = []
eval_pcc_plot = []

perturbation_value = 0.1
# The scale variable defines the scale of the changes in de. This is defined empirically for each problem investigated
scale = 1.0/100

number_of_tries = 0
number_of_accepted = 0

# calculate initial error
gerror_initial, mse = cumulative_error(peptides, y, lamb, weights)

# for each temperature
for t in T:

    # for each iteration
    for i in range(0, iters):


        # get 2 random weight indexes
        weight_index_1 = np.random.randint(len(weights))
        weight_index_2 = np.random.randint(len(weights))

        # ensure they are different
        while weight_index_1 == weight_index_2:
            weight_index_2 = np.random.randint(len(weights))

        # store original weight values
        original_weight_1 = weights[weight_index_1]
        original_weight_2 = weights[weight_index_2]


        # apply random perturbation to both weights
        perturbation = np.random.uniform(0, perturbation_value)
        weights[weight_index_1] += perturbation
        weights[weight_index_2] -= perturbation


        # calculate new error
        gerror_new, mse = cumulative_error(peptides, y, lamb, weights)

        # compute error difference
        de = (gerror_new - gerror_initial) * scale

        # compute acceptance probability
        if ( de < 0):
            p = 1
        else:
            p = np.exp(-de/t)

        # throw coin
        coin = np.random.uniform(0.0, 1.0, 1)[0]

        # weight change is accepted
        if ( coin < p ):
            gerror_initial = gerror_new
            gerror_plot.append(gerror_new)
            mse_plot.append(mse)
            number_of_accepted += 1
        # weight change is declined, restore previous weights
        else:
            weights[weight_index_1] = original_weight_1
            weights[weight_index_2] = original_weight_2

        number_of_tries += 1

        # define size of move so that on avarage 50% are accepted
        if number_of_tries == 100:

            if float(number_of_accepted)/number_of_tries > 0.5:
                perturbation_value *= 1.1
            else:
                perturbation_value *= 0.9

            number_of_tries = 0
            number_of_accepted = 0

        # predict on training data
        train_pred = predict(peptides, weights)
        train_mse = cal_mse(y, train_pred)
        train_pcc = pearsonr(y, train_pred)[0]
        train_mse_plot.append(train_mse)
        train_pcc_plot.append(train_pcc)

        # predict on evaluation data
        eval_pred = predict(evaluation_peptides, weights)
        eval_mse = cal_mse(evaluation_targets, eval_pred)
        eval_pcc = pearsonr(evaluation_targets, eval_pred)[0]
        eval_mse_plot.append(eval_mse)
        eval_pcc_plot.append(eval_pcc)

    print ("#t:", t, gerror_new, perturbation_value, train_mse, train_pcc, eval_mse, eval_pcc)


#
# our matrices are vectors of dictionaries
def vector_to_matrix(vector, alphabet):

    rows = int(len(vector)/len(alphabet))

    matrix = [0] * rows

    offset = 0

    for i in range(0, rows):

        matrix[i] = {}

        for j in range(0, 20):

            matrix[i][alphabet[j]] = vector[j+offset]

        offset += len(alphabet)

    return matrix


# ### Matrix to Psi-Blast
def to_psi_blast(matrix):

    # print to user

    header = ["", "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

    print('{:>4} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}'.format(*header))

    letter_order = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

    for i, row in enumerate(matrix):

        scores = []

        scores.append(str(i+1) + " A")

        for letter in letter_order:

            score = row[letter]

            scores.append(round(score, 4))

        print('{:>4} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}'.format(*scores))


# ### Print
matrix = vector_to_matrix(weights, alphabet)
to_psi_blast(matrix)
