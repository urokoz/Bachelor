#!/usr/bin/env python
# coding: utf-8

# # A very first Python program
#
# In this program we will:
#
# 1.   Load a peptide-target list using NumPy.
# 2.   Keep only 9mer peptides.
# 3.   Discard peptides using a target value threshold.
# 4.   Print the filtered result.
# 5.   Download the notebook as a Python file (.py).
# 6.   Implement a command line parser in the .py file.
#
# In the code below, some parts have been blanked out with XX. Fill in the codee these places.
#
# Note that if you are using python2 you might want to remove the "()" parentheses in the print commands

# ## Import NumPy

# In[21]:


import numpy as np
from argparse import ArgumentParser


# ## DEFINE THE PATH TO YOUR COURSE DATA DIRECTORY

# In[5]:


data_dir = "/home/mathias/bioinfo_algos/data/"


# ## Load peptides-targets data

# **Specify file path**

# In[6]:

parser = ArgumentParser(description="A very first Python program")
parser.add_argument("-t", action="store", dest="threshold", type=float, default=0.5, help="Target value filtering threshold (default: 0.5)")
parser.add_argument("-f", action="store", dest="peptides_targets_file", type=str, help="Peptides-Targets file")
parser.add_argument("-l", action="store", dest="peptide_length", type=int, help="Desired length of the peptides")

args = parser.parse_args()
threshold = args.threshold
peptides_targets_file = args.peptides_targets_file
peptide_length = args.peptide_length


# **Load the file with the numpy text parser *np.loadtxt* and reshape it into a numpy array of shape (-1, 2)**
#
# *This means "give me all the rows you want, but I demand only 2 columns", and ensures the [PEPTIDE, TARGET_VALUE] structure we want*

# In[8]:


peptides_targets = np.loadtxt(peptides_targets_file, dtype=str)


# **Check the shape of your newly created numpy array**

# In[9]:


print(peptides_targets.shape)


# ## Store peptides in vector

# **Fish out the peptides from the loaded file using array indexing and slicing**

# In[10]:


peptides = peptides_targets[:, 0]


# **Check that *peptides* has the same data type as *peptides_targets***

# In[11]:


print(type(peptides), type(peptides_targets))


# ## Store targets in vector

# **Fish out the target values from the loaded file using array indexing and slicing**
#
# *Remember that we used a text parser to load the data. So, we need to cast this value to float somehow*
#
# *Hint, use a similar command as the one for fishing out the peptides, and add .astype(float)*

# In[12]:


targets = peptides_targets[:, 1].astype(float)


# ## Keep 9mers only

# **Declare two Python lists to store peptides and targets**
#

# In[13]:


peptides_9mer = []
targets_9mer = []


# **Iterate over the elements of the peptides list and keep peptides with length == 9 using the .append command**

# In[14]:


for i in range(0, len(peptides)):

    if len(peptides[i]) == peptide_length:

        peptides_9mer.append(peptides[i])

        targets_9mer.append(targets[i])


# ## Remove peptides with target value < threshold



# **Declare python list to store the indexes of the elements to be removed**

# In[16]:


to_remove = []


# **Iterate over the 9mer peptides, check which target values < threshold, and store the indexes in the to_remove  array**

# In[17]:


for i in range(0, len(peptides_9mer)):

        if targets_9mer[i] < threshold:

            to_remove.append([i])


# **Use the *delete* NumPy function to remove the peptides**

# In[18]:


peptides_9mer_t = np.delete(peptides_9mer, to_remove)
targets_9mer_t = np.delete(targets_9mer, to_remove)


# **Check that no elements with target < threshold are present in the target values array**

# In[19]:


error = False

for i in range(0, len(peptides_9mer_t)):

        if targets_9mer_t[i] < threshold:

            error = True

            break

if error:

    print("Something went wrong")

else:

    print("Success")


# ## Print the final, filtered peptide-target pairs

# **Ensure that this output is consistent with the data filtering steps you have made!**

# In[20]:


for i in range(0, len(peptides_9mer_t)):

    print(peptides_9mer_t[i], targets_9mer_t[i])


# ## Adding a command line parser

# In[ ]:


################################
# ADDING A COMMAND LINE PARSER #
################################

# For this step, we need first to import an argument parser
# to do this, add the following line just below the numpy import:

# from argparse import ArgumentParser

# We will now create an argument parser that will receive as arguments two values
# 1) the peptides-targets file to open (-f option)
# 2) the threshold to be applied in the target value filtering step (-t option)
# To achieve this, add the following lines below the ArgumentParser import line:


# parser = ArgumentParser(description="A very first Python program")
# parser.add_argument("-t", action="store", dest="threshold", type=float, default=0.5, help="Target value filtering threshold (default: 0.5)")
# parser.add_argument("-f", action="store", dest="peptides_targets_file", type=str, help="Peptides-Targets file")
# args = parser.parse_args()
# threshold = args.threshold
# peptides_targets_file = args.peptides_targets_file


# After adding these lines, you will now be able to call this python program
# from the terminal while specifying these arguments:


# python Python_Intro.py -t some_threshold -f file_with_peptides_and_targets

# Note you can also parse switches with the ArgumentParser, i.e
# parser.add_argument('-w', action='store_true', default=False, dest='sequence_weighting', help='Use sequence weighting')


# REMEMBER!
# 1) The argument parser needs to be declared on the beginning of the script, right after the imports
# 2) In order for this program to work properly after adding the parser, you must now comment or delete
#    the previous declarations of the variables "threshold" and "peptides_target_file"


# # Now download this as a Python file (File -> Save as .py), and continue working with this file offline
#
# ## Modify the code to include a command line parser to allow the program to accept three options
#
# 1. -t THRESHOLD          Target value filtering threshold (default: 0.5)
#
# 2. -f PEPTIDES_TARGETS_FILE Peptides-Targets file
#
# 3. -l PEPLEN_THRESHOLD   Peptide length (default: 9)

# In[ ]:
