#!/usr/bin/env python
import pickle

with open('/home/mathias/Bachelor/Bachelor/code/specialkursus/data/bg_HLA_dict.pkl', 'rb') as f:
    loaded_dict = pickle.load(f)

for key, val in sorted(loaded_dict.items()):
    print(key, len(val))
