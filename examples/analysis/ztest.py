"""
# -*- coding: utf-8 -*-
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
import numpy as np

filename = 'chromo2019_11_13_time_15_35'

velname = filename + '_vels_forces.txt'
pdbname = filename + '.pdb'


from Python_functions.import_vel_force import read_file
from Python_functions.import_pdb       import import_trajectories_from_pdb

#vel_data,  time = read_file(velname)
data, time = import_trajectories_from_pdb(pdbname)

from scipy.spatial.distance import pdist, squareform

d= squareform( pdist(data[151, 162: , :]) )
ind1 , ind2 = np.where(d<2)
unique = (ind1<ind2)
ind1 = ind1[unique]
ind2 = ind2[unique]

for i1, i2 in zip(ind1,ind2):
    if abs(i1-i2)!=1:
        print(i1,i2)
d = d<2 
fig, ax = plt.subplots()

ax.matshow(d , cmap=plt.cm.Blues)
plt.show()
print(d.shape)
