"""
# -*- coding: utf-8 -*-
Spyder Editor

This is a temporary script file.
"""
from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
import sys

filename = 'chromo2019_11_18_time_23_15'
filename = 'chromo2019_11_22_time_11_40'

filename = 'chromo2019_11_23_time_07_29'
filename = str(sys.argv[1]) 

velname = filename + '_vels_forces.txt'
pdbname = filename + '.pdb'

from Python_functions.import_pdb       import import_trajectories_from_pdb


data, time = import_trajectories_from_pdb(pdbname)
print(data.shape)

#data = data[:,2162:,:]

from scipy.spatial.distance import pdist, squareform

fs, dim1, dim2 = data.shape
import matplotlib.animation

fig, ax = plt.subplots()


cutoff = 4

init_t = 250
fin_t  = 2000

eye = np.eye(dim1)

d = squareform( pdist(data[init_t, :, :]) ) < cutoff
d = d.astype(int)

title = "CM: "+filename
label = "\ntime: "+str(init_t)+" to "+str(fs)+"  cut off: "+str(cutoff)+" Ang"
ax.set_title(title)
ax.set_xlabel(label)
#im = ax.matshow(d , cmap=plt.cm.Blues)
im = ax.matshow(d , cmap=plt.cm.jet)
plt.colorbar(im);

time_bin=10

def update(i):
    global d
    if i==0:
        d = squareform( pdist(data[init_t, :, :]) ) < cutoff
        d = d.astype(int)
        d = d.astype(float)
        return d
    else :
        c = np.zeros((dim1,dim1))
        for j in range(1,time_bin+1):
            b = squareform( pdist(data[init_t+(i-1)*time_bin+j, :, :]) ) < cutoff
            c += b.astype(int)
        c = c.astype(bool)
        c = c.astype(int)
        d += c - eye
        return d/np.amax(d)

def animate(i):
    im.set_data(update(i))
    print(i*time_bin+init_t, end='\r')
    sys.stdout.flush()
    return
    

ani = matplotlib.animation.FuncAnimation(fig, animate, frames = (fin_t-init_t)//time_bin, interval=1)
ani.save(filename+'_cm.gif', writer='imagemagick', fps=fs//10)

