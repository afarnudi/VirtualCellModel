"""
# -*- coding: utf-8 -*-
Spyder Editor

This is a temporary script file.
"""

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


from scipy.spatial.distance import pdist, squareform

fs, dim1, dim2 = data.shape
import matplotlib.animation

fig, ax = plt.subplots()


cutoff = 2.5
init_t = 1000
fin_t  = 1100

eye = np.eye(dim1)

d = squareform( pdist(data[init_t, :, :]) ) < cutoff
d = d.astype(int)
d = d.astype(float)
#d += np.ones_like(d)


title = "CM: "+filename
label = "\ntime: "+str(init_t)+" to "+str(fin_t)+"  cut off: "+str(cutoff)+" Ang"
ax.set_title(title)
ax.set_xlabel(label)
#im = ax.matshow(d , cmap=plt.cm.Blues)
im = ax.matshow(d , cmap=plt.cm.jet)
plt.colorbar(im);

def update(i):
    global d
    if i==0:
        d = squareform( pdist(data[init_t, :, :]) ) < cutoff
        d = d.astype(int)
        d = d.astype(float)
        d += np.ones_like(d)
        d = d/np.amax(d)
        return d
    else :
        b = squareform( pdist(data[init_t+i, :, :]) ) < cutoff
        b = b.astype(int)
        d += b #- eye
        print(i+init_t,end='\r')
        return d#np.log(d)/np.amax(np.log(d))

def animate(i):
    im.set_data(update(i))
    return
    

ani = matplotlib.animation.FuncAnimation(fig, animate, frames = fin_t-init_t, interval=2)
ani.save(filename+'_cm.gif', writer='imagemagick', fps=(fin_t-init_t)//10)

