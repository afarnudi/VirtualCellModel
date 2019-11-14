"""
# -*- coding: utf-8 -*-
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
import numpy as np

filename = 'chromo2019_11_14_time_11_41'

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
init_t = 0



d = squareform( pdist(data[init_t, :, :]) ) < cutoff
d = d.astype(int)

title = "CM: "+filename
label = "\ntime: "+str(init_t)+" to "+str(fs)+"  cut off: "+str(cutoff)+" Ang"
ax.set_title(title)
ax.set_xlabel(label)
im = ax.matshow(d , cmap=plt.cm.Blues)
def update(i):
    global d
    if i==0:
        d = squareform( pdist(data[init_t, :, :]) ) < cutoff
        d = d.astype(int)
    else :
        b = squareform( pdist(data[init_t+i, :, :]) ) < cutoff
        d += b.astype(int)
    return d

def animate(i):
    im.set_data(update(i))
    return

ani = matplotlib.animation.FuncAnimation(fig, animate, frames = fs-init_t)
ani.save(filename+'_cm.gif', writer='imagemagick', fps=fs//10)


