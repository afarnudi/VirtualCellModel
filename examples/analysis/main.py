# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
import numpy as np

from Python_functions.import_pdb      import import_trajectories_from_pdb
from Python_functions.import_vel_force import read_file
#from Python_functions.gyration_radius import gyration
from Python_functions.RMSD            import RMSD
from Python_functions.R2N             import R2N
import sys

input_path = '../../../Membrae/DerivedData/Membrae/Build/Products/Debug/Results/'
output_path = input_path

pdbname = 'chromo2019_11_18_time_23_15'+'.pdb'
pdbname = 'chromo2019_11_22_time_11_40'+'.pdb'

pdbname = 'chromo2019_11_30_time_08_38.pdb'

filename = input_path + str(sys.argv[1]) 
velname = filename + '_vels_forces.txt'
pdbname = filename + '.pdb'

data, time = import_trajectories_from_pdb(pdbname)
vel_data,  time = read_file(velname)

fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(10,4))

m=10.

energy = np.sum(np.sum(np.square(vel_data[:,:,:3]), axis=2), axis=1)
print(vel_data.shape) 

energy *= 0.5*m 
num_particles = len(vel_data[0])
avgvel = np.mean(np.sum(np.sum(np.square(vel_data[2000:,:,:3]), axis=2), axis=1), axis=0) 
avgvel = m*0.5*avgvel

avgvels = np.ones_like(time)
avgvels.fill(avgvel) 

T = avgvel*120.273/(num_particles*3*.5)

ax2.plot(time, energy)
#ax2.plot(time, energy)
ax2.plot(time[2000:], avgvels[2000:], label="T = "+str(T))
ax2.legend()


mem_nodes = 162
ini_t = 400
fin_t = 800


RMSD(data[ini_t:fin_t,mem_nodes:1000+mem_nodes,:]     , time[ini_t:fin_t], 'ch0',ax1)
RMSD(data[ini_t:fin_t,mem_nodes+1000:2000+mem_nodes,:], time[ini_t:fin_t], 'ch1',ax1)
RMSD(data[ini_t:fin_t,mem_nodes+2000:3000+mem_nodes,:], time[ini_t:fin_t], 'ch2',ax1)
RMSD(data[ini_t:fin_t,mem_nodes+3000:,:],               time[ini_t:fin_t], 'ch3',ax1)


fig.subplots_adjust(wspace=.3)


plt.show()
