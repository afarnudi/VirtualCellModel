# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
#import numpy as np

from Python_functions.import_pdb      import import_trajectories_from_pdb
from Python_functions.gyration_radius import gyration
from Python_functions.RMSD            import RMSD
from Python_functions.R2N             import R2N


#pdbname = 'chromo2019_10_28_time_10_48'+'.pdb'
pdbname = 'chromo2019_11_09_time_10_23'+'.pdb'
#pdbname = 'chromo2019_11_10_time_00_41'+'.pdb'
#data, time = import_trajectories_from_pdb(pdbname, 'mem0')
data, time = import_trajectories_from_pdb(pdbname)

fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(10,4))

mem_nodes = 162

RMSD(data[:,mem_nodes:1000+mem_nodes,:]     , time[:], 'ch0',ax1)
RMSD(data[:,mem_nodes+1000:2000+mem_nodes,:], time[:], 'ch1',ax1)
RMSD(data[:,mem_nodes+2000:3000+mem_nodes,:], time[:], 'ch2',ax1)
RMSD(data[:,mem_nodes+3000:,:],               time[:], 'ch3',ax1)


R2N (data[100:,mem_nodes:1000+mem_nodes,:],      'ch0',ax2)
R2N (data[100:,mem_nodes+1000:2000+mem_nodes,:], 'ch1',ax2)
R2N (data[100:,mem_nodes+2000:3000+mem_nodes,:], 'ch2',ax2)
R2N (data[100:,mem_nodes+3000:,:],               'ch3',ax2)

left, right = ax2.get_xlim()
ax2.set_xlim(50, right)
down, up = ax2.get_ylim()
ax2.set_ylim(10, up)


fig2, ax3 = plt.subplots()
rg = gyration(data, time)
ax3.scatter(time,rg, label='mem0')
ax3.set_title("Radius of Gyration of structures")
ax3.set_xlabel("time (Ps)")
ax3.set_ylabel("Rg")
ax3.grid(True)
ax3.legend()

plt.show()