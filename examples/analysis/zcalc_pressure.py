"""
# -*- coding: utf-8 -*-
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
import numpy as np

filename = 'chromo2019_11_09_time_10_23'

velname = filename + '_vels_forces.txt'
pdbname = filename + '.pdb'


from Python_functions.import_vel_force import read_file
from Python_functions.import_pdb       import import_trajectories_from_pdb

vel_data,  time = read_file(velname)
traj_data, time = import_trajectories_from_pdb(pdbname)


ts, atoms, inds = traj_data.shape

mass = np.ones(atoms)
mem_nodes = 162
for i in range(mem_nodes):
    mass[i]=200

mass = np.asarray([mass])
vs = vel_data[:,:,:3]
fs = vel_data[:,:,3:-1]

unit1=24 # 1Kcal = 24 gr.ang^2/Ps^2.mol
p=np.sum(np.sum(np.square(vs),axis=2)*mass, axis=1) + unit1*np.sum(np.sum(traj_data*fs, axis=2), axis=1)
p = p/3


from Python_functions.gyration_radius import gyration
rg = gyration(traj_data, time)
rg = rg * 1.8

rg3= 4*np.pi*rg*rg*rg/3
p=p/rg3

plt.plot(time,rg, label='Rg mem')

plt.plot(time,p, label='Pressure')
plt.legend()
plt.grid(True)
plt.xlabel('time (Ps)')
plt.show()

