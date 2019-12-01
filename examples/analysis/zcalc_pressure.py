"""
# -*- coding: utf-8 -*-
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
import numpy as np
import sys

filename = 'chromo2019_11_21_time_10_33'
filename = 'chromo2019_11_22_time_11_40'
filename = 'chromo2019_11_23_time_07_29'
filename = str(sys.argv[1])
#filename = 'chromo2019_11_18_time_11_11'

velname = filename + '_vels_forces.txt'
pdbname = filename + '.pdb'

from Python_functions.import_vel_force import read_file
from Python_functions.import_pdb       import import_trajectories_from_pdb

vel_data,  time = read_file(velname)
traj_data, time = import_trajectories_from_pdb(pdbname)
ts, atoms, inds = traj_data.shape


def calc_pressure(data, Vels, mnodes=162):
    mass = np.ones(atoms)

    for i in range(mnodes):
        mass[i]=200
    mass = np.asarray([mass])

    vs = Vels[:,:,:3]
    fs = Vels[:,:,3:-1]

    unit1=4.182 
    unit2=0.174 # 1Kcal = 24 gr.ang^2/Ps^2.mol
    p=unit2*np.sum(np.sum(np.square(vs),axis=2)*mass, axis=1) + unit1*np.sum(np.sum(data*fs, axis=2), axis=1)
    return p/3


def calc_chromo_pressure(data, Vels, mnodes=162):
    mass = np.asarray([np.ones(atoms-mnodes)])
    vs = Vels[:,:,:3]
    fs = Vels[:,:,3:-1]
    unit1=24 # 1Kcal = 24 gr.ang^2/Ps^2.mol
    p=np.sum(np.sum(np.square(vs),axis=2)*mass, axis=1) + unit1*np.sum(np.sum(data*fs, axis=2), axis=1)
    return p/3

fig, axes = plt.subplots(2, 2,figsize=(10,8))

mem_nodes = 162


from Python_functions.gyration_radius import gyration

memrg = gyration(traj_data[:,:mem_nodes,:], time)
memrg3= 4*np.pi*memrg*memrg*memrg/3


axes[0,0].plot(time,memrg, label='mem', color='orange')
axes[0,0].set_title("Membrane Radius")
axes[0,0].set_ylabel("R_g in Angstrom")

p = calc_pressure(traj_data, vel_data)
p=p/memrg3

pc = calc_chromo_pressure(traj_data[:,mem_nodes:,:], vel_data[:,mem_nodes:,:])
pc = pc/memrg3

axes[0,1].plot(time,p, label='mem + chromos', color='green')
axes[0,1].set_title('Pressure')
axes[0,1].set_ylabel('P (1.6*GPa)')

ch0rg = gyration(traj_data[:,mem_nodes:1000 + mem_nodes,:], time)
ch1rg = gyration(traj_data[:,1000 + mem_nodes:2000 + mem_nodes,:], time)
ch2rg = gyration(traj_data[:,2000 + mem_nodes:3000 + mem_nodes,:], time)
ch3rg = gyration(traj_data[:,3000 + mem_nodes:,:], time)

axes[1,0].plot(time,ch0rg, label='ch0')
axes[1,0].plot(time,ch1rg, label='ch1')
axes[1,0].plot(time,ch2rg, label='ch2')
axes[1,0].plot(time,ch3rg, label='ch3')
axes[1,0].set_title("Radius of Gyration")
axes[1,0].set_ylabel("R_g in Angstrom")

ch_rdius = 1
mem_radius = 4
ch_volume = ch_rdius*ch_rdius*ch_rdius*4000
eff_radius = memrg - mem_radius
Mem_volume = eff_radius*eff_radius*eff_radius
volume_fraction = ch_volume/Mem_volume

axes[1,1].plot(time, volume_fraction, label='chromos')
axes[1,1].set_title('Volume Fraction')
# =============================================================================



for ax in [axes[0,0], axes[0,1], axes[1,0], axes[1,1]]:
    ax.legend()
    ax.grid(True)
    ax.set_xlabel('time (Ps)')
fig.subplots_adjust(hspace=.5)
plt.show()


