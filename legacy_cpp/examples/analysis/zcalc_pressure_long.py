"""
# -*- coding: utf-8 -*-
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
import numpy as np
import sys

input_path = '../../../Membrae/DerivedData/Membrae/Build/Products/Debug/Results/'
output_path = input_path

filename = 'chromo2019_11_21_time_10_33'
filename = 'chromo2019_11_22_time_11_40'
filename = 'chromo2019_11_23_time_07_29'

filename = input_path + str(sys.argv[1])

velname = filename + '_vels_forces.txt'
pdbname = filename + '.pdb'

from Python_functions.import_vel_force import read_file
from Python_functions.import_pdb       import import_trajectories_from_pdb

vel_data,  time = read_file(velname)
traj_data, time = import_trajectories_from_pdb(pdbname)
ts, atoms, inds = traj_data.shape

print(vel_data.shape)
mem_nodes = 162

mass = np.ones(atoms)
for i in range(mem_nodes):
    mass[i]=45

from Python_functions.temperature import calc_temp



fig, axes = plt.subplots(2, 2,figsize=(10,8))


from Python_functions.gyration_radius import gyration

memrg = gyration(traj_data[:,:mem_nodes,:], time)
memrg3= 4*np.pi*memrg*memrg*memrg/3

axes[0,0].plot(time, vel_data[:,0,6]/1000., label='calculated volume')#1000 Ang^3 == 1 nm^3
axes[0,0].plot(time, memrg3/1000.0,'.',ms=1, label='gyration volume', color='orange')
axes[0,0].plot(time,memrg, label='mem', color='orange')
axes[0,0].set_title('Radius of Gyration')
# =============================================================================
from Python_functions.pressure import calc_pressure
p = calc_pressure(traj_data, vel_data, mass)
p=p/memrg3

from Python_functions.pressure import calc_chromo_pressure
pc = calc_chromo_pressure(traj_data[:,mem_nodes:,:], vel_data[:,mem_nodes:,:])
pc = pc/memrg3

axes[0,1].plot(time,p, label='mem + chromos', color='green')
axes[0,1].set_title('Pressure')
#axes[0,1].plot(time,pc, label='chromo Pressure', color='blue')
# =============================================================================

ch0rg = gyration(traj_data[:,mem_nodes:1000 + mem_nodes,:], time)
ch1rg = gyration(traj_data[:,1000 + mem_nodes:2000 + mem_nodes,:], time)
ch2rg = gyration(traj_data[:,2000 + mem_nodes:3000 + mem_nodes,:], time)
ch3rg = gyration(traj_data[:,3000 + mem_nodes:,:], time)

chrg = gyration(traj_data[:,mem_nodes:,:], time)



temperature = calc_temp(traj_data, vel_data, mass)
axes[1,0].plot(time,temperature,'r', label='T = 300')
avgT = np.mean(temperature[1000:])
averageT = np.zeros_like(temperature[1000:])
averageT.fill(avgT)
strK = 'K = '+ str(avgT/300.)
axes[1,0].plot(time[1000:], averageT,'b-', label=strK[:9] )



ch_rdius = 1
mem_radius = 6.33
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


