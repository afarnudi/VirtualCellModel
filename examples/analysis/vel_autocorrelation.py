# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
import numpy as np
import sys

from Python_functions.import_vel_force import autocorr
from Python_functions.import_vel_force import read_file

input_path = '../../../Membrae/DerivedData/Membrae/Build/Products/Debug/Results/'
output_path = input_path

filename = 'chromo2019_11_09_time_10_23'+'_vels_forces.txt'
filename = input_path + str(sys.argv[1])
filename = filename + '.pdb'
data, time =read_file(filename)

mem_num_nodes=162
time_shift=1000
dt=0


vx = data[time_shift:,1000+mem_num_nodes,0]
vy = data[time_shift:,1000+mem_num_nodes,1]
vz = data[time_shift:,1000+mem_num_nodes,2]

V=np.sqrt(np.square(vx) + np.square(vy) + np.square(vz))


if dt !=0:
    cor = autocorr(vx)
    plt.scatter(time[time_shift:time_shift+dt],cor[:dt], label='V_x')
    cor = autocorr(vy)
    plt.scatter(time[time_shift:time_shift+dt],cor[:dt], label='V_y')
    cor = autocorr(vz)
    plt.scatter(time[time_shift:time_shift+dt],cor[:dt], label='V_z')
else:
    cor = autocorr(vx)
    plt.scatter(time[time_shift:time_shift+cor[:].size],cor[:], label='V_x')
    cor = autocorr(vy)
    plt.scatter(time[time_shift:time_shift+cor[:].size],cor[:], label='V_y')
    cor = autocorr(vz)
    plt.scatter(time[time_shift:time_shift+cor[:].size],cor[:], label='V_z')

plt.title("AutoCorrelation")
plt.xlabel("time in Ps")
plt.ylabel("C(t)")
plt.grid(True)
plt.legend()
plt.show()