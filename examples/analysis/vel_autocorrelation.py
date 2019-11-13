# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
import numpy as np

from Python_functions.import_vel_force import autocorr
from Python_functions.import_vel_force import read_file

filename = 'chromo2019_11_09_time_10_23'+'_vels_forces.txt'

data, time =read_file(filename)

mem_num_nodes=162
time_shift=1000
dt=10
vx = data[time_shift:,1000+mem_num_nodes,0]
vy = data[time_shift:,1000+mem_num_nodes,1]
vz = data[time_shift:,1000+mem_num_nodes,2]

V=np.sqrt(np.square(vx) + np.square(vy) + np.square(vz))

cor = autocorr(vx)
plt.scatter(time[time_shift:time_shift+dt],cor[:dt], label='V_x')


cor = autocorr(vy)
plt.scatter(time[time_shift:time_shift+dt],cor[:dt], label='V_y')

cor = autocorr(vz)
plt.scatter(time[time_shift:time_shift+dt],cor[:dt], label='V_z')


plt.title("AutoCorrelation")
plt.xlabel("time in Ps")
plt.ylabel("C(t)")
plt.grid(True)
plt.legend()
plt.show()