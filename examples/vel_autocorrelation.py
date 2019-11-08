# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
import numpy as np


def autocorr(x):
    m = np.mean(x)
    N = x.size//2
    results = np.zeros(N)
    results[0] = 1
    for tau in range(1,N):
        for i, j in zip(x[:-tau],x[tau:]):
            results[tau] += (i-m)*(j-m)
        results[tau] = results[tau]/(np.sum(np.square(x)) -2*m*np.sum(x) - m*m )
    return results
    
    

filename = 'chromo2019_10_27_time_11_18vels.txt'
data = np.loadtxt(fname = filename, delimiter = '\t')


time = data[:,0]
time = time[::3]



vs = data[:,1]
vx = vs[::3]
vy = vs[1::3]
vz = vs[2::3]


dtime =30
shift_data = 600

cor = autocorr(vx[shift_data:])
plt.scatter(time[:dtime],cor[:dtime], label='V_x')

cor = autocorr(vy[shift_data:])
plt.scatter(time[:dtime],cor[:dtime], label='V_y')

cor = autocorr(vz[shift_data:])
plt.scatter(time[:dtime],cor[:dtime], label='V_z')


plt.title("AutoCorrelation")
plt.xlabel("time in Ps")
plt.ylabel("C(t)")
plt.grid(True)
plt.legend()
