#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  9 14:18:50 2019

@author: ali
"""
import numpy as np

def R2N(data, lab, ax):
    """
    Calculates the root mean square distance of the data averaged over the time period.
    inputs:
        data:
            a 3d numpy array: [time step, atom index, coordinates (3)]
        lab:
            label of the data to be displayed in the graph
        ax:
            axis to be ploted on.
    """
    r2n=[]
    N=len(data[0])//2
    for j in range(1,N+1,10):
        r2n=np.append(r2n,np.mean(np.sqrt(np.mean(np.sum(np.square(data[:,:-j,:]-data[:,j:,:]),axis=2),axis=1))))
        
    ax.scatter(range(1,N+1,10),r2n, label=lab, alpha=0.3, s=15)
    ax.set_title("$R^2$ of $bps$ separated along the chain\n(averaged over time)")
    ax.set_xlabel("$N_{bp}$")
    ax.set_ylabel("$R^2(N_{bp})$")
    ax.grid(True, which='both')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.legend()
    print('{}: R2N finished.'.format(lab))
    return r2n
