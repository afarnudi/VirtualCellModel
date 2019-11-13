#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  9 14:18:50 2019

@author: ali
"""
import numpy as np

def R2N(data, lab, ax):
    r2n=[]
    N=len(data[0])//2
    for j in range(1,N+1,10):
        r2n=np.append(r2n,np.mean(np.sqrt(np.mean(np.sum(np.square(data[:,:-j,:]-data[:,j:,:]),axis=2),axis=1))))
        
    ax.scatter(range(1,N+1,10),r2n, label=lab)
    ax.set_title("R2 of bps along the chain of structures (averaged over time)")
    ax.set_xlabel("N_bb")
    ax.set_ylabel("R2(N_bb)")
    ax.grid(True, which='both')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.legend()
    print('{}: R2N finished.'.format(lab))
