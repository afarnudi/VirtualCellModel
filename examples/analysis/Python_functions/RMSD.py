#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  9 14:18:13 2019

@author: ali
"""
import numpy as np
def RMSD(data, time, lab, ax):
    """
    Calculate the root mean square distance of the nodes relative to the initial frame.
    """
    rmsd=np.sqrt(np.mean(np.sum(np.square(data[:,:,:]-data[0,:,:]),axis=2),axis=1))
    ax.plot(time,rmsd, label=lab)
    ax.set_title("RMSD of structures")
    ax.set_xlabel("time (Ps)")
    ax.set_ylabel("RMSD")
    ax.grid(True)
    ax.legend()
    print('{}: RMSD finished.'.format(lab))
