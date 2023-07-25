import matplotlib.pyplot as plt
import numpy as np
import sys

def calc_pressure(data, vels, mass, unit=24):
    """
    unit=24 # 1Kcal = 24 gr.ang^2/Ps^2.mol
    """

    vs = vels[:,:,:3]
    fs = vels[:,:,3:6]
    return (np.sum(np.sum(np.square(vs),axis=2)*mass, axis=1) + unit*np.sum(np.sum(data*fs, axis=2), axis=1) )/3.

def calc_chromo_pressure(data, Vels, mass=1, unit=24):
    """
    unit :: 1Kcal = 24 gr.ang^2/Ps^2.mol
    """
    vs = Vels[:,:,:3]
    fs = Vels[:,:,3:6]
    return ( np.sum(np.sum(np.square(vs),axis=2)*mass, axis=1) + unit*np.sum(np.sum(data*fs, axis=2), axis=1) )/3.
