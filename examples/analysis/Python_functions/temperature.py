import matplotlib.pyplot as plt
import numpy as np
import sys


def calc_temp(data, Vels, mass):
    vs = Vels[:,:,:3]
    K=np.sum(np.sum(np.square(vs),axis=2)*mass, axis=1)
    K=K*10 # to sort out the units for the mv^2 == gr/mol * (Ang/picosec)^2 to KJ/mol
    KT = K/(3*data[1].size)
    return KT

