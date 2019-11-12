# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
import numpy as np



filename = 'CMs/'+ 'chromo2019_11_09_time_10_23'+'3'+'.txt'
contact_matrix = np.loadtxt(fname = filename, delimiter = ' ')

fig, ax = plt.subplots()


ax.matshow(contact_matrix , cmap=plt.cm.Blues)
