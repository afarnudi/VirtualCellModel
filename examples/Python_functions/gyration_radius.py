#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  9 17:57:36 2019

@author: ali
"""
import numpy as np
def gyration(data,time):
    return np.sqrt(np.mean(np.sum(np.square(data-np.broadcast_to(np.mean(data,axis=1).reshape(len(data),1,3),(len(data),len(data[0]),3))),axis=2),axis=1))