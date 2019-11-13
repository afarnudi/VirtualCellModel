#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  9 14:17:41 2019

@author: ali
"""
import numpy as np

def concat_data(words, xyz):
    temp =np.array([]).reshape(0,3)
    for i in [6,7,8]:
        try:
            temp = np.append(temp, float(words[i]) )
        except:
            split = words[i].split('-')
            if len(split)== 3:
                temp = np.append(temp, -float(split[1]) )
                temp = np.append(temp, -float(split[2]) )
            elif len(split)== 4:
                temp = np.append(temp, -float(split[1]) )
                temp = np.append(temp, -float(split[2]) )
                temp = np.append(temp, -float(split[3]) )
                break
            else:
                temp = np.append(temp, float(split[0]) )
                temp = np.append(temp, -float(split[1]) )
            if i == 6:
                temp = np.append(temp, float(words[7]) )
                break
            else:
                break
    if len(xyz)==0:
        xyz=np.append(xyz,temp)
        xyz=[xyz]
    else:
        xyz=np.concatenate((xyz, [temp]))
    return xyz

def import_trajectories_from_pdb(filename, label = None):
    name = filename
#    filename = '../'+filename
    
    try:
        data = np.load(filename[:-3]+'npy')
        time = np.load(filename[:-4]+'_t.npy')
        
        print('{} imported successfully from the "npy"s.'.format(name))
        return data, time
        
    except:
        print('No npy files located. importing from the pdb.')
    
    time = []
    energy = []
    data=[]
    with open(filename, "r") as f:
        for line in f:
            words = line.split()
            if words[0] == 'MODEL':
                xyz=[]
            elif words[0] == 'REMARK':
                time = np.append(time, float(words[2].replace('time=','')) )
                energy = np.append(energy, float(words[4].replace('energy=','')) )
            elif words[0]=='ATOM':
                if label is None:
                    xyz=concat_data(words,xyz)
                elif len(label) == 3:
                    lab = words[2]
                    if lab[:3] == label:
                        xyz=concat_data(words,xyz)
                    else:
                        continue
                else:
                    if words[2] == label:
                        xyz=concat_data(words,xyz)
                    else:
                        continue
                

            elif words[0]=='ENDMDL':
                if len(data)==0:
                    data=np.asarray([xyz])
                else:
                    data=np.concatenate((data,[xyz]))
    np.save(filename[:-4],data)
    np.save(filename[:-4]+'_t',time)
    print('{} imported successfully.'.format(name))
    return data, time