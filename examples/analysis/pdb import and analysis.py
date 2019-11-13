# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
import numpy as np

def import_trajectories_from_pdb(filenmae):
    time = []
    energy = []
    data=[]
    with open(filenmae, "r") as f:
        for line in f:
            words = line.split()
            if words[0] == 'MODEL':
                xyz=[]
            elif words[0] == 'REMARK':
                time = np.append(time, float(words[2].replace('time=','')) )
                energy = np.append(energy, float(words[4].replace('energy=','')) )
            elif words[0]=='ATOM':
                temp =np.array([]).reshape(0,3)
                for i in [6,7,8]:
                    temp = np.append(temp, float(words[i]) )
                if len(xyz)==0:
                    xyz=np.append(xyz,temp)
                    xyz=[xyz]
                else:
                    xyz=np.concatenate((xyz, [temp]))

            elif words[0]=='ENDMDL':
                if len(data)==0:
                    data=np.asarray([xyz])
                else:
                    data=np.concatenate((data,[xyz]))
    return data, time

def RMSD(d,atoms_begin, atoms_end,t1,t2, lab, ax):
    rmsd=[]
    t=range(t1,t2)
    for i in t:
        rmsd=np.append(rmsd,np.sqrt(np.sum(np.square(d[i,atoms_begin:atoms_end,:]-d[t1,atoms_begin:atoms_end,:]))/(atoms_end-atoms_begin+1)))
    ax.scatter(time[t1:t2],rmsd, label=lab)
    ax.set_title("RMSD of structures")
    ax.set_xlabel("time (Ps)")
    ax.set_ylabel("RMSD")
    ax.grid(True)
    ax.legend()
    

def R2N(d,atoms_begin, atoms_end, t1, t2, lab,ax):
    r2n=[]
    
    t=range(t1,t2)
    N = (atoms_end+1-atoms_begin)//2
    
    for i in t:
        r2=[]
        for j in range(1,N):
            r2=np.append(r2,np.sqrt(np.sum(np.square(d[i,atoms_begin:atoms_end-j,:]-d[i,atoms_begin+j:atoms_end,:]))/(atoms_end-atoms_begin+1-j)))
        r2n=np.append(r2n,np.mean(r2))
    ax.scatter(time[t1:t2],r2n, label=lab)
    ax.set_title("R2 of bps along the chain of structures")
    ax.set_xlabel("time (Ps)")
    ax.set_ylabel("R2(N_bb)")
    ax.grid(True)
    ax.legend()




pdbname = 'chromo2019_10_28_time_10_48.pdb'
data, time = import_trajectories_from_pdb(pdbname)

frame_count = len(time)

fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(10,4))
RMSD(data,0  ,299,0,frame_count, 'ch0',ax1)
RMSD(data,300,599,0,frame_count, 'ch1',ax1)
R2N (data,0  ,299,0,frame_count, 'ch0',ax2)
R2N (data,300,599,0,frame_count, 'ch1',ax2)
