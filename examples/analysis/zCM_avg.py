"""
# -*- coding: utf-8 -*-
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
import numpy as np
import sys

pdbname = str(sys.argv[1])+ '.pdb'

from Python_functions.import_pdb       import import_trajectories_from_pdb

def findab(fname):  
    label=[]
    with open(fname, "r") as f:
        for line in f:
            words = line.split()
            if words[0] == 'MODEL':
                pass
            elif words[0] == 'REMARK':
                pass
            elif words[0]=='ATOM':
                label.append(words[2])
            elif words[0]=='ENDMDL':
                return label


data, time = import_trajectories_from_pdb(pdbname)
print(data.shape)


from scipy.spatial.distance import pdist, squareform

fs, dim1, dim2 = data.shape
#import matplotlib.animation

fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(14,6))


cutoff = 2.5
init_t = 1000
fin_t  = 2000

interval = 1

mem=162
chi=0
chf=chi+1000+3000


listlab = findab(pdbname)
labels = np.asarray(listlab)

label_name = 'ch'+str(chi//1000)
As, = np.where(labels==label_name + '0')
Bs, = np.where(labels==label_name + '1')


axnum=7
eye = np.eye(chf-chi)
for i in range(chf-chi-axnum):
    for j in range(1,axnum+1):
        eye[i][i+j]=1
        eye[i+j][i]=1

dim=chf-chi
if axnum>=2:
    count=0
    for i in range(dim-axnum, dim):
        for j in range(1, axnum-count ):
            eye[i][i+j] = 1
            eye[i+j][i] = 1
        count += 1



d = squareform( pdist(data[init_t, mem+chi:mem+chf, :]) ) < cutoff
d = d.astype(int)
d = d.astype(float)
d += np.ones_like(d)

for i in range(init_t,fin_t,interval):
    b = squareform( pdist(data[i, mem+chi:mem+chf, :]) ) < cutoff
    b = b.astype(int)
    d += b# - eye
    d = abs(d)
    print(i,end="\r")


d = np.log(d)
#d = np.nan_to_num(d)
d = d/(np.amax(d))

title = "CM: "+pdbname[:-4]
label = "\ntime: "+str(init_t)+" to "+str(fs)+"  cut off: "+str(cutoff)+" Ang"
ax1.set_title(title)
ax1.set_xlabel(label)
#im = ax.matshow(d , cmap=plt.cm.Blues)
im = ax1.matshow(d , cmap=plt.cm.jet)



dord = np.zeros_like(d)
count = 0
for i in As.astype(int):
    dord[count] = d[i-mem-chi]
    count+=1

for i in Bs.astype(int):
    dord[count] = d[i-mem-chi]
    count+=1

title = "Orderd:"
label = "\ntime: "+str(init_t)+" to "+str(fs)+"  cut off: "+str(cutoff)+" Ang"
ax2.set_title(title)
ax2.set_xlabel(label)
#im = ax.matshow(d , cmap=plt.cm.Blues)
im = ax2.matshow(dord , cmap=plt.cm.jet)
plt.colorbar(im);

plt.show()

