# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
import numpy as np

from Python_functions.import_pdb      import import_trajectories_from_pdb
#from Python_functions.gyration_radius import gyration
from Python_functions.RMSD            import RMSD
from Python_functions.R2N             import R2N
import sys

input_path = '../../../Membrae/DerivedData/Membrae/Build/Products/Debug/Results/'
output_path = input_path

soft_pdbs = ['chromo2019_12_06_time_22_53.pdb']
hard_pdbs = ['chromo2019_12_18_time_23_20.pdb','chromo2019_12_19_time_17_48.pdb']

def add_input_path(namelist, inputpath):
    temp=[]
    for name in namelist:
        temp.append(inputpath + name)
    return temp

soft_pdbs=add_input_path(soft_pdbs, input_path)
hard_pdbs=add_input_path(hard_pdbs, input_path)

fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(12,5.5))
fig2, ax3 = plt.subplots(1,1)
mem_nodes = 162
ini_t = int(400/0.4)
fin_t = int(800/0.4)

data, time = import_trajectories_from_pdb(soft_pdbs[0])

sr2n = [R2N (data[ini_t:fin_t,mem_nodes:1000+mem_nodes,:],    None  ,ax2)]
sr2n = np.concatenate((sr2n, [R2N (data[ini_t:fin_t,mem_nodes+1000:2000+mem_nodes,:],None , ax2)]))
sr2n = np.concatenate((sr2n, [R2N (data[ini_t:fin_t,mem_nodes+2000:3000+mem_nodes,:],None , ax2)]))
sr2n = np.concatenate((sr2n, [R2N (data[ini_t:fin_t,mem_nodes+3000:,:],              None, ax2)]))

data, time = import_trajectories_from_pdb(hard_pdbs[0])

hr2n = [R2N (data[ini_t:fin_t,mem_nodes:1000+mem_nodes,:],      None,ax2)]
hr2n = np.concatenate((hr2n, [R2N (data[ini_t:fin_t,mem_nodes+1000:2000+mem_nodes,:], None,ax2)]))
hr2n = np.concatenate((hr2n, [R2N (data[ini_t:fin_t,mem_nodes+2000:3000+mem_nodes,:], None,ax2)]))
hr2n = np.concatenate((hr2n, [R2N (data[ini_t:fin_t,mem_nodes+3000:,:],               None,ax2)]))

for name in hard_pdbs[1:]:
    data, time = import_trajectories_from_pdb(name)
    hr2n = np.concatenate((hr2n, [R2N (data[ini_t:fin_t,mem_nodes:1000+mem_nodes,:],      None,ax2)]))
    hr2n = np.concatenate((hr2n, [R2N (data[ini_t:fin_t,mem_nodes+1000:2000+mem_nodes,:], None,ax2)]))
    hr2n = np.concatenate((hr2n, [R2N (data[ini_t:fin_t,mem_nodes+2000:3000+mem_nodes,:], None,ax2)]))
    hr2n = np.concatenate((hr2n, [R2N (data[ini_t:fin_t,mem_nodes+3000:,:],               None,ax2)]))

sr2ntomum = np.square(sr2n.mean(axis=0))
sr2ntomum = np.mean(np.square(sr2n),axis=0)
hr2ntomum = np.square(hr2n.mean(axis=0))
hr2ntomum = np.mean(np.square(hr2n),axis=0)

n, N = sr2n.shape
N = N*10
Nbp = []
for  i in range(1,N+1,10):
    Nbp.append(i*0.063)

ax3.errorbar(range(1,N+1,10),sr2n.mean(axis=0)/32.1, label='avg soft', marker='s', c='black', ms=4, yerr=sr2n.std(axis=0), ecolor='r' )

ax3.errorbar(range(1,N+1,10),hr2n.mean(axis=0)/28.8, label='hvg soft', marker='s', c='blue', ms=4, yerr=hr2n.std(axis=0), ecolor='r' )

ax3.set_yscale('log')
ax3.set_xscale('log')
ax3.legend()

sr2ntomum100 = sr2ntomum*0.01
hr2ntomum100 = hr2ntomum*0.01

ax1.errorbar(Nbp,sr2ntomum100, label='Savg-100nm', marker='s', c='black', ms=2, yerr=np.std(np.square(sr2n),axis=0)*0.01, ecolor='g', ls='none', zorder=8 )
ax1.errorbar(Nbp,hr2ntomum100, label='Havg-100nm', marker='s', c='b'    , ms=2, yerr=np.std(np.square(hr2n),axis=0)*0.01, ecolor='g', ls='none', zorder=7 )


sr2ntomum30 = sr2ntomum*0.03*0.03
#ax1.errorbar(Nbp,sr2ntomum30, label='Savg-30nm', marker='s', c='orange', ms=2, yerr=sr2n.std(axis=0)*0.03*0.03, ecolor='g', ls='none', zorder=2, alpha=0.7 )

sr2ntomum280 = sr2ntomum*0.28*0.28
#ax1.errorbar(Nbp,sr2ntomum280, label='Savg-280nm', marker='s', c='m', ms=2, yerr=sr2n.std(axis=0)*0.28*0.28, ecolor='g', ls='none', zorder=3 ,alpha=0.7)

xwlcfunc = np.logspace(-2,2,100)
def WLC(x,lk):
    return ((lk*lk)/2)*(2*x/lk+np.exp(-2*x/lk)-1)*10
ywlcfunc = WLC(xwlcfunc,0.3)

def RING(x):
    return np.power (x, 2./3.)
yring = RING(xwlcfunc)
ax1.plot(xwlcfunc, ywlcfunc, 'black', label= 'WLC', zorder =9)
ax1.plot(xwlcfunc[32:-12], yring[32:-12], 'r', label= 'RING', zorder =10)

snucleus_diameter = np.ones_like(xwlcfunc)
snucleus_diameter.fill(9.67)
hnucleus_diameter = np.ones_like(xwlcfunc)
hnucleus_diameter.fill(8.29)

ax1.plot(xwlcfunc,snucleus_diameter,'--', label='$R^2_g$ soft')
ax1.plot(xwlcfunc,hnucleus_diameter,'--', label='$R^2_g$ hard')


ax1.grid(True, which='both')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.legend()
left, right = ax1.get_xlim()
ax1.set_xlim(0.025, 50)
down, up = ax1.get_ylim()
ax1.set_ylim(0.025, 50)
ax1.set_xlabel('$N_{bp}(Mbp)$', size=20)
ax1.set_ylabel('$[R^2(N_{bp})][\mu m^2]$', size=20)


left, right = ax2.get_xlim()
ax2.set_xlim(50, right)
down, up = ax2.get_ylim()
ax2.set_ylim(10, up)

#fig.subplots_adjust(wspace=.3)
fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=0.3)

#print("r2n:")
#for i in sr2n.mean(axis=0):
#    print(i)
#print('error:')
#for i in np.square(sr2n).std(axis=0):
#    print(i)

fig.savefig('foo.jpg')
plt.show()
