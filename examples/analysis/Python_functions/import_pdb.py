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
        data             = np.load(filename[:-3]+'npy')
        time             = np.load(filename[:-4]+'_t.npy')
        energy           = np.load(filename[:-4]+'_en.npy')
        
        print('{} imported successfully from the "npy"s.'.format(name))
        return data, time, energy
        
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
                en =[]
            elif words[0] == 'REMARK':
                time = np.append(time, float(words[2].replace('time=','')) )
                en = np.append(en, float(words[4].replace('energy=','')) )
                en = np.append(en, float(words[6].replace('energy=','')) )
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
                    data   = np.asarray([xyz])
                    energy = np.asarray([en])
                else:
                    data   = np.concatenate((data,[xyz]))
                    energy = np.concatenate((energy,[en]))
    np.save(filename[:-4],data)
    np.save(filename[:-4]+'_t',time)
    np.save(filename[:-4]+'_en',energy)
    
    print('{} imported successfully.'.format(name))
    return data, time, energy



def import_pdb_with_biopandas(fname, label = None):
    from biopandas.pdb import PandasPdb
    cite = '''
    @article{raschkas2017biopandas,
             doi = {10.21105/joss.00279},
             url = {http://dx.doi.org/10.21105/joss.00279},
             year  = {2017},
             month = {jun},
             publisher = {The Open Journal},
             volume = {2},
             number = {14},
             author = {Sebastian Raschka},
             title = {BioPandas: Working with molecular structures in pandas DataFrames},
             journal = {The Journal of Open Source Software}
             }
    '''
    print(cite)
    #import numpy as np
    #from biopandas.pdb import PandasPdb
    ppdb = PandasPdb()
    ppdb.read_pdb(fname)
    properties = ppdb.df['ATOM'].dtypes.index
    
    num_of_atoms = ppdb.df['ATOM'][properties[1]].max()
    frame_nums   = ppdb.df['ATOM'][properties[1]].size//num_of_atoms
    coords = ppdb.df['ATOM'][properties[11:14]].to_numpy().reshape((frame_nums,num_of_atoms,3))
    
    
    properties = ppdb.df['OTHERS'].dtypes.index
    rows = str(ppdb.df['OTHERS'][properties[1]][1::3].values).split()
    
    time   = np.asarray([float(x.replace('time=',''  ) ) for x in rows[2::8]])
    energy = np.asarray([float(x.replace('energy=','') ) for x in rows[4::8]])
    return coords, time, energy



def read_frame(datamap):
        import sys
        global linehist
        try:
            firstcoords = np.zeros(3)
            
            line  = datamap.readline()
            words = line.split()
            firstcoords[0]=float(line[30:38])
            firstcoords[1]=float(line[38:46])
            firstcoords[2]=float(line[46:54])
            
            linehist = [line,line]
            line  = datamap.readline()
            words = line.split()
            while words[0].decode('utf-8') == 'ATOM':
                linehist[0] = linehist[1]
                linehist[1] = line
                firstcoords = np.append(firstcoords, [float(line[30:38]),float(line[38:46]),float(line[46:54]) ]) 
                line   = datamap.readline()
                words  = line.split()
                
            
            
            return firstcoords
            
        except:
            print('reached EOF')
            raise

def import_pdb_with_mmap(fname, label = None):
    import sys
    import mmap
    
    energies = np.zeros(2)
    time     = np.zeros(1)
    
    try:
        coords   = np.load(fname[:-4]+'_coords.npy')
        time     = np.load(fname[:-4]+'_time.npy')
        energies = np.load(fname[:-4]+'_ens.npy')
        
        print('Loaded from .npy files.')
        
        return coords, time, energies
    except:
        pass
    
    
    try:
        f=open(fname)
        num_lines = sum(1 for line in f)
        f.close()
        print('pdb opend successfully. It contains {} lines.'.format(num_lines))
    except:
        print('Cannot find the pdb.')
        
    
    
    with open(fname, "r+") as f:
        try:
            datamap = mmap.mmap(f.fileno(),0, offset=0, prot=mmap.PROT_READ)
            
        except:
            print("Unexpected error:", sys.exc_info()[0])
            raise
        
        line  = datamap.readline()
        line  = datamap.readline()
        words = line.split()
            
        time[0]     = float(words[2][5:])
        energies[0] = float(words[4][7:])
        energies[1] = float(words[6][7:])
            
        tempcoords = read_frame(datamap)
        
        num_particles = len(tempcoords)//3
        num_frames = num_lines//(num_particles+3)
        
        coords = np.zeros((num_frames, num_particles, 3))
        coords[0,:,:] = tempcoords.reshape((num_particles,3))
        print('pdb contains {} frames with each frame displaying {} particles'.format(num_frames,num_particles))
            
        
        line  = datamap.readline()
        line  = datamap.readline()
        words = line.split()
            
        
        
        
        counter = 0
        while True:

            try:
                time     = np.append(time,      float(words[2][5:]))
                energies = np.append(energies, [float(words[4][7:]), float(words[6][7:])])
                counter+=1
                temp = read_frame(datamap)
                line  = datamap.readline()
                line  = datamap.readline()
                words = line.split()
            except:
                break
            
            coords[counter,:,:] = temp.reshape((num_particles,3))
            
            
        
        energies = energies.reshape((num_frames,2))
        print('last row of coordinates read:\n{}'.format(coords[-1,-1,-3:]))
        print('finished')
        
        return coords, time, energies

