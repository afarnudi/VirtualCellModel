# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np

def autocorr(x):
    cor = np.correlate(x-np.mean(x),x-np.mean(x), 'same')
    cor2 = cor[len(cor)//2:]
    return cor2/cor2[0]

def autocorr2(x):
    m = np.mean(x)
    N = x.size//2
    results = np.zeros(N)
    results[0] = 1
    for tau in range(1,N):
        for i, j in zip(x[:-tau],x[tau:]):
            results[tau] += (i-m)*(j-m)
        results[tau] = results[tau]/(np.sum(np.square(x)) -2*m*np.sum(x) - m*m )
    return results
 
    
def import_data(filename):
    time=[]
    data=[]
    d = []
    with open(filename, "r") as f:
        for line in f:
            words = line.split()
            if words[0]=='time:':
                time = np.append(time, float(words[1]))
                if len(time) == 2:
                    data = np.asarray([d])
                    d=[]
                elif len(time)>2:
                    data = np.concatenate((data, [d]))
                    d=[]
            else:
                temp=[]
                for i in  words[1:]:
                    temp = np.append(temp,float(i))
                if len(d)==0:
                    d = np.asarray([temp])
                else:
                    d = np.concatenate((d,[temp] ))
    data = np.concatenate((data,[d]))
    return data, time

def read_file(filename):
    try:
        data = np.load(filename[:-3]+'npy')
        time = np.load(filename[:-4]+'_t.npy')
        print('{} imported successfully from the "npy"s.'.format(filename))
        return data, time
    except:
        print('No npy files located. importing from the txt.')
        
        data, time = import_data(filename)
        np.save(filename[:-4],data)
        np.save(filename[:-4]+'_t',time)
        print('{} imported successfully.'.format(filename))
        return data, time