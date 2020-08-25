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

def autocorr3(x):
    N = x.size//2
    results = np.zeros(N)
    results[0] = 1
    for tau in range(1,N):
        results[tau] = np.sum(x[tau:]*x[:-tau])/len(x[:-tau])-(np.sum(x[tau:])/len(x[:-tau]))*(np.sum(x[:-tau])/len(x[:-tau]))
        results[tau] = results[tau]/(np.sqrt(np.var(x[:-tau])*np.var(x[tau:])))
    return results    

    
def import_data(filename):
    time=[]
    data=[]
    d = []
    
    properties = []
    props = []
    with open(filename, "r") as f:
        for line in f:
            words = line.split()
            if words[0]=='time:':
                time = np.append(time, float(words[1]))
                props = np.append(props, float(words[8].replace('Volume=','')) )
                props = np.append(props, float(words[10].replace('Area=','')) )

                if len(time) == 2:
                    data = np.asarray([d])
                    d=[]
                    
                    properties = np.asarray([props[:2]])
                    properties = np.concatenate((properties,[props[2:]] ))
                    props = []
                    
                    
                elif len(time)>2:
                    data = np.concatenate((data, [d]))
                    d=[]
                    
                    properties = np.concatenate((properties,[props] ))
                    props = []
                    
            else:
                temp=[]
                for i in  words[1:]:
                    temp = np.append(temp,float(i))
                if len(d)==0:
                    d = np.asarray([temp])
                else:
                    d = np.concatenate((d,[temp] ))
    data = np.concatenate((data,[d]))
    return data, time, properties

def read_file(filename):
    try:
        data       = np.load(filename[:-3]+'npy')
        time       = np.load(filename[:-4]+'_t.npy')
        properties = np.load(filename[:-4]+'_props.npy')
        print('{} imported successfully from the "npy"s.'.format(filename))
        return data, time, properties
    except:
        print('No npy files located. importing from the txt.')
        
        data, time, properties = import_data(filename)
        np.save(filename[:-4],data)
        np.save(filename[:-4]+'_t',time)
        np.save(filename[:-4]+'_props',properties)
        print('{} imported successfully.'.format(filename))
        return data, time, properties

def import_vels_with_mmap(fname, label = None, buffer=1):
    import sys
    import mmap
    
    time       = np.zeros(1)
    properties = np.zeros(3)
    
    try:
        data       = np.load(fname[:-16]+'_vfs.npy')
        time       = np.load(fname[:-16]+'_time.npy')
        properties = np.load(fname[:-16]+'_props.npy')
        
        print('Loaded from .npy files.')
        
        return data, time, properties
    except:
        pass
    
    try:
        f=open(fname)
        num_lines = sum(1 for line in f)
        f.close()
        print('.txt file opend successfully. It contains {} lines.'.format(num_lines))
    except:
        print('Cannot find the txt file.')
        
    
    
    with open(fname, "r+") as f:
        try:
            datamap = mmap.mmap(f.fileno(),0, offset=0, prot=mmap.PROT_READ)
            
        except:
            print("Unexpected error:", sys.exc_info()[0])
            raise
        
        line  = datamap.readline()
        words = line.split()
            
        time[0]       = float(words[1])
        # print(words[8][7:])
        # print(words[10][5:])
        properties[0] = float(words[8][7:])
        properties[1] = float(words[10][5:])
        try:
            properties[2] = float(words[12][15:])
        except:
            pass
        import sys
        try:
            
            line  = datamap.readline()
            words = line.split()
            data = np.asarray([float(x) for x in words[1:7]])
            num_data  =  6 #len(words)-1
            
            line  = datamap.readline()
            words = line.split()
            while words[0].decode('utf-8') != 'time:':
                
                data = np.append(data, np.asarray([float(x) for x in words[1:7]]) )
                line   = datamap.readline()
                words  = line.split()
                
            
            
            rows = data.size//num_data
            tempdata=data.reshape(rows,num_data)
            
        except:
            print('reached EOF')
            rows = data.size//num_data
            tempdata=data.reshape(rows,num_data)
            raise
        
        num_particles = tempdata.shape[0]
        num_frames = num_lines//(num_particles+1)
        
        data = np.zeros((num_frames, tempdata.shape[0], tempdata.shape[1]))
        data[0,:,:] = tempdata
        print(line)
        print('txt contains {} frames for {} particles'.format(num_frames,num_particles))
            
          
        
        counter = 0
        while True:
            counter+=1
            time = np.append(time, float(words[1]) )
            properties = np.append(properties,float(words[8][7:]))
            properties = np.append(properties,float(words[10][5:]))
            try:
                properties = np.append(properties,float(words[12][15:]))
            except:
                properties = np.append(properties,0.0)
            try:
                line  = datamap.readline()
                words = line.split()
                tempdata = np.asarray([float(x) for x in words[1:7]])
            
                line  = datamap.readline()
                words = line.split()
                while words[0].decode('utf-8') != 'time:':
                
                    tempdata = np.append(tempdata, np.asarray([float(x) for x in words[1:7]]) )
                    line   = datamap.readline()
                    words  = line.split()
                
                data[counter,:,:] = tempdata.reshape(rows,num_data)
            
            except:
                data[counter,:,:] = tempdata.reshape(rows,num_data)
                break
        
        print('last row of coordinates read:\n{}'.format(data[-1,-1,:]))
        print('finished')
        return data, time, properties.reshape(num_frames,3)

