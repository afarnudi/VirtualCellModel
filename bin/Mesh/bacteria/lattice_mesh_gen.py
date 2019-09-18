#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 10:40:44 2019

@author: alifarnudi
"""
import numpy as np
#import matplotlib.pyplot as plt


s=2
a=np.sqrt(3)
d=a*np.sqrt(3)/2
N=20
n=int(N/2)
l=n*a
w=n*d

line=1
index=1
with open("membrane_lattice.geo","w") as file:
    for i in range(N):
        if line==1:
            file.write("\n")
            line=0
        else:
            line=1
            file.write("\n")
        for j in range(N+1):
            if line==0:
                file.write("Point({})={{ {:2.2f},0,{:2.2f},{} }};\n".format(index,-l+j*a,-w+i*d,s))
                index+=1
            else:
                file.write("Point({})={{ {:2.2f},0,{:2.2f},{} }};\n".format(index,-l+j*a+0.5*a,-w+i*d,s))
                index+=1
    file.write("\n")
    index=1
    line=0
    for i in range(1,N):
        if line==1:
            line=0
        else:
            line=1
        for j in range(1,N+2):
            file.write("Line({})={{ {:2.0f},{:2.0f} }};\n".format(index,(i-1)*(N+1)+j,i*(N+1)+j))
            index+=1
            if j<N+1:
                file.write("Line({})={{ {:2.0f},{:2.0f} }};\n".format(index,(i-1)*(N+1)+j,(i-1)*(N+1)+j+1))
                index+=1
                if line==1:
                    file.write("Line({})={{ {:2.0f},{:2.0f} }};\n".format(index,(i-1)*(N+1)+j+1,i*(N+1)+j))
                    index+=1
                else:
                    file.write("Line({})={{ {:2.0f},{:2.0f} }};\n".format(index,(i-1)*(N+1)+j,i*(N+1)+j+1))
                    index+=1
    for j in range(N):
        file.write("Line({})={{ {:2.0f},{:2.0f} }};\n".format(index,N*N+j,N*N+j+1))
        index+=1
    line=-1
    
    b=3
    a=1
    c=0
    d=(N-1)*(3*N+1)+1
    index=1
    for i in range(1,N):
        if line==1:
            line=-1
        else:
            line=1
        for j in range(N):
            file.write("Curve Loop({})={{ {:2.0f},{:2.0f},{:2.0f} }};\n".format(index,-a,a+1,line*(a+2)))
            file.write("Plane Surface({})={{{:2.0f}}};\n".format(index,index))
            file.write("Physical Surface({})={{{:2.0f}}};\n".format(index,index))
            a+=3
            index+=1
            if i <N-1:
                file.write("Curve Loop({})={{ {:2.0f},{:2.0f},{:2.0f} }};\n".format(index,-b,b+int(-line*0.5+1.5),-(b-line+(N*3+1))))
                file.write("Plane Surface({})={{{:2.0f}}};\n".format(index,index))
                file.write("Physical Surface({})={{{:2.0f}}};\n".format(index,index))
                b+=3
                index+=1
            else:
                file.write("Curve Loop({})={{ {:2.0f},{:2.0f},{:2.0f} }};\n".format(index,-((N-2)*(3*N+1)+3+c),(N-2)*(3*N+1)+3+c+1,-d))
                file.write("Plane Surface({})={{{:2.0f}}};\n".format(index,index))
                file.write("Physical Surface({})={{{:2.0f}}};\n".format(index,index))
                c+=3
                d+=1
                index+=1
        if line==1:
            a+=2
            b-=1
        else:
            b+=3
    
    
    
    
    
    
    
    
    
