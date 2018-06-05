#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  5 07:42:29 2018
this is to calculate the auto-correlation function of single rolling atoms
@author: jiahaoz
"""
import numpy as np
import matplotlib.pyplot as plt
def autocorri(data,step):
    le=len(data);
    nextstart=step;
    averagetimes=le-nextstart;
    sum=0;
    for i in range(averagetimes):
        sum=sum+np.dot(data[i,0:3],data[i+step,0:3]);
    return sum/averagetimes;
atomlist=np.loadtxt("cadata.txt");
ensemble=len(atomlist);
timestep=0.005;
yall=[];
for j in atomlist:
    file="move"+str(int(j))+".txt";
    data=np.loadtxt(file,delimiter=' ');
    le=len(data);
    print j
    x=range(1,le/2);
    y=[];
    for i in x:
        y.append(autocorri(data,i));
    yall.append(y);
y=np.sum(np.array(yall),axis=0);
y=[i/ensemble for i in y];
x=[i*timestep for i in x];
plt.plot(x,y);
plt.xlabel("time/ps");
ylabel("auto-corrilation(A^2)")