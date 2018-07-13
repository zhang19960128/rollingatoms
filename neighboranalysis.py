#!usr/bin/env python2
import numpy as np
import math
def wash(pointone,pointtwo,pb):
    dis=pb[1]-pb[0];
    if(pointtwo-pointone>dis/2):
        return pointtwo-dis;
    elif(pointtwo-pointone<-1*dis/2):
        return pointtwo+dis;
    else:
        return pointtwo+0.0;
def index3Dto1D(index,p):
    re=(index[0]+p)%p+p*((index[1]+p)%p)+p**2*((index[2]+p)%p)
    return int(re)
def searchAone(index,p):
    '''index3D=[nx,ny,nz]'''
    nz=int(math.floor(index/(p**2)));
    ny=int(math.floor((index-nz*p**2)/p));
    nx=int(math.floor(index-nz*p**2-ny*p));
    '''first layer'''
    layone=[];
    for i in range(-1,2,1):
        for j in range(-1,2,1):
            for k in range(-1,2,1):
                if(math.fabs(i)+math.fabs(j)+math.fabs(k)==1):
                    layone.append(index3Dto1D([nx+i,ny+j,nz+k],p));

    return layone;
def searchAtwo(index,p):
    '''second layer'''
    nz=int(math.floor(index/(p**2)));
    ny=int(math.floor((index-nz*p**2)/p));
    nx=int(math.floor(index-nz*p**2-ny*p));
    laytwo=[];
    for i in range(-1,2,1):
        for j in range(-1,2,1):
            for k in range(-1,2,1):
                if(math.fabs(i)+math.fabs(j)+math.fabs(k)==2):
                    laytwo.append(index3Dto1D([nx+i,ny+j,nz+k],p));
    return laytwo;
def searchAthree(index,p):
    '''third layer '''
    nz=int(math.floor(index/(p**2)));
    ny=int(math.floor((index-nz*p**2)/p));
    nx=int(math.floor(index-nz*p**2-ny*p));
    laythree=[];
    for i in range(-1,2,1):
        for j in range(-1,2,1):
            for k in range(-1,2,1):
                if(math.fabs(i)+math.fabs(j)+math.fabs(k)==3):
                    laythree.append(index3Dto1D([nx+i,ny+j,nz+k],p));
    return laythree;
def varofpoint(index,raw,length):
    linenum=9+index;
    linepx=5;
    linepy=6;
    linepz=7;
    dataall=[];
    px=[];
    py=[];
    pz=[];
    for i in range(linenum,length,5009):
        line=raw[i];
        line=line.split();
        for j in range(3):
            line[j]=float(line[j]);
        dataall.append(line);
    for i in range(linepx,length,5009):
        line=raw[i];
        line=line.split();
        for j in range(2):
            line[j]=float(line[j]);
        px.append(line);
    for i in range(linepy,length,5009):
        line=raw[i];
        line=line.split();
        for j in range(2):
            line[j]=float(line[j]);
        py.append(line);
    for i in range(linepz,length,5009):
        line=raw[i];
        line=line.split();
        for j in range(2):
            line[j]=float(line[j]);
        pz.append(line);
    for i in range(len(dataall)):
        dataall[i][0]=wash(dataall[0][0],dataall[i][0],px[i]);
        dataall[i][1]=wash(dataall[0][1],dataall[i][1],py[i]);
        dataall[i][2]=wash(dataall[0][2],dataall[i][2],pz[i]);
    dataall=np.array(dataall);
    return [np.var(dataall[0:-1,0]),np.var(dataall[0:-1,1]),np.var(dataall[0:-1,2])]
calist=open("cadata.txt","r");
ca=[];
for i in calist.readlines():
    ca.append(int(i));
for i in ca:
    nei=searchAone(i,10);
    common=set(calist)-(set(calist)-set(nei));
    print common
