#メインチューブの半径の2倍の内側での速度の平均

import numpy as np
import math
import sympy as sym
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import R2D2
import cv2
import sys
import os
import pandas as pd

print('roll number')
roll=input()
roll=int(roll)

print("input caseid id")
caseid=0
caseid=input()
caseid='d'+caseid.zfill(3)

datadir="../run/"+caseid+"/data/"

d=R2D2.R2D2_data(datadir)
for key in d.p:
    exec('%s=%s%s%s' % (key, 'd.p["',key,'"]'))

print('input max step number')
nd=input()
nd=int(nd)

plt.close("all")

t0=d.read_time(0,silent=True)
d.read_qq_2d(0,silent=True)

dx=(xmax-xmin)/ix
dy=(ymax-ymin)/jx

t=d.read_time(nd,silent=True)
d.read_qq_2d(nd,silent=True)

bx=d.q2['bx']
by=d.q2['by']

#平行移動
bx=np.roll(bx,roll,axis=1)
by=np.roll(by,roll,axis=1)

#psiの計算
psi=np.zeros((ix,jx))
psi[0,0]=0.0
i=0
for j in range(1,jx-1):
    psi[i,j]=psi[i,j-1]-(bx[i,j-1]*dy+4*bx[i,j]*dy+bx[i,j+1]*dy)/3
for i in range(1,ix):
    j=0
    psi[i,j]=psi[i-1,j]+by[i,j]*dx
    for j in range(1,jx-1):
        psi[i,j]=psi[i,j-1]-(bx[i,j-1]*dy+4*bx[i,j]*dy+bx[i,j+1]*dy)/3

# 2値化を行う
psi_thr=np.zeros((ix,jx))
psi_max=np.max(psi)
psi_min=100.0
while True:
    center=(psi_max+psi_min)*0.5
    for j in range(0,jx):
        for i in range(0,ix):
            if (psi[i,j] > center):
                psi_thr[i,j]=110.0
            else:
                psi_thr[i,j]=0.0

    ret,binary_psi = cv2.threshold(psi_thr, 100, 255, cv2.THRESH_BINARY)

    psi_label=binary_psi.astype('uint8')
    nlabels,labellmages,stats,centroids=cv2.connectedComponentsWithStats(psi_label)
    if (nlabels-1 > 1):
        psi_min=center
    else:
        psi_max=center
    if (abs(psi_max-psi_min) < 0.5) and (nlabels-1==1):
        break

n0=0
hight=np.zeros(nd)
velocity=np.zeros(nd)
#convection_around=np.zeros(nd)
for m in range(1,5):
    globals()['convection_around'+str(m)]=np.zeros(nd)
for n in range(n0,nd):
    d=R2D2.R2D2_data(datadir)
    print(n)
    t=d.read_time(n,silent=True)
    d.read_qq_2d(n,silent=True)

    bx=d.q2['bx']
    by=d.q2['by']
    vx=d.q2['vx']

    #平行移動
    bx=np.roll(bx,roll,axis=1)
    by=np.roll(by,roll,axis=1)
    vx=np.roll(vx,roll,axis=1)

    #psiの計算
    psi=np.zeros((ix,jx))
    psi[0,0]=0.0
    i=0
    for j in range(1,jx-1):
        psi[i,j]=psi[i,j-1]-(bx[i,j-1]*dy+4*bx[i,j]*dy+bx[i,j+1]*dy)/3
    for i in range(1,ix):
        j=0
        psi[i,j]=psi[i-1,j]+by[i,j]*dx
        for j in range(1,jx-1):
            psi[i,j]=psi[i,j-1]-(bx[i,j-1]*dy+4*bx[i,j]*dy+bx[i,j+1]*dy)/3
    
    # 2値化を行う
    psi_thr=np.zeros((ix,jx))
    for j in range(0,jx):
        for i in range(0,ix):
            if (psi[i,j] > center):
                psi_thr[i,j]=110.0
            else:
                psi_thr[i,j]=0.0

    ret,binary_psi = cv2.threshold(psi_thr, 100, 255, cv2.THRESH_BINARY)

    psi_label=binary_psi.astype('uint8')
    nlabels,labellmages,stats,centroids=cv2.connectedComponentsWithStats(psi_label)

    if (nlabels-1 > 1):
        area_max=0
        for k in range(1,nlabels):
            if (stats[k,4] > area_max):
                area_max=stats[k,4]
                kk=k
        stats[1,4]=stats[kk,4]
        centroids[1,:]=centroids[kk,:]
    
    #メインチューブの有効半径
    r_eff=0
    r_eff=np.sqrt(stats[1,4]*dx*dy/math.pi)
    #メインチューブの中心
    cy=centroids[1,0]
    cx=centroids[1,1]

    #高さの計算
    hight[n]=cx*dx*1.e-8
    print('hight[',n,'] =',hight[n])
    #速度の計算
    if (n != 0):
        velocity[n]=(hight[n]-hight[n-1])/(8*60*60)*1.e+3
        #velocity[n]=vx[int(cx),int(cy)]*1.e-5
        print('velocity[',n,'] =',velocity[n])
    
    #半径n倍内の熱対流速度の計算
    d=R2D2.R2D2_data("../run/d004/data/")
    t=d.read_time(300+n,silent=True)
    d.read_qq_2d(300+n,silent=True)
    convection_total=0
    count=0
    rr=np.zeros((ix,jx))
    for i in range(0,ix):
        for j in range(0,jx):
            rr[i,j]=np.sqrt((cx*dx-i*dx)**2+(cy*dy-j*dy)**2)
    
    vx=d.q2['vx']
    vx=np.roll(vx,roll,axis=1)

    vx=pd.DataFrame(vx)
    rr=pd.DataFrame(rr)
    for m in range(1,5):
        globals()['convection_around'+str(m)][n]=vx[rr < m*r_eff].mean().mean()*1.e-5
        #print('convection[',n,'] =',convection_around[n])
        #globals()['convection_around'+str(m)][n]=convection_around

##熱対流の計算
"""
convection_rms=np.zeros((nd,ix))
convection_mean=np.zeros(ix)
d=R2D2.R2D2_data("../run/d004/data/")
for n in range(0,nd):
    t=d.read_time(300+n,silent=True)
    d.read_qq_2d(300+n,silent=True)
    vx=d.q2['vx']
    convection=0
    count=0
    for i in range(0,ix):
        for j in range(0,jx):
            if (vx[i,j] > 0):
                convection=convection+vx[i,j]*vx[i,j]
                count=count+1
        convection_rms[n,i]=np.sqrt(convection/count)
        convection=0
        count=0
convection_mean=np.mean(convection_rms,axis=0)*1.e-5
"""
##plot
time=np.linspace(0,nd*8,nd)
fig=plt.figure(figsize=(16,8))
ax1=fig.add_subplot(1,2,1)
#ax1.plot(((xmin-rsun)*1.e-8)+hight,velocity,label="rise velocity")
#ax1.plot(((xmin-rsun)*1.e-8)+hight,convection_around,label="around velocity")
#ax1.plot((x-rsun)*1.e-8,convection_mean,label="convection velocity")
ax1.plot(time,velocity,label="rise velocity")
ax1.plot(time,convection_around1,label="around velocity(r)")
ax1.plot(time,convection_around2,label="around velocity(2r)")
ax1.plot(time,convection_around3,label="around velocity(3r)")
ax1.plot(time,convection_around4,label="around velocity(4r)")
#ax1.plot((x-rsun)*1.e-8,convection_mean,label="convection velocity")
ax1.hlines(0,0,nd*8,colors='black',linestyle='dashed')
ax1.set_xlim(0,nd*8)
ax1.set_xlabel('time [h]')
ax1.set_ylabel('vx [km/s]')
plt.legend()

ax2=fig.add_subplot(1,2,2)
ax2.plot(time,((xmin-rsun)*1.e-8)+hight,label='hight')
ax2.set_xlim(0,nd*8)
ax2.set_ylim((xmin-rsun)*1.e-8,(xmax-rsun)*1.e-8)
ax2.set_xlabel('time [h]')
ax2.set_ylabel('hight [Mm]')
plt.legend()





