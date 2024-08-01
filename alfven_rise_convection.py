#アルフベン速度、磁束管の上昇速度、熱対流速度の比較
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

print("input caseid id")
caseid=0
caseid=input()
caseid='d'+caseid.zfill(3)

datadir="../run/"+caseid+"/data/"

d=R2D2.R2D2_data(datadir)
for key in d.p:
    exec('%s=%s%s%s' % (key, 'd.p["',key,'"]'))

try:
    n0
except NameError:
    n0=0
if n0>d.p["nd"]:
    n0=d.p["nd"]

print("Maximum time step=",nd,"time="\
        ,dtout*float(nd)/3600./24.,"[day]")

plt.close("all")

t0=d.read_time(0,silent=True)
d.read_qq_2d(0,silent=True)

dx=(xmax-xmin)/ix
dy=(ymax-ymin)/jx

##メインチューブの定義
print('input max step number')
nd=input()
nd=int(nd)

print('roll number')
roll=input()
roll=int(roll)

t=d.read_time(nd-1,silent=True)
d.read_qq_2d(nd-1,silent=True)

bx=d.q2['bx']
by=d.q2['by']
bz=d.q2['bz']

#平行移動
bx=np.roll(bx,roll,axis=1)
by=np.roll(by,roll,axis=1)
bz=np.roll(bz,roll,axis=1)

#psiの計算
#psi_x=np.zeros((ix,jx))
#psi_y=np.zeros((ix,jx))
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

#2値化を行う
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

##各ステップでの計算
hight=np.zeros(nd)
velocity=np.zeros(nd)
convection_rms=np.zeros((nd,ix))
convection_mean=np.zeros(ix)
alfven=np.zeros(nd)
bz_center=np.zeros(nd)
rho=np.zeros(nd)

for n in range(n0,nd):
    d=R2D2.R2D2_data(datadir)
    for key in d.p:
        exec('%s=%s%s%s' % (key, 'd.p["',key,'"]'))

    print(n)
    t=d.read_time(n,silent=True)
    d.read_qq_2d(n,silent=True)
    
    bx=d.q2['bx']
    by=d.q2['by']
    bz=d.q2['bz']
    vx=d.q2['vx']

    #平行移動
    bx=np.roll(bx,roll,axis=1)
    by=np.roll(by,roll,axis=1)
    bz=np.roll(bz,roll,axis=1)
    vx=np.roll(vx,roll,axis=1)

    #psiの計算
    #psi_x=np.zeros((ix,jx))
    #psi_y=np.zeros((ix,jx))
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
        for m in range(1,nlabels):
            if (stats[m,4] > area_max):
                area_max=stats[m,4]
                mm=m
        stats[1,4]=stats[mm,4]
        centroids[1,:]=centroids[mm,:]

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
        
    
    #alfven速度の計算
    total_alfven=0
    count=0
    ro=d.q2['ro']
    Ro0,tmp=np.meshgrid(ro0,y,indexing='ij')
    Ro=ro+Ro0
    Ro=np.roll(Ro,roll,axis=1)
    #"""
    for j in range(0,jx):
        for i in range(0,ix):
            if (psi[i,j] > center):
                total_alfven=total_alfven+bz[i,j]/np.sqrt(4*math.pi*Ro[i,j])
                count=count+1
    alfven[n]=total_alfven/count*1.e-5
    #"""
    #alfven[n]=bz[int(cx),int(cy)]/np.sqrt(4*math.pi*(Ro0[int(cx),int(cy)]+ro[int(cx),int(cy)]))*1.e-5
    print('alfven[',n,'] =',alfven[n])

    """
    bz_center[n]=bz[int(cx),int(cy)]
    rho[n]=Ro0[int(cx),int(cy)]+ro[int(cx),int(cy)]
    """
    #熱対流の速度の計算
    
    d=R2D2.R2D2_data("../run/d004/data/")
    t=d.read_time(300+n,silent=True)
    d.read_qq_2d(300+n,silent=True)

    vx=d.q2['vx']
    convection=0
    count=0
    """
    for j in range(0,jx):
        for i in range(0,ix):
            if (vx[i,j] > 0):
                convection=convection+vx[i,j]*vx[i,j]
                count=count+1
    convection_rms[n]=np.sqrt(convection/count)*1.e-5
    """
    #convection_rms[n]=vx[int(cx),int(cy)]*1.e-5

    for i in range(0,ix):
        for j in range(0,jx):
            if (vx[i,j] > 0):
                convection=convection+vx[i,j]*vx[i,j]
                count=count+1
        convection_rms[n,i]=np.sqrt(convection/count)
        convection=0
        count=0
    #print('convection_rms[',n,'] =',convection_rms[n])
convection_mean=np.mean(convection_rms,axis=0)*1.e-5

##plot
#time=np.linspace(0,n*8,n+1)
fig=plt.figure(figsize=(8,8))
ax1=fig.add_subplot(1,1,1)
ax1.plot(((xmin-rsun)*1.e-8)+hight,velocity,marker='D',label="rise velocity")
ax1.plot(((xmin-rsun)*1.e-8)+hight,alfven,marker='D',label="alfven velocity")
ax1.plot((x-rsun)*1.e-8,convection_mean,label="convection velocity")
plt.legend()
ax1.set_xlabel('x [Mm]')
ax1.set_ylabel('vx [km/s]')
"""
fig=plt.figure(figsize=(8,8))
ax2=fig.add_subplot(2,1,1)
ax3=fig.add_subplot(2,1,2)
ax2.plot(((xmin-rsun)*1.e-8)+hight,bz_center)
ax3.plot(((xmin-rsun)*1.e-8)+hight,rho)
"""    

