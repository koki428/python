##浮力、ローレンツ力、トータルの力の比較
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

print("input caseid id (3 digit)")
caseid=0
caseid=input()
caseid="d"+caseid.zfill(3)

datadir="../run/"+caseid+"/data/"

d=R2D2.R2D2_data(datadir)
for key in d.p:
    exec('%s=%s%s%s' % (key, 'd.p["',key,'"]'))

try:
    n0
except NameError:
    n0=0
if  n0>d.p["nd"]:
    n0=d.p["nd"]

print("Maximum time step=",nd,"time="\
        ,dtout*float(nd)/3600./24.,"[day]")

plt.close("all")

t0=d.read_time(0,silent=True)
d.read_qq_2d(0,silent=True)

dx=(xmax-xmin)/ix
dy=(ymax-ymin)/jx
dz=(zmax-zmin)/kx

print('input max step number')
nd=input()
nd=int(nd)

##メインチューブの定義
t=d.read_time(nd,silent=True)
d.read_qq_2d(nd,silent=True)

bx=d.q2['bx']
by=d.q2['by']
bz=d.q2['bz']

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
hight=np.zeros(nd+1)
total_mean=np.zeros(nd+1)
bouyancy_mean=np.zeros(nd+1)
lorentz_mean=np.zeros(nd+1)
inertia_mean=np.zeros(nd+1)

for n in range(n0,nd+1):
    d=R2D2.R2D2_data(datadir)
    print(n)
    t=d.read_time(n,silent=True)
    d.read_qq_2d(n,silent=True)

    #gx=np.zeros(ix)

    bx=d.q2['bx']
    by=d.q2['by']
    bz=d.q2['bz']
    vx=d.q2['vx']
    vy=d.q2['vy']
    vz=d.q2['vz']
    ro=d.q2['ro']
    se=d.q2['se']
    #pr=d.q2['pr']
    #Ro0,tmp=np.meshgrid(ro0,y,indexing='ij')
    #Pr0,tmp=np.meshgrid(pr0,y,indexing='ij')
    Gx,tmp=np.meshgrid(gx,y,indexing='ij')
    DprDro,tmp=np.meshgrid(dprdro,y,indexing='ij')
    DprDse,tmp=np.meshgrid(dprdse,y,indexing='ij')
    
    #Ro=ro+Ro0
    #Pr=pr+Pr0

    #平行移動
    bx=np.roll(bx,roll,axis=1)
    by=np.roll(by,roll,axis=1)
    bz=np.roll(bz,roll,axis=1)
    vx=np.roll(vx,roll,axis=1)
    vy=np.roll(vy,roll,axis=1)
    vz=np.roll(vz,roll,axis=1)
    ro=np.roll(ro,roll,axis=1)
    se=np.roll(se,roll,axis=1)
    #pr=np.roll(pr,roll,axis=1)
    Gx=np.roll(Gx,roll,axis=1)
    DprDro=np.roll(DprDro,roll,axis=1)
    DprDse=np.roll(DprDse,roll,axis=1)

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

    cy=centroids[1,0]
    cx=centroids[1,1]

    psi=pd.DataFrame(psi)
    #高さの計算
    hight[n]=cx*dx*1.e-8
    #print('hight[',n,'] =',hight[n])

    #buoyancyの計算
    mag_bouyancy=np.zeros((ix,jx))
    pr=DprDro*ro+DprDse*se
    for i in range(2,ix-2):
        for j in range(2,jx-2):
            mag_bouyancy[i,j]=-(-pr[i+2,j]+8*pr[i+1,j]-8*pr[i-1,j]+pr[i-2,j])/(12*dx)-ro[i,j]*Gx[i,j]

    #mag_bouyancy=pd.DataFrame(mag_bouyancy)
    #bouyancy[n]=mag_bouyancy[psi > center].mean().mean()
    #print('bouyancy[',n,'] =',bouyancy[n])

    #lorentz forceの計算
    mag_lorentz=np.zeros((ix,jx))
    j_x=np.zeros((ix,jx))
    j_y=np.zeros((ix,jx))
    j_z=np.zeros((ix,jx))
    for i in range(2,ix-2):
        for j in range(2,jx-2):
            j_x[i,j]=(-bz[i,j+2]+8*bz[i,j+1]-8*bz[i,j-1]+bz[i,j-2])/(12*dy)
            j_y[i,j]=-(-bz[i+2,j]+8*bz[i+1,j]-8*bz[i-1,j]+bz[i-2,j])/(12*dx)
            j_z[i,j]=(-by[i+2,j]+8*by[i+1,j]-8*by[i-1,j]+by[i-2,j])/(12*dx)\
                    -(-bx[i,j+2]+8*bx[i,j+1]-8*bx[i,j-1]+bx[i,j-2])/(12*dy)

            mag_lorentz[i,j]=1/(4*math.pi)*(j_y[i,j]*bz[i,j]-j_z[i,j]*by[i,j])
    
    #lorentz_force=pd.DataFrame(lorentz_force)
    #lorentz[n]=lorentz_force[psi > center].mean().mean()
    #print('lorentz[',n,'] =',lorentz[n])

    #inertia forceの計算
    mag_inertia=np.zeros((ix,jx))
    for i in range(2,ix-2):
        for j in range(2,jx-2):
            mag_inertia[i,j]=-ro[i,j]*(vx[i,j]*(-vx[i+2,j]+8*vx[i+1,j]-8*vx[i+1,j]+vx[i-2,j])/(12*dx)\
                                      +vy[i,j]*(-vx[i,j+2]+8*vx[i,j+1]-8*vx[i,j-1]+vx[i,j-2])/(12*dy))

    #inertia_force=pd.DataFrame(inertia_force)
    #inertia[n]=inertia_force[psi > center].mean().mean()

    #total[n]=inertia[n]+bouyancy[n]+lorentz[n]
    #print("total[",n,"] =",total[n])

    """
    ##熱対流のみの運動方程式の各要素の計算
    d=R2D2.R2D2_data("../run/d004/data/")
    t=d.read_time(300+n,silent=True)
    d.read_qq_2d(300+n,silent=True)

    #ro0=np.zeros(ix)
    #pr0=np.zeros(ix)
    gx=np.zeros(ix)

    bx=d.q2['bx']
    by=d.q2['by']
    bz=d.q2['bz']
    vx=d.q2['vx']
    vy=d.q2['vy']
    vz=d.q2['vz']
    ro=d.q2['ro']
    pr=d.q2['pr']
    #ro0,tmp=np.meshgrid(ro0,y,indexing='ij')
    #pr0,tmp=np.meshgrid(pr0,y,indexing='ij')
    gx,tmp=np.meshgrid(gx,y,indexing='ij')

    #Ro=ro+ro0
    #Pr=pr+pr0

    #平行移動
    bx=np.roll(bx,roll,axis=1)
    by=np.roll(by,roll,axis=1)
    bz=np.roll(bz,roll,axis=1)
    vx=np.roll(vx,roll,axis=1)
    vy=np.roll(vy,roll,axis=1)
    vz=np.roll(vz,roll,axis=1)
    ro=np.roll(ro,roll,axis=1)
    pr=np.roll(pr,roll,axis=1)
    gx=np.roll(gx,roll,axis=1)

    #浮力の計算
    th_bouyancy=np.zeros((ix,jx))
    for i in range(2,ix-2):
        for j in range(2,jx-2):
            th_bouyancy[i,j]=-(-pr[i+2,j]+8*pr[i+1,j]-8*pr[i-1,j]+pr[i-2,j])/(12*dx)-ro[i,j]*gx[i,j]
    
    #ローレンツ力の計算
    th_lorentz=np.zeros((ix,jx))
    j_x=np.zeros((ix,jx))
    j_y=np.zeros((ix,jx))
    j_z=np.zeros((ix,jx))
    for i in range(2,ix-2):
        for j in range(2,jx-2):
            j_x[i,j]=(-bz[i,j+2]+8*bz[i,j+1]-8*bz[i,j-1]+bz[i,j-2])/(12*dy)
            j_y[i,j]=-(-bz[i+2,j]+8*bz[i+1,j]-8*bz[i-1,j]+bz[i-2,j])/(12*dx)
            j_z[i,j]=(-by[i+2,j]+8*by[i+1,j]-8*by[i-1,j]+by[i-2,j])/(12*dx)\
                    -(-bx[i,j+2]+8*bx[i,j+1]-8*bx[i,j-1]+bx[i,j-2])/(12*dy)

            th_lorentz[i,j]=1/(4*math.pi)*(j_y[i,j]*bz[i,j]-j_z[i,j]*by[i,j])
    
    #慣性力の計算
    th_inertia=np.zeros((ix,jx))
    for i in range(2,ix-2):
        for j in range(2,jx-2):
            th_inertia[i,j]=-ro[i,j]*(vx[i,j]*(-vx[i+2,j]+8*vx[i+1,j]-8*vx[i+1,j]+vx[i-2,j])/(12*dx)\
                                     +vy[i,j]*(-vx[i,j+2]+8*vx[i,j+1]-8*vx[i,j-1]+vx[i,j-2])/(12*dy))
    """
    th_bouyancy=np.zeros((ix,jx))
    th_lorentz=np.zeros((ix,jx))
    th_inertia=np.zeros((ix,jx))
    bouyancy=np.zeros((ix,jx))
    lorentz=np.zeros((ix,jx))
    inertia=np.zeros((ix,jx))
    total=np.zeros((ix,jx))

    bouyancy=mag_bouyancy-th_bouyancy
    lorentz=mag_lorentz-th_lorentz
    inertia=mag_inertia-th_inertia

    total=bouyancy+lorentz+inertia

    bouyancy=pd.DataFrame(bouyancy)
    lorentz=pd.DataFrame(lorentz)
    inertia=pd.DataFrame(inertia)
    total=pd.DataFrame(total)

    bouyancy_mean[n]=bouyancy[psi > center].mean().mean()
    lorentz_mean[n]=lorentz[psi > center].mean().mean()
    inertia_mean[n]=inertia[psi > center].mean().mean()
    total_mean[n]=total[psi > center].mean().mean()

    print('bouyancy =',bouyancy_mean[n])
    print('lorentz =',lorentz_mean[n])
    print('inertial =',inertia_mean[n])
    print('total =',total_mean[n])



##plot
time=np.linspace(0,nd,nd+1)
fig=plt.figure(figsize=(12,6))
plt.rcParams['mathtext.fontset']='cm'

ax1=fig.add_subplot(1,2,1)
ax1.plot(time,bouyancy_mean,label=r"$-\frac{\partial p}{\partial x}-\rho g$ (bouyancy)")
ax1.plot(time,lorentz_mean,label=r"$\frac{1}{4\pi}[(\boldsymbol{\nabla} \times \boldsymbol{B}) \times \boldsymbol{B}]_{x}$ (lorentz)")
ax1.plot(time,inertia_mean,label=r'$\rho(\boldsymbol{v} \cdot \boldsymbol{\nabla})v_{x}$ (inertial)')
ax1.plot(time,total_mean,label="total")
ax1.set_xlabel('Time [h]')
ax1.set_ylabel('Force $[g/cm^2/s^2]$')
ax1.set_ylim(-0.1,0.1)
plt.legend()

ax2=fig.add_subplot(1,2,2)
ax2.plot(time,(xmin-rsun)*1.e-8+hight,label='hight')
ax2.set_xlabel('Time [h]')
ax2.set_ylabel('Hight [Mm]')
ax2.set_ylim((xmin-rsun)*1.e-8,(xmax-rsun)*1.e-8)
#ax1.set_xlabel('x [Mm]')
#ax1.set_ylabel('vx [km/s]')




    

