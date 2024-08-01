#磁気エネルギーと運動エネルギーに高さ分布

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import math
import R2D2
import sys
import os
import cv2
import pandas as pd
import sympy as sym

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

print('Maximum time step =',nd,'time ='\
        ,dtout*float(nd)/3600./24.,'[day]')

plt.close("all")

t0=d.read_time(0,silent=True)
d.read_qq_2d(0,silent=True)

dx=(xmax-xmin)/ix
dy=(ymax-ymin)/jx

n0=0

print('input max step number')
nd=input()
nd=int(nd)

##メインチューブの定義
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
xp=np.zeros(nd)
bb=np.zeros(nd)
p_center=np.zeros(nd)
B_energy=np.zeros(nd)
Convection_energy=np.zeros(nd)

n0=0
for n in range(n0,nd):
    print(n)
    d=R2D2.R2D2_data(datadir)
    t=d.read_time(n,silent=True)
    d.read_qq_2d(n,silent=True)

    bx=d.q2['bx']
    by=d.q2['by']
    bz=d.q2['bz']
    pr=d.q2['pr']

    Pr0,tmp=np.meshgrid(pr0,y,indexing='ij')
    pre=Pr0+pr

    #平行移動
    bx=np.roll(bx,roll,axis=1)
    by=np.roll(by,roll,axis=1)
    bz=np.roll(bz,roll,axis=1)
    pre=np.roll(pre,roll,axis=1)

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
        for m in range(1,nlabels):
            if (stats[m,4] > area_max):
                area_max=stats[m,4]
                mm=m
        stats[1,4]=stats[mm,4]
        centroids[1,:]=centroids[mm,:]

    cy=centroids[1,0]
    cx=centroids[1,1]

    #cx=int(cx)
    #cy=int(cy)

    #idx=np.unravel_index(np.argmax(bz),bz.shape)
    #cx=int(idx[0])
    #cy=int(idx[1])

    #有効半径の計算
    r_eff=np.sqrt(stats[1,4]*dx*dy/math.pi)

    #メインチューブの中心から各点までの距離の計算
    rr=np.zeros((ix,jx))
    for i in range(0,ix):
        for j in range(0,jx):
            rr[i,j]=np.sqrt((cx*dx-i*dx)**2+(cy*dy-j*dy)**2)

    #高さの計算
    hight[n]=cx*dx*1.e-8

    p_center[n]=pre[int(cx),int(cy)]
    xp[n]=p_center[n]/p_center[0]

    bb[n]=np.sqrt(bx[int(cx),int(cy)]**2+by[int(cx),int(cy)]**2+bz[int(cx),int(cy)]**2)

    #メインチューブ内の磁気エネルギーの計算
    b_energy=(bx**2+by**2+bz**2)/(8*math.pi)
    b_energy=pd.DataFrame(b_energy)
    rr=pd.DataFrame(rr)

    B_energy[n]=b_energy[rr < r_eff].mean().mean()

    #メインチューブ内の熱対流の運動エネルギーの計算
    d=R2D2.R2D2_data('../run/d004/data/')
    t=d.read_time(300+n,silent=True)
    d.read_qq_2d(300+n,silent=True)

    vx=d.q2['vx']
    vy=d.q2['vy']
    ro=d.q2['ro']
    Ro0,tmp=np.meshgrid(ro0,y,indexing='ij')
    Ro=ro+Ro0

    #平行移動
    vx=np.roll(vx,roll,axis=1)
    vy=np.roll(vy,roll,axis=1)
    Ro=np.roll(Ro,roll,axis=1)
    

    convection_energy=Ro*(vx**2+vy**2)/2
    convection_energy=pd.DataFrame(convection_energy)
    Convection_energy[n]=convection_energy[rr < r_eff].mean().mean()




xp_re=list(reversed(xp))


##熱対流の運動エネルギーの計算
#kinetic=np.zeros(ix)
kinetic=np.zeros(nd)
kinetic_mean=np.zeros(ix)
kinetic_same_hight=np.zeros((nd,ix))
d=R2D2.R2D2_data("../run/d004/data/")
t=d.read_time(300+n,silent=True)
d.read_qq_2d(300+n,silent=True)

for n in range(0,nd):
    t=d.read_time(300+n,silent=True)
    d.read_qq_2d(300+n,silent=True)

    vx=d.q2['vx']
    vy=d.q2['vy']
    ro=d.q2['ro']
    Ro0,tmp=np.meshgrid(ro0,y,indexing='ij')
    Ro=ro+Ro0

    #平行移動
    vx=np.roll(vx,roll,axis=1)
    vy=np.roll(vy,roll,axis=1)
    Ro=np.roll(Ro,roll,axis=1)

    #vv=vx[cx,cy]**2+vy[cx,cy]**2
    #kinetic[n]=Ro[cx,cy]*vv/2
    #"""
    kinetic=0
    count=0
    kinetic=np.mean(Ro*(vx**2+vy**2)/2,axis=1)
    if (n==0):
        kinetic_mean=kinetic
    else:
        kinetic_mean=(kinetic_mean+kinetic)/2
    #"""
    #kinetic=0
"""
    for i in range(0,ix):
        for j in range(0,jx):
            kinetic=kinetic+Ro[i,j]*vv[i,j]/2
        kinetic_same_hight[n,i]=kinetic/jx
        kinetic=0
kinetic_mean=np.mean(kinetic_same_hight,axis=0)
"""
"""
convection_rms=np.zeros((nd,ix))
convection_mean=np.zeros(ix)
d=R2D2.R2D2_data("../run/d004/data/")
for n in range(0,nd):
    t=d.read_time(300+n,silent=True)
    d.read_qq_2d(300+n,silent=True)
    vx=d.q2['vx']
    vy=d.q2['vy']
    convection=0
    count=0
    for i in range(0,ix):
        for j in range(0,jx):
            convection=convection+vx[i,j]
            count=count+1
        convection_rms[n,i]=convection/count
        convection=0
        count=0
convection_mean=np.mean(convection_rms,axis=0)
"""
gam=5/3
##plot
time=np.linspace(0,nd*8,n+1)
fig=plt.figure(figsize=(8,8))

ax1=fig.add_subplot(1,1,1)
#plt.yscale('log')
ax1.plot(((xmin-rsun)*1.e-8)+hight,(bb[0]*xp**(1/gam))**2/(8*math.pi),label="magnetic")
#ax1.plot(((xmin-rsun)*1.e-8)+hight,alfven,marker='D',label="alfven velocity")
ax1.plot((x-rsun)*1.e-8,kinetic_mean,label="convection")
#ax1.plot(((xmin-rsun)*1.e-8)+hight,kinetic,label="kinetic")
ax1.set_yscale('log')
ax1.set_xlabel('x [Mm]')
ax1.set_ylabel('energy $[g/cm/s^2]$')
plt.legend()

"""
ax2=fig.add_subplot(1,2,2)
#ax2.plot(time,((xmin-rsun)*1.e-8)+hight,label="hight")
ax2.plot(time,B_energy,label="B")
ax2.plot(time,Convection_energy,label="convection")
ax2.set_yscale('log')
ax2.set_xlabel('x [Mm]')
ax2.set_ylabel('energy $[g/cm/s^2]$')
plt.legend()
"""
plt.tight_layout()

