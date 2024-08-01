#熱対流有と無の上昇速度の比較

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

##熱対流有の計算
print("input convection caseid id")
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
nd1=input()
nd1=int(nd1)

t=d.read_time(nd1-1,silent=True)
d.read_qq_2d(nd1-1,silent=True)

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
hight1=np.zeros(nd1)
velocity1=np.zeros(nd1)

for n in range(n0,nd1):
    print(n)
    t=d.read_time(n,silent=True)
    d.read_qq_2d(n,silent=True)
    
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
    hight1[n]=cx*dx*1.e-8
    print('hight1[',n,'] =',hight1[n])
    #速度の計算
    if (n != 0):
        velocity1[n]=(hight1[n]-hight1[n-1])/(8*60*60)*1.e+3
        #velocity[n]=vx[int(cx),int(cy)]*1.e-5
        print('velocity1[',n,'] =',velocity1[n])


##熱対流なしの計算    
print("input no convection caseid id")
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
nd2=input()
nd2=int(nd2)

t=d.read_time(nd2-1,silent=True)
d.read_qq_2d(nd2-1,silent=True)

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
hight2=np.zeros(nd2)
velocity2=np.zeros(nd2)

for n in range(n0,nd2):
    print(n)
    t=d.read_time(n,silent=True)
    d.read_qq_2d(n,silent=True)
    
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
    hight2[n]=cx*dx*1.e-8
    print('hight2[',n,'] =',hight2[n])
    #速度の計算
    if (n != 0):
        velocity2[n]=(hight2[n]-hight2[n-1])/(8*60*60)*1.e+3
        #velocity[n]=vx[int(cx),int(cy)]*1.e-5
        print('velocity2[',n,'] =',velocity2[n])

##熱対流の計算
convection_rms=np.zeros((nd1,ix))
convection_mean=np.zeros(ix)
d=R2D2.R2D2_data("../run/d004/data/")
t=d.read_time(300+n,silent=True)
d.read_qq_2d(300+n,silent=True)
for n in range(0,nd1):
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

##plot
fig=plt.figure(figsize=(8,8))
ax1=fig.add_subplot(1,1,1)
ax1.plot(((xmin-rsun)*1.e-8)+hight1,velocity1,label="velocity(convection)")
ax1.plot(((xmin-rsun)*1.e-8)+hight2,velocity2,label="velocity(no convection)")
ax1.plot((x-rsun)*1.e-8,convection_mean,label="convection velocity")
plt.legend()
ax1.set_xlabel('x [Mm]')
ax1.set_ylabel('vx [km/s]')