#磁束管の上昇速度をプロット
import numpy as np
import matplotlib.pyplot as plt
import math
import R2D2
import sys
import os
import cv2

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
if n0 > d.p["nd"]:
    n0=d.p["nd"]

print('Maximum time step =',nd,'time ='\
        ,dtout*float(nd)/3600./24.,'[day]')

plt.close("all")

t0=d.read_time(0,silent=True)
d.read_qq_2d(0,silent=True)

dx=(xmax-xmin)/ix
dy=(ymax-ymin)/jx

#メインチューブの定義
print('input start step number')
n0=input()
n0=int(n0)

print('input max step number')
nd=input()
nd=int(nd)

t=d.read_time(nd,silent=True)
d.read_qq_2d(nd,silent=True)

bx=d.q2['bx']
by=d.q2['by']
bz=d.q2['bz']

#psiの計算
psi_x=np.zeros((ix,jx))
psi_y=np.zeros((ix,jx))
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

#各ステップの計算
hight=np.zeros(nd+1)
velocity=np.zeros(nd+1)

#n0ステップの計算
print(n0)
t=d.read_time(n0,silent=True)
d.read_qq_2d(n0,silent=True)

bx=d.q2['bx']
by=d.q2['by']
bz=d.q2['bz']

#psiの計算
psi_x=np.zeros((ix,jx))
psi_y=np.zeros((ix,jx))
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
hight[n0]=cx*dx*1.e-8
print('hight[',n0,'] =',hight[n0])

#n0+1~ndステップの計算
for n in range(n0+1,nd+1):
    print(n)
    t=d.read_time(n,silent=True)
    d.read_qq_2d(n,silent=True)
    
    bx=d.q2['bx']
    by=d.q2['by']
    bz=d.q2['bz']

    #psiの計算
    psi_x=np.zeros((ix,jx))
    psi_y=np.zeros((ix,jx))
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
    velocity[n]=(hight[n]-hight[n-1])/8*60*60
    #velocity[n]=vx[int(cx),int(cy)]*1.e-8
    print('velocity[',n,'] =',velocity[n])

#plot
time=np.linspace(0,480,nd+1)
fig=plt.figure(figsize=(8,8))

ax1=fig.add_subplot(2,1,1)
ax1.plot(time,velocity,marker='D',linestyle='None')
ax1.set_xlim(0,480)
#ax1.set_xlabel('Time (h)')
ax1.set_ylabel('Rise velocity (Mm/s)',fontsize=22)

ax2=fig.add_subplot(2,1,2)
ax2.plot(time,hight-hight[n0],marker='D',linestyle='None')
ax2.set_xlim(0,480)
ax2.set_xlabel('Time (h)',fontsize=20)
ax2.set_ylabel('Hight (Mm)',fontsize=20)
