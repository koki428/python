#磁束の保持率の比較
import numpy as np
import math
import sympy as sym
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import R2D2
import cv2
import sys
import os
from decimal import Decimal, ROUND_HALF_UP, ROUND_HALF_EVEN

"""
print('compare case count')
count=0
count=input()
count=int(count)
"""
flux=np.zeros(8)
ratio=np.zeros(8)

print('roll number')
roll=input()
roll=int(roll)

for n in range(0,8):
    print("input caseid id",n+1)
    caseid=0
    caseid=input()
    caseid='d'+caseid.zfill(3)

    datadir="../run/"+caseid+"/data/"
    pngdir="../figs/"+caseid+"/png/"

    d=R2D2.R2D2_data(datadir)
    for key in d.p:
        exec('%s=%s%s%s' % (key, 'd.p["',key,'"]'))

    try:
        n0
    except NameError:
        n0=0
    if n0>d.p["nd"]:
        n0=d.p["nd"]

    plt.close("all")

    t0=d.read_time(0,silent=True)

    dx=(xmax-xmin)/ix
    dy=(ymax-ymin)/jx

        
    print('input max step number')
    nd=input()
    nd=int(nd)

    if (nd != 0):
        ##t＝0の計算
        t=d.read_time(0,silent=True)
        d.read_qq_2d(0,silent=True)
        bz=d.q2['bz']

        #平行移動
        bz=np.roll(bz,roll,axis=1)

        #t=0での磁束の計算
        initial_flux=0.0
        for j in range(0,jx):
            for i in range(0,ix):
                initial_flux=initial_flux+bz[i,j]*dx*dy

        ##t=ndの計算
        t=d.read_time(nd-1,silent=True)
        d.read_qq_2d(nd-1,silent=True)
        bx=d.q2['bx']
        by=d.q2['by']
        bz=d.q2['bz']

        #平行移動
        bx=np.roll(bx,roll,axis=1)
        by=np.roll(by,roll,axis=1)
        bz=np.roll(bz,roll,axis=1)

        #t=ndでのpsiの計算
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

        #2値化
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
            if (abs(psi_max-psi_min) < 1.e+1):
                break
        
        #t=nでの磁束の計算
        for j in range(0,jx):
            for i in range(0,ix):
                if (psi[i,j] > center):
                    flux[n]=flux[n]+bz[i,j]*dx*dy
        
        ##保持率の計算
        ratio[n]=flux[n]/initial_flux*100
        print('caseid',caseid,'ratio =',ratio[n])
    if (nd == 0):
        print('caseid',caseid,'ratio =',ratio[n])

##plot
x=[0.00,0.10,0.15,0.20,0.25,0.30,0.35,0.40]
#x=[1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7]
y=np.zeros(8)
yy=np.zeros(len(y))
for n in range(0,8):
    y[n]=ratio[n]
    yy[n]=Decimal(str(y[n])).quantize(Decimal('0.1'),ROUND_HALF_UP)
fig=plt.figure(figsize=(10,6))
ax=fig.add_subplot(1,1,1)
ax=plot(x,y,marker='D',linestyle='None')
plt.xlim(0,0.45)
plt.ylim(0,100)
plt.ylabel('flux retained (%)',fontsize=22)
plt.xlabel('$\lambda$',fontsize=22)
#plt.xlabel('x $10^5$ [G]',fontsize=22)
for n in range(0,len(y)):
    plt.text(x[n],y[n],yy[n],fontsize=16)