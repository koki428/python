#上部に達したステップ数を出力

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
from decimal import Decimal, ROUND_HALF_UP, ROUND_HALF_EVEN

print('roll number')
roll=input()
roll=int(roll)

print('input caseid id')
caseid=0
caseid=input()
caseid="d"+caseid.zfill(3)

datadir="../run/"+caseid+"/data/"
pngdir="../figs/"+caseid+"/png/"

d=R2D2.R2D2_data(datadir)
for key in d.p:
    exec('%s=%s%s%s' % (key, 'd.p["',key,'"]'))

print("Maximum time step= ",nd," time ="\
            ,dtout*float(nd)/3600./24.," [day]")

try:
    n0
except NameError:
    n0=0
if n0>d.p["nd"]:
    n0=d.p["nd"]

plt.close('all')

t0=d.read_time(0,silent=True)

dx=(xmax-xmin)/ix
dy=(ymax-ymin)/jx

##t=nの計算
found=False
for n in range(n0,nd+1):
    #print(n)
    t=d.read_time(n,silent=True)
    d.read_qq_2d(n,silent=True)
    bx=d.q2['bx']
    by=d.q2['by']
    # bz=d.q2['bz']

    #平行移動
    bx=np.roll(bx,roll,axis=1)
    by=np.roll(by,roll,axis=1)
    # bz=np.roll(bz,roll,axis=1)

    #t=nでのpsiの計算
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
    
    #中心の高さの計算
    cx=centroids[1,1]
    hight=(xmin-rsun+cx*dx)*1.e-8

    rr=np.sqrt(stats[1,4]*dx*dy/math.pi)*1.e-8

    #上部の到達したか判定
    if np.any(psi_thr[475,:] > 0):
        print('step',n)
        found=True
    
    if found:
        break

