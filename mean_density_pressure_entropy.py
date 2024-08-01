##メインチューブ内の密度、圧力、エントロピーを規格化

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

print('roll nmber')
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
ro_norm_mean=np.zeros(nd+1)
pr_norm_mean=np.zeros(nd+1)
se_norm_mean=np.zeros(nd+1)
hight=np.zeros(nd+1)

for n in range(n0,nd+1):
    print(n)
    #磁場なしの計算
    d=R2D2.R2D2_data("../run/d004/data/")
    t=d.read_time(300+n,silent=True)
    d.read_qq_2d(300+n,silent=True)

    ro=d.q2['ro']
    se=d.q2['se']
    DprDro,tmp=np.meshgrid(dprdro,y,indexing='ij')
    DprDse,tmp=np.meshgrid(dprdse,y,indexing='ij')

    ro=np.roll(ro,roll,axis=1)
    se=np.roll(se,roll,axis=1)
    DprDro=np.roll(DprDro,roll,axis=1)
    DprDse=np.roll(DprDse,roll,axis=1)

    #圧力の計算
    pr=DprDro*ro+DprDse*se

    #密度、圧力、エントロピーの2乗
    ro_sq=ro**2
    se_sq=se**2
    pr_sq=pr**2

    #密度、圧力、エントロピーの水平方向平均の計算
    ro=pd.DataFrame(ro)
    se=pd.DataFrame(se)
    pr=pd.DataFrame(pr)

    ro_hd=ro.mean(axis=1)
    se_hd=se.mean(axis=1)
    pr_hd=pr.mean(axis=1)

    #密度、圧力、エントロピーのrmsの計算
    ro_sq=pd.DataFrame(ro_sq)
    # se_sq=pd.DataFrame(se_sq)
    # pr_sq=pd.DataFrame(pr_sq)

    ro_rms=np.sqrt(ro_sq.mean(axis=1))
    # se_rms=np.sqrt(se_sq.mean(axis=0))
    # pr_rms=np.sqrt(pr_sq.mean(axis=0))

    #磁場ありの計算
    d=R2D2.R2D2_data(datadir)
    t=d.read_time(n,silent=True)
    d.read_qq_2d(n,silent=True)

    bx=d.q2['bx']
    by=d.q2['by']
    ro=d.q2['ro']
    se=d.q2['se']
    DprDro,tmp=np.meshgrid(dprdro,y,indexing='ij')
    DprDse,tmp=np.meshgrid(dprdse,y,indexing='ij')
    
    bx=np.roll(bx,roll,axis=1)
    by=np.roll(by,roll,axis=1)
    ro=np.roll(ro,roll,axis=1)
    se=np.roll(se,roll,axis=1)
    DprDro=np.roll(DprDro,roll,axis=1)
    DprDse=np.roll(DprDse,roll,axis=1)

    #圧力の計算
    pr=DprDro*ro+DprDse*se

    #規格化した密度、圧力、エントロピーの計算
    ro_norm=np.zeros((ix,jx))
    pr_norm=np.zeros((ix,jx))
    se_norm=np.zeros((ix,jx))
    for i in range(0,ix):
        for j in range(0,jx):
            ro_norm[i,j]=(ro[i,j]-ro_hd[i])/ro_rms[i]
            pr_norm[i,j]=(pr[i,j]-pr_hd[i])/(ro_rms[i]*DprDro[i,j])
            se_norm[i,j]=(se[i,j]-se_hd[i])/ro_rms[i]*(-DprDse[i,j]/DprDro[i,j])

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

    hight[n]=cx*dx*1.e-8

    #規格化した密度、圧力、エントロピーの平均
    ro_norm=pd.DataFrame(ro_norm)
    pr_norm=pd.DataFrame(pr_norm)
    se_norm=pd.DataFrame(se_norm)

    ro_norm_mean[n]=ro_norm[psi > center].mean().mean()
    pr_norm_mean[n]=pr_norm[psi > center].mean().mean()
    se_norm_mean[n]=se_norm[psi > center].mean().mean()

    print('ro normalization =',ro_norm_mean[n])
    print('pr normalization =',pr_norm_mean[n])
    print('se normalization =',se_norm_mean[n])

#plot
time=np.linspace(0,nd,nd+1)
fig=plt.figure(figsize=(12,6))
plt.rcParams['mathtext.fontset']='cm'

ax1=fig.add_subplot(1,2,1)
ax1.plot(time,ro_norm_mean,label=r"$\tilde{\rho}$")
ax1.plot(time,pr_norm_mean,label=r"$\tilde{p}$")
ax1.plot(time,se_norm_mean,label=r'$\tilde{s}$')
# ax1.plot(time,total_mean,label="total")
ax1.set_xlabel('Time [h]',fontsize=15)
ax1.set_ylabel(r'$\tilde{\rho},\tilde{p},\tilde{s}$',fontsize=22)
# ax1.set_ylim(-0.1,0.1)
plt.legend()

ax2=fig.add_subplot(1,2,2)
ax2.plot(time,(xmin-rsun)*1.e-8+hight,label='hight')
ax2.set_xlabel('Time [h]',fontsize=15)
ax2.set_ylabel('Hight [Mm]',fontsize=15)
ax2.set_ylim((xmin-rsun)*1.e-8,(xmax-rsun)*1.e-8)
    
    
            

    



