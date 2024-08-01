#メインチューブ内の運動方程式の各成分の２次元マップの画像

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
import matplotlib.gridspec as gridspec
from scipy.ndimage import gaussian_filter

print('roll number')
roll=input()
roll=int(roll)

print("input caseid id")
caseid=0
caseid=input()
caseid='d'+caseid.zfill(3)

datadir="../run/"+caseid+"/data/"
pngdir="../figs/"+caseid+"/png/"

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
dz=(zmax-zmin)/kx

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


fig=plt.figure(figsize=(12,9))
#fig=plt.rcParams['figure.figsize']=(12,12)
gs=gridspec.GridSpec(4,2,height_ratios=(1,1,1,1),width_ratios=(3,1))
n0=0
for n in range(n0,nd+1):
    
    #ro0=np.zeros(ix)
    #pr0=np.zeros(ix)
    #gx=np.zeros(ix)

    d=R2D2.R2D2_data(datadir)
    print(n)
    t=d.read_time(n,silent=True)
    d.read_qq_2d(n,silent=True)

    bx=d.q2['bx']
    by=d.q2['by']
    bz=d.q2['bz']
    vx=d.q2['vx']
    vy=d.q2['vy']
    vz=d.q2['vz']
    ro=d.q2['ro']
    se=d.q2['se']
    #pr=d.q2['pr']
    #print(np.shape(ro0))
    #ro0,tmp=np.meshgrid(ro0,y,indexing='ij')
    #pr0,tmp=np.meshgrid(pr0,y,indexing='ij')
    Gx,tmp=np.meshgrid(gx,y,indexing='ij')
    DprDro,tmp=np.meshgrid(dprdro,y,indexing='ij')
    DprDse,tmp=np.meshgrid(dprdse,y,indexing='ij')

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
    #pr=np.roll(pr,roll,axis=1)
    se=np.roll(se,roll,axis=1)
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
    
    #メインチューブの有効半径
    r_eff=0
    r_eff=np.sqrt(stats[1,4]*dx*dy/math.pi)
    #メインチューブの中心
    cy=centroids[1,0]
    cx=centroids[1,1]
    
    #print(gx)
    ##運動方程式の各要素の計算
    #浮力の計算
    mag_bouyancy=np.zeros((ix,jx))
    pr=DprDro*ro+DprDse*se

    for i in range(2,ix-2):
        for j in range(2,jx-2):
            mag_bouyancy[i,j]=-(-pr[i+2,j]+8*pr[i+1,j]-8*pr[i-1,j]+pr[i-2,j])/(12*dx)-ro[i,j]*Gx[i,j]
            #mag_bouyancy[i,j]=-ro[i,j]*gx[i]
            #mag_bouyancy[i,j]=-(pr[i+1,j]-pr[i-1,j])/(2*dx)-ro[i,j]*gx[i]
    
    #ローレンツ力の計算
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
    
    #慣性力の計算
    mag_inertia=np.zeros((ix,jx))
    for i in range(2,ix-2):
        for j in range(2,jx-2):
            mag_inertia[i,j]=-ro[i,j]*(vx[i,j]*(-vx[i+2,j]+8*vx[i+1,j]-8*vx[i+1,j]+vx[i-2,j])/(12*dx)\
                                      +vy[i,j]*(-vx[i,j+2]+8*vx[i,j+1]-8*vx[i,j-1]+vx[i,j-2])/(12*dy))
    
    th_bouyancy=np.zeros((ix,jx))
    th_lorentz=np.zeros((ix,jx))
    th_inertia=np.zeros((ix,jx))
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
    rr=np.zeros((ix,jx))
    for i in range(0,ix):
        for j in range(0,jx):
            rr[i,j]=np.sqrt((cx*dx-i*dx)**2+(cy*dy-j*dy)**2)
    
    mag_bouyancy=pd.DataFrame(mag_bouyancy)
    mag_lorentz=pd.DataFrame(mag_lorentz)
    mag_inertia=pd.DataFrame(mag_inertia)
    """
    th_bouyancy=pd.DataFrame(th_bouyancy)
    th_lorentz=pd.DataFrame(th_lorentz)
    th_inertia=pd.DataFrame(th_inertia)
    """
    # mag_bouyancy[psi < center]=0
    # mag_lorentz[psi < center]=0
    # mag_inertia[psi < center]=0
    """
    th_bouyancy[psi < center]=0
    th_lorentz[psi < center]=0
    th_inertia[psi < center]=0
    """
    

    ##プロット
    shading="gouraud"
    #ax1=fig.add_subplot(3,2,1)
    ax1=fig.add_subplot(gs[0,0])
    p=ax1.pcolormesh(y*1.e-8,(x-rsun)*1.e-8,mag_bouyancy-th_bouyancy,vmax=np.max(abs(mag_bouyancy)),vmin=-np.max(abs(mag_bouyancy)),cmap='seismic',shading=shading)
    ax1.set_title('bouyancy')
    ax1.set_xlabel('Y [Mm]')
    ax1.set_ylabel('X [Mm]')
    fig.colorbar(p,ax=ax1,label='$[g/cm^2/s^2]$')

    #ax2=fig.add_subplot(3,2,3)
    ax2=fig.add_subplot(gs[1,0])
    p=ax2.pcolormesh(y*1.e-8,(x-rsun)*1.e-8,mag_lorentz-th_lorentz,vmax=np.max(abs(mag_lorentz)),vmin=-np.max(abs(mag_lorentz)),cmap='seismic',shading=shading)
    ax2.set_title('lorentz')
    ax2.set_xlabel('Y [Mm]')
    ax2.set_ylabel('X [Mm]')
    fig.colorbar(p,ax=ax2,label='$[g/cm^2/s^2]$')

    ax3=fig.add_subplot(gs[2,0])
    p=ax3.pcolormesh(y*1.e-8,(x-rsun)*1.e-8,mag_inertia-th_inertia,vmax=np.max(abs(mag_inertia)),vmin=-np.max(abs(mag_inertia)),cmap='seismic',shading=shading)
    ax3.set_title('inertial')
    ax3.set_xlabel('Y [Mm]')
    ax3.set_ylabel('X [Mm]')
    fig.colorbar(p,ax=ax3,label='$[g/cm^2/s^2]$')

    #ax3=fig.add_subplot(3,2,5)
    ax4=fig.add_subplot(gs[3,0])
    p=ax4.pcolormesh(y*1.e-8,(x-rsun)*1.e-8,(mag_bouyancy-th_bouyancy)+(mag_lorentz-th_lorentz)+(mag_inertia-th_inertia),vmax=np.max(abs(mag_bouyancy+mag_lorentz+mag_inertia)),vmin=-np.max(abs(mag_bouyancy+mag_lorentz+mag_inertia)),cmap='seismic',shading=shading)
    ax4.set_title('total')
    ax4.set_xlabel('Y [Mm]')
    ax4.set_ylabel('X [Mm]')
    fig.colorbar(p,ax=ax4,label='$[g/cm^2/s^2]$')

    #ax4=fig.add_subplot(3,2,2,aspect='equal')
    ax5=fig.add_subplot(gs[0,1],aspect='equal')
    ax5.pcolormesh(y*1.e-8,(x-rsun)*1.e-8,mag_bouyancy-th_bouyancy,vmax=np.max(abs(mag_bouyancy)),vmin=-np.max(abs(mag_bouyancy)),cmap='seismic',shading=shading)
    ax5.set_title('bouyancy')
    ax5.set_xlabel('Y [Mm]')
    ax5.set_ylabel('X [Mm]')
    #fig.colorbar(p,ax=ax1,label='$[g/cm^2/s^2]$')
    ax5.set_xlim((cy*dy-3*r_eff)*1.e-8,(cy*dy+3*r_eff)*1.e-8)
    ax5.set_ylim((xmin+cx*dx-rsun-3*r_eff)*1.e-8,(xmin+cx*dx-rsun+3*r_eff)*1.e-8)

    ax6=fig.add_subplot(gs[1,1],aspect='equal')
    ax6.pcolormesh(y*1.e-8,(x-rsun)*1.e-8,mag_lorentz-th_lorentz,vmax=np.max(abs(mag_lorentz)),vmin=-np.max(abs(mag_lorentz)),cmap='seismic',shading=shading)
    ax6.set_title('lorentz')
    ax6.set_xlabel('Y [Mm]')
    ax6.set_ylabel('X [Mm]')
    ax6.set_xlim((cy*dy-3*r_eff)*1.e-8,(cy*dy+3*r_eff)*1.e-8)
    ax6.set_ylim((xmin+cx*dx-rsun-3*r_eff)*1.e-8,(xmin+cx*dx-rsun+3*r_eff)*1.e-8)

    ax7=fig.add_subplot(gs[2,1],aspect='equal')
    ax7.pcolormesh(y*1.e-8,(x-rsun)*1.e-8,mag_inertia-th_inertia,vmax=np.max(abs(mag_inertia)),vmin=-np.max(abs(mag_inertia)),cmap='seismic',shading=shading)
    ax7.set_title('inertial')
    ax7.set_xlabel('Y [Mm]')
    ax7.set_ylabel('X [Mm]')
    ax7.set_xlim((cy*dy-3*r_eff)*1.e-8,(cy*dy+3*r_eff)*1.e-8)
    ax7.set_ylim((xmin+cx*dx-rsun-3*r_eff)*1.e-8,(xmin+cx*dx-rsun+3*r_eff)*1.e-8)

    ax8=fig.add_subplot(gs[3,1],aspect='equal')
    ax8.pcolormesh(y*1.e-8,(x-rsun)*1.e-8,(mag_bouyancy-th_bouyancy)+(mag_lorentz-th_lorentz)+(mag_inertia-th_inertia),vmax=np.max(abs(mag_bouyancy+mag_lorentz+mag_inertia)),vmin=-np.max(abs(mag_bouyancy+mag_lorentz+mag_inertia)),cmap='seismic',shading=shading)
    ax8.set_title('total')
    ax8.set_xlabel('Y [Mm]')
    ax8.set_ylabel('X [Mm]')
    ax8.set_xlim((cy*dy-3*r_eff)*1.e-8,(cy*dy+3*r_eff)*1.e-8)
    ax8.set_ylim((xmin+cx*dx-rsun-3*r_eff)*1.e-8,(xmin+cx*dx-rsun+3*r_eff)*1.e-8)

    if (n == n0):
        fig.tight_layout(pad=1)
    
    plt.pause(0.1)
    plt.savefig(pngdir+'force_2D_main_v4'+'{0:08d}'.format(n)+'.png')

    if (n != nd):
        clf()

    

