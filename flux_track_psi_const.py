#最後の時間でのpsiで閾値を決定し、それで追跡
import numpy as np
import math
import sympy as sym
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import R2D2
import cv2
import sys
import os
#import matplotlib
#matplotlib.use("Agg")

print("input caseid id")
caseid=0
caseid=input()
caseid='d'+caseid.zfill(3)

datadir="../run/"+caseid+"/data/"
pngdir="../figs/"+caseid+"/png/"
os.makedirs(pngdir,exist_ok=True)

d=R2D2.R2D2_data(datadir)
for key in d.p:
    exec('%s=%s%s%s' % (key, 'd.p["',key,'"]'))

try:
    n0
except NameError:
    n0=0
if n0>d.p["nd"]:
    n0=d.p["nd"]

print("Maximum time step =",nd,"time="\
        ,dtout*float(nd)/3600./24.,"[day]")

plt.close("all")

print('input start step number')
n0=input()
n0=int(n0)

t0=d.read_time(n0,silent=True)

dx=(xmax-xmin)/jx
dy=(ymax-ymin)/ix

print('input max step number')
nd=input()
nd=int(nd)

##t＝ndの計算
t=d.read_time(nd,silent=True)
d.read_qq_2d(nd,silent=True)
bx=d.q2['bx']
by=d.q2['by']
bz=d.q2['bz']

##t=ndでのpsiの計算
psi=np.zeros((ix,jx))
psi[0,0]=0.0
i=0
for j in range(1,jx):
    psi[i,j]=psi[i,j-1]-bx[i,j]*dy
for i in range(1,ix):
    j=0
    psi[i,j]=psi[i-1,j]+by[i,j]*dx
    for j in range(1,jx):
        psi[i,j]=psi[i,j-1]-bx[i,j]*dy

##2値化
psi_thr=np.zeros((ix,jx))
psi_max=np.max(psi)
psi_min=100.0
#print('psi_max =',np.max(psi))
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
pp=1.2
center=center*pp

##t=ndでの磁束の計算
"""
initial_flux=0.0

for j in range(0,jx):
    for i in range(0,ix):
        initial_flux=initial_flux+bz[i,j]*dx*dy
#print('step 0 flux = ',initial_flux)
"""
##t=nの計算
#print('input step number')
#n=input()
#n=int(n)

lfac=1.e-8
fig=plt.figure(figsize=(12,6))

flux=np.zeros(nd+1)
ratio=np.zeros(nd+1)
area=np.zeros(nd+1)
b_total=np.zeros(nd+1)
b_min=np.zeros(nd+1)

for n in range(n0,nd+1):
    print('step',n)
    t=d.read_time(n,silent=True)
    d.read_qq_2d(n,silent=True)
    bx=d.q2['bx']
    by=d.q2['by']
    bz=d.q2['bz']

    #t=nでのpsiの計算
    psi=np.zeros((ix,jx))
    psi[0,0]=0.0
    i=0
    for j in range(1,jx):
        psi[i,j]=psi[i,j-1]-bx[i,j]*dy
    for i in range(1,ix):
        j=0
        psi[i,j]=psi[i-1,j]+by[i,j]*dx
        for j in range(1,jx):
            psi[i,j]=psi[i,j-1]-bx[i,j]*dy


    #2値化
    psi_thr=np.zeros((ix,jx))
    #psi_max=np.max(psi)
    #psi_min=100.0
    #print('psi_max =',np.max(psi))
    #while True:
        #center=(psi_max+psi_min)*0.5
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
        for m in range(1,nlabels-1):
            if (stats[1,4] < stats[m,4]):
                stats[1,4]=stats[m,4]

    area[n]=stats[1,4]
    #print(area[n])

    #t=nでの磁束の計算
    #flux=0.0
    #count=0
    b_min[n]=1.e+5
    for j in range(0,jx):
        for i in range(0,ix):
            if (psi[i,j] > center):
                flux[n]=flux[n]+bz[i,j]*dx*dy
                if (b_min[n] > bz[i,j]):
                    b_min[n]=bz[i,j]
    #print(count)

    #print('step',n,'flux =',flux[n])
    
    """
    ##初期の磁束の計算
    for j in range(0,jx):
        for i in range(0,ix):
            flux[n0]=flux[n0]+bz[i,j]*dx*dy

    ##保持率の計算
    ratio[n]=flux[n]/flux[n0]*100
    print('ratio =',ratio[n])
    """
    """
    ##Bの総和の計算
    for j in range(0,jx):
        for i in range(0,ix):
            if (psi[i,j] > center):
                b_total[n]=b_total[n]+bz[i,j]
    
    flux[n]=area[n]*b_total[n]

    ratio[n]=flux[n]/flux[n0]*100
    """

    #メインチューブの追跡画像
    ax1=fig.add_subplot(1,1,1,aspect='equal')
    ax2=ax1.pcolormesh(y*lfac,(x-rsun)*lfac,psi,vmin=0,cmap='gist_stern')
    ax3=plt.contour(y*lfac,(x-rsun)*lfac,psi,colors='g',levels=[center])
    plt.ylabel('X [Mm]')
    plt.xlabel('Y [Mm]')
    divider=make_axes_locatable(ax1)
    ax_cb=divider.new_horizontal(size='5%',pad=0.1)
    fig.add_axes(ax_cb)
    fig.colorbar(ax2,label='$\psi$ [G$\cdot$cm]',cax=ax_cb)

    pp=str(pp)
    if(n == n0):
        fig.tight_layout(pad=0.8)

    plt.pause(0.1)
    plt.savefig(pngdir+"py_main_track_psi_const_"+pp+_+"{0:08d}".format(n)+".png")

    if (n != nd):
        clf()

#print(np.max(ratio[nd]))
#print(np.min(ratio))

##初期の磁束の計算
flux_initial=0
for j in range(0,jx):
    for i in range(0,ix):
        flux_initial=flux_initial+bz[i,j]*dx*dy

##保持率の計算
for n in range(0,nd+1):
    ratio[n]=flux[n]/flux_initial*100
    print('step',n,'ratio =',ratio[n])

##plot
x=np.linspace(0,nd*8,nd+1)
fig=plt.figure(figsize=(12,6))
ax1=fig.add_subplot(1,1,1)
ax1=plot(x,ratio,marker='D',linestyle='None')
#plt.ylim(90,102)
plt.xlabel('time (h)',fontsize=22)
plt.ylabel('flux retained (%)',fontsize=22)
"""
ax2=fig.add_subplot(3,1,2)
ax2=plot(x,b_min,marker='D',linestyle='None')
#plt.ylim(90,102)
#plt.ylabel('flux retained (%)',fontsize=22)
plt.xlabel('time (h)',fontsize=22)

ax3=fig.add_subplot(3,1,3)
ax3=plot(x,area,marker='D',linestyle='None')
plt.xlabel('time (h)',fontsize=22)
"""
