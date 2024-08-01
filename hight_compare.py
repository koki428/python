#cheung2006 高さと物理量の関係

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

print('input start step number')
n0=input()
n0=int(n0)

#メインチューブの定義
print('input max step number')
nd=input()
nd=int(nd)

t=d.read_time(nd,silent=True)
d.read_qq_2d(nd,silent=True)

bx=d.q2['bx']
by=d.q2['by']
bz=d.q2['bz']

#psi cal
psi_x=np.zeros((ix,jx))
psi_y=np.zeros((ix,jx))
psi=np.zeros((ix,jx))

#psi_x[0,0]=0.0
#j=0
#for i in range(1,ix):
#    psi_x[i,j]=psi_x[i-1,j]+by[i,j]*dx
#for j in range(1,jx):
#    i=0
#    psi_x[i,j]=psi_x[i,j-1]-bx[i,j]*dy
#    for i in range(1,ix):
#        psi_x[i,j]=psi_x[i-1,j]+by[i,j]*dx

psi[0,0]=0.0
i=0
for j in range(1,jx):
    psi[i,j]=psi[i,j-1]-bx[i,j]*dy
for i in range(1,ix):
    j=0
    psi[i,j]=psi[i-1,j]+by[i,j]*dx
    for j in range(1,jx):
        psi[i,j]=psi[i,j-1]-bx[i,j]*dy

#for j in range(0,jx):
#    for i in range(0,ix):
#        psi[i,j]=(psi_x[i,j]+psi_y[i,j])*0.5

# 2値化を行う
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
    if (abs(psi_max-psi_min) < 0.5) and (nlabels-1==1):
        break

center=center*1.25

#物理量の計算
rad=np.zeros(nd+1-n0)
pre=np.zeros(nd+1-n0)
tem=np.zeros(nd+1-n0)
bet=np.zeros(nd+1-n0)
bzz=np.zeros(nd+1-n0)
bb=np.zeros(nd+1-n0)
xp=np.zeros(nd+1-n0)

for n in range(0,nd+1-n0):
    print(n+n0)
    t=d.read_time(n+n0,silent=True)
    d.read_qq_2d(n+n0,silent=True)

    bx=d.q2['bx']
    by=d.q2['by']
    bz=d.q2['bz']
    pr=d.q2['pr']
    te=d.q2['te']

    #psi cal
    psi_x=np.zeros((ix,jx))
    psi_y=np.zeros((ix,jx))
    psi=np.zeros((ix,jx))

    #psi_x[0,0]=0.0
    #j=0
    #for i in range(1,ix):
    #    psi_x[i,j]=psi_x[i-1,j]+by[i,j]*dx
    #for j in range(1,jx):
    #    i=0
    #    psi_x[i,j]=psi_x[i,j-1]-bx[i,j]*dy
    #    for i in range(1,ix):
    #        psi_x[i,j]=psi_x[i-1,j]+by[i,j]*dx

    psi[0,0]=0.0
    i=0
    for j in range(1,jx):
        psi[i,j]=psi[i,j-1]-bx[i,j]*dy
    for i in range(1,ix):
        j=0
        psi[i,j]=psi[i-1,j]+by[i,j]*dx
        for j in range(1,jx):
            psi[i,j]=psi[i,j-1]-bx[i,j]*dy

    #for j in range(0,jx):
    #    for i in range(0,ix):
    #        psi[i,j]=(psi_x[i,j]+psi_y[i,j])*0.5
    """
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
        if (abs(psi_max-psi_min) < 1.e+8) and (nlabels-1==1):
            break
    """
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

    rad[n]=np.sqrt(stats[1,4]*dx*dy/math.pi)

    cx=int(centroids[1,0])
    cy=int(centroids[1,1])

    print('(cx,cy) = (',cy,',',cx,')')
    #print('bz_max =',np.max(bz))
    print('bz_cen =',bz[cy,cx]) 
    
    idx=np.unravel_index(np.argmax(bz),bz.shape)
    cx=int(idx[0])
    cy=int(idx[1])

    print('(cxx,cyy) = (',cx,',',cy,')')
    print('bz_max =',np.max(bz))

    Pr0,tmp=np.meshgrid(pr0,y,indexing='ij')
    pre[n]=Pr0[cx,cy]+pr[cx,cy]

    Te0,tmp=np.meshgrid(te0,y,indexing='ij')
    tem[n]=Te0[cx,cy]+te[cx,cy]

    #idx=np.unravel_index(np.argmax(bz),bz.shape)
    #cxx=int(idx[0])
    #cyy=int(idx[1])
    #print(idx[1])
    bb[n]=np.sqrt(bx[cx,cy]**2+by[cx,cy]**2+bz[cx,cy]**2)
    #bb[n]=np.sqrt(bx[cy,cx]**2+by[cy,cx]**2+bz[cy,cx]**2)
    #bb[n]=np.max(bz)
    #idx=np.nuravel_index(np.argmax(bz),bz.shape)

    bet[n]=8*math.pi*pre[n]/(bb[n]**2)

    #bzz[n]=abs(bz[cy,cx])

    xp[n]=pre[n]/pre[0]

    #print(1/bb[n]**2-1/(bb[0]*xp[n]**(1/gam))**2)

#for n in range(1,13):
#    rad[n]=rad[13]
#    tem[n]=tem[13]
#    bet[n]=bet[13]
#    bzz[n]=bzz[13]
#    bb[n]=bb[13]
#    xp[n]=xp[13]

#プロット
gam=5/3
fig=plt.figure(figsize=(12,10))

ax1=fig.add_subplot(2,2,1)
plt.yscale("log")
#plt.xscale("log")
ax1.plot(xp,bet,marker='D',linestyle='None',label="simulation")
ax1.plot(xp,bet[0]*xp**((gam-2)/gam),label=r"$\beta_{0}\chi^{\gamma-2/\gamma}_{p}$")
#ax1.set_xlabel('Pressure Contract')
ax1.set_ylabel('β',fontsize=20)
ax1.legend(fontsize=16)

ax2=fig.add_subplot(2,2,2)
plt.yscale("log")
#plt.xscale("log")
ax2.plot(xp,bb,marker='D',linestyle='None',label="simulation")
ax2.plot(xp,bb[0]*xp**(1/gam),label=r"$B_{0}\chi^{1/\gamma}_{p}$")
#ax2.set_xlabel('Pressure Contract')
ax2.set_ylabel('|B|',fontsize=20)
ax2.legend(fontsize=16)

ax3=fig.add_subplot(2,2,3)
plt.yscale("log")
#plt.xscale("log")
ax3.plot(xp,tem,marker='D',linestyle='None',label="simulation")
ax3.plot(xp,tem[0]*xp**((gam-1)/gam),label=r"$T_{0}\chi^{\gamma-1/\gamma}_{p}$")
ax3.set_xlabel('Pressure Contract',fontsize=20)
ax3.set_ylabel('Temperature',fontsize=20)
ax3.legend(fontsize=16)

ax4=fig.add_subplot(2,2,4)
plt.yscale("log")
#plt.xscale("log")
ax4.plot(xp,rad,marker='D',linestyle='None',label="simulation")
ax4.plot(xp,rad[0]*xp**(-1/(2*gam)),label=r"$R_{0}\chi^{-1/2\gamma}_{p}$")
ax4.set_xlabel('Pressure Contract',fontsize=20)
ax4.set_ylabel('Radius',fontsize=20)
ax4.legend(fontsize=16)

"""
ax5=fig.add_subplot(2,2,3)
plt.yscale("log")
ax5.plot(xp,pre*8*math.pi,marker='D',linestyle='None')
ax5.set_xlabel('Pressure Contract')
ax5.set_ylabel('Pressure x 8π')

ax6=fig.add_subplot(2,2,4)
plt.yscale("log")
ax6.plot(xp,abs(1/bb**2-1/(bb[0]*xp**(1/gam))**2),marker='D',linestyle='None')
ax6.set_xlabel('Pressure Contract')
#ax6.set_ylabel('Pressure x 8π')
"""
