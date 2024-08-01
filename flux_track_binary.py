import numpy as np
import math
import sympy as sym
import matplotlib.pyplot as plt
import R2D2
import cv2
import sys
import os

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

print("Maximum time step=",nd,"time="\
        ,dtout*float(nd)/3600./24.,"[day]")

plt.close("all")

t0=d.read_time(0,silent=True)

dx=(xmax-xmin)/jx
dy=(ymax-ymin)/ix

print('input step number')
n=input()
n=int(n)

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

#for j in range(1,jx):
#    psi[i,j]=psi[i,j-1]+by[i,j]*dx
#for i in range(1,ix):
#    j=0
#    psi[i,j]=psi[i-1,j]-bx[i,j]*dy
#    for j in range(1,jx):
#        psi[i,j]=psi[i,j-1]+by[i,j]*dx

#2値化
psi_thr=np.zeros((ix,jx))
psi_max=np.max(psi)
psi_min=100
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
    if (abs(psi_max-psi_min) < 1.e+8):
        break

##plot
fig=plt.figure(figsize=(8,10))
lfac=1.e-8
ax=fig.add_subplot(1,1,1,aspect='equal')
ax=pcolormesh(y*lfac,(x-rsun)*lfac,binary_psi,cmap='gray')
plt.ylabel('X [Mm]')
plt.xlabel('Y [Mm]')

plt.savefig(pngdir+'py'+'_psi_binary_track'+'{0:08d}'.format(n)+'.png')
