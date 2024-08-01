import numpy as np
import math
import sympy as sym
import matplotlib.pyplot as plt
import R2D2
import sys
import os

print("input caseid id")
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

print("Maximum time step=",nd,"time="\
        ,dtout*float(nd)/3600./24.,"[day]")

plt.close("all")

t0=d.read_time(0,silent=True)

dx=(xmax-xmin)/ix
dy=(ymax-ymin)/jx
fig=plt.figure(figsize=(8,10))
lfac=1.e-8

for n in range(n0,100):
    print(n)

    t=d.read_time(n,silent=True)
    d.read_qq_2d(n,silent=True)

    bx=d.q2['bx']
    by=d.q2['by']
   # bz=d.q2['bz']

    #psi cal
    psi_x=np.zeros((ix,jx))
    psi_y=np.zeros((ix,jx))
    psi=np.zeros((ix,jx))

    psi_x[0,0]=0.0
    j=0
    for i in range(1,ix):
        psi_x[i,j]=psi_x[i-1,j]+by[i,j]*dx
    for j in range(1,jx):
        i=0
        psi_x[i,j]=psi_x[i,j-1]-bx[i,j]*dy
        for i in range(1,ix):
            psi_x[i,j]=psi_x[i-1,j]+by[i,j]*dx

    psi_y[0,0]=0.0
    i=0
    for j in range(1,jx):
        psi_y[i,j]=psi_y[i,j-1]-bx[i,j]*dy
    for i in range(1,ix):
        j=0
        psi_y[i,j]=psi_y[i-1,j]+by[i,j]*dx
        for j in range(1,jx):
            psi_y[i,j]=psi_y[i,j-1]-bx[i,j]*dy

    for j in range(0,jx):
        for i in range(0,ix):
            psi[i,j]=(psi_x[i,j]+psi_y[i,j])*0.5

    ax4=fig.add_subplot(111,aspect='equal')
    ax4=ax4.pcolormesh(y*lfac,(x-rsun)*lfac,psi,vmin=0,cmap='gist_stern')
    #ax44=plt.contour(y*lfac,(x-rsun)*lfac,psi,colors='g',levels=[0.5e+12,1.e+12,1.5e+12,2.e+12,2.5e+12,3.e+12,3.5e+12,4.e+12])
    #plt.clabel(ax44,colors='white')
    fig.colorbar(ax4)

    if (n==n0):
        fig.tight_layout(pad=0.1)

    plt.pause(0.1)
    plt.savefig(pngdir+"psi_"+"{0:08d}".format(n)+".png")
    if(n != 99):
        clf()

