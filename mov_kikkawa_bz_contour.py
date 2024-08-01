import numpy as np
import matplotlib.pyplot as plt
import math
import R2D2
import sys
import os

try:
    caseid
except NameError:
    print("input caseid id (3 digit)")
    caseid=0
    caseid=input()
    caseid="d"+caseid.zfill(3)

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
if  n0>d.p["nd"]:
    n0=d.p["nd"]

print("Maximum time step=",nd,"time="\
        ,dtout*float(nd)/3600./24.,"[day]")

plt.close("all")

t0=d.read_time(0,silent=True)

yran=ymax-ymin
xran=min(xmax-xmin,yran)

#xsize=12
#ysize=xsize*yran/yran/2
xsize=8
ysize=10
fig=plt.figure(num=1,figsize=(xsize,ysize))

for n in range(n0,nd+1):
    print(n)

    t=d.read_time(n,silent=True)

    d.read_qq_2d(n,silent=True)

    shading='groroud'

    lfac=1.e-8

    ax1=fig.add_subplot(111,aspect='equal')

    bz=abs(d.q2['bz'])#/np.sqrt(2*4*math.pi*abs(d.q2['pr']))

   # ax=ax1.pcolormesh(y*1.e-8,(x-rsun)*1.e-8,log10(bz),vmax=4,vmin=0,cmap='gist_stern',shading=shading)
    ax=ax1.pcolormesh(y*1.e-8,(x-rsun)*1.e-8,bz,cmap='gist_stern',shading=shading)
   # ax1.contour(y*1.e-8,(x-rsun)*1.e-8,d.q2['tu'],levels=[1.])
    ax2=plt.contour(y*lfac,(x-rsun)*lfac,bz,colors='g')
    ax2.clabel(fmt='%1.1f',fontsize=8,colors='white')
    fig.colorbar(ax)

    if(n == n0):
        fig.tight_layout(pad=0.1)

    plt.pause(0.1)
    plt.savefig(pngdir+"py"+"_k.k"+"_bz_contour"+"{0:08d}".format(n)+".png")

    if(n != nd):
        clf()
