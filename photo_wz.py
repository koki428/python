import numpy as np
import matplotlib.pyplot as plt
import math
import R2D2
import sys
import os

print("input caseid")
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
if nd>d.p["nd"]:
    n0=d.p["nd"]

print("Maximum time step=",nd,"time="\
        ,dtout*float(nd)/3600./24.,"[day]")

plt.close("all")

t0=d.read_time(0,silent=True)

xsize=8
ysize=10
fig=plt.figure(num=1,figsize=(xsize,ysize))

dx=(xmax-xmin)/ix
dy=(ymax-ymin)/jx

#print("input shot number")
n=input('input shot number : ')
#print(n)
n=int(n)

t=d.read_time(n,silent=True)
d.read_qq_2d(n,silent=True)

shading='gouraud'
lfac=1.e-8

ax=fig.add_subplot(111,aspect='equal')
#bz=abs(d.q2['bz'])
vydx=np.array(np.gradient(d.q2['vy'],dx,axis=1))
vxdy=np.array(np.gradient(d.q2['vx'],dy,axis=1))
wz=vydx-vxdy

ax=ax.pcolormesh(y*lfac,(x-rsun)*lfac,wz,cmap='gray',vmax=0.00005,vmin=-0.00005,shading=shading)
fig.colorbar(ax)

plt.savefig(pngdir+'py'+'_k.k'+'_wz'+'{0:08d}'.format(n)+'.png')
