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

#print("input shot number")
n=input('input step number : ')
#print(n)
n=int(n)

t=d.read_time(n,silent=True)
d.read_qq_2d(n,silent=True)

shading='gouraud'
lfac=1.e-8

ax=fig.add_subplot(111,aspect='equal')
bz=np.sqrt((d.q2['bz'])**2)
ax=ax.pcolormesh(y*lfac,(x-rsun)*lfac,bz,cmap='gist_stern',shading=shading)
fig.colorbar(ax)

plt.savefig(pngdir+'py'+'_k.k'+'_bz'+'{0:08d}'.format(n)+'.png')
