import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import R2D2
import sys
import os

#try:
#    caseid
#except NameError:
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
ysize=8
fig=plt.figure(num=1,figsize=(xsize,ysize))

print('start step number')
n0=input()
n0=int(n0)

print('final step number')
nd=input()
nd=int(nd)

dx=(xmax-xmin)/jx
dy=(ymax-ymin)/ix

for n in range(n0,nd+1):
    print(n)

    t=d.read_time(n,silent=True)

    d.read_qq_2d(n,silent=True)

    shading='gouraud'
   # shading='flat'

    lfac=1.e-8

    ax1=fig.add_subplot(111,aspect='equal')

#    bz=np.sqrt(d.q2['bz']**2)
    vydx=np.array(np.gradient(d.q2['vy'],dx,axis=1))
    vxdy=np.array(np.gradient(d.q2['vx'],dy,axis=1))
    wz=vydx-vxdy

   # ax=ax1.pcolormesh(y*1.e-8,(x-rsun)*1.e-8,log10(bz),vmax=4,vmin=0,cmap='gist_stern',shading=shading)
   # ax1.contour(y*1.e-8,(x-rsun)*1.e-8,d.q2['tu'],levels=[1.])
    ax=ax1.pcolormesh(y*lfac,(x-rsun)*lfac,wz,cmap='gray',vmax=0.00004,vmin=-0.00004,shading=shading)
    plt.ylabel('X [Mm]')
    plt.xlabel('Y [Mm]')
    divider=make_axes_locatable(ax1)
    ax_cb=divider.new_horizontal(size='5%',pad=0.1)
    fig.add_axes(ax_cb)
    fig.colorbar(ax,label='$\omega$ [/s]',cax=ax_cb)

    if(n == n0):
        fig.tight_layout(pad=0.1)

    plt.pause(0.1)
    plt.savefig(pngdir+"py"+"_k.k"+"_wz"+"{0:08d}".format(n)+".png")

    if(n != nd):
        clf()
