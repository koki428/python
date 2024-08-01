#磁場の画像を作成（binary）
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
caseid = 0
caseid = input()
caseid = "d"+caseid.zfill(3)
    
datadir="../run/"+caseid+"/data/"
pngdir="../figs/"+caseid+"/png/"
os.makedirs(pngdir,exist_ok=True)

d = R2D2.R2D2_data(datadir)
for key in d.p:
    exec('%s = %s%s%s' % (key, 'd.p["',key,'"]'))
    
try:
    n0
except NameError:
    n0 = 0
if  n0 > d.p["nd"]:
    n0 = d.p["nd"]

print("Maximum time step= ",nd," time ="\
          ,dtout*float(nd)/3600./24.," [day]")

plt.close('all')

# read time
t0 = d.read_time(0,silent=True)

yran = ymax - ymin
xran = min(xmax-xmin,yran)

xsize = 12
ysize = 6
fig = plt.figure(num=1,figsize=(xsize,ysize))

for n in range(n0,nd):
#for n in range(0,1):
    print(n)
    ##############################
    # read time
    t = d.read_time(n,silent=True)
        

    ##############################
    # read value

    d.read_qq_2d(n,silent=True)
    ##############################

    #shading = "flat"
    shading = "gouraud"

    lfac = 1.e-8
    
    ax1 = fig.add_subplot(1,1,1,aspect='equal')
    #ax2 = fig.add_subplot(122,aspect='equal')

    sem = d.q2['se'].mean(axis=1)
    sem2, tmp = np.meshgrid(sem,y,indexing='ij')
    serms = np.sqrt(((d.q2['se'] - sem2)**2).mean(axis=1))
    serms2, tmp = np.meshgrid(serms,y,indexing='ij')
    
    
    #ax1.pcolormesh(y*1.e-8,(x-rsun)*1.e-8,(d.q2['se']-sem2)/serms2,vmin=-4,vmax=4)
    #ax1.contour(y*1.e-8,(x-rsun)*1.e-8,d.q2['tu'],levels=[1.])

    """
    ax=ax1.pcolormesh(y*1.e-8,(x-rsun)*1.e-8,d.q2['vx'],vmax=1.e+5,vmin=-1.e+5,cmap='coolwarm',shading=shading)
    ax1.contour(y*1.e-8,(x-rsun)*1.e-8,d.q2['tu'],levels=[1.])
    divider=make_axes_locatable(ax1)
    ax_cb=divider.new_horizontal(size="5%", pad=0.1)
    fig.add_axes(ax_cb)
    fig.colorbar(ax,cax=ax_cb)
    """

    bz = np.sqrt(d.q2['bz']**2)

    ax1.pcolormesh(y*1.e-8,(x-rsun)*1.e-8,bz,cmap='binary',shading=shading)
    ax1.contour(y*1.e-8,(x-rsun)*1.e-8,d.q2['tu'],levels=[1.])
    
#    ax1.annotate(s="t="+"{:.2f}".format((t-t0)/60/60)+" [hour]"\
#                     ,xy=[0.05,0.03],xycoords="figure fraction"\
#                     ,fontsize=18,color='black')#,bbox=bbox_props)

    if(n == n0):
        fig.tight_layout(pad=0.5)
    
    plt.pause(0.1)
    plt.savefig(pngdir+"py_bz"+'{0:08d}'.format(n)+".png")

    if(n != nd):
        clf()
