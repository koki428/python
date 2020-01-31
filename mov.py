import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import R2D2
import sys
import os

try:
    caseid
except NameError:
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


plt.rcParams['font.size'] = 16
# read time
t0 = d.read_time(0,silent=True)

yran = ymax - ymin
xran = min(xmax-xmin,yran)

xsize = 9
ysize = xsize*(yran + xran)/2/yran
fig = plt.figure(num=1,figsize=(xsize,ysize))

grid = GridSpec(2,2,height_ratios=[yran,xran])


for n in range(n0,nd+1):
#for n in range(0,1):
    print(n)
    ##############################
    # read time
    t = d.read_time(n,silent=True)
        
    ##############################
    # read time
    d.read_qq_tau(n*int(ifac),silent=True)

    ##############################
    # read value

    d.read_vc(n,silent=True)
    ##############################

    shading = "flat"
    #shading = "groroud"

    lfac = 1.e-8
    
    ax1 = fig.add_subplot(grid[0,0],aspect='equal')
    ax2 = fig.add_subplot(grid[0,1],aspect='equal')
    ax3 = fig.add_subplot(grid[1,0],aspect='equal')
    ax4 = fig.add_subplot(grid[1,1],aspect='equal')
    
    ax1.tick_params(labelbottom=False)
    if deep_flag == 1:
        d.read_qq_select(xmax,n,silent=True)
        in0 = d.qs["vx"].copy()
        in0s = np.roll(in0,[jx//2-jc,kx//2-kc],axis=[0,1])
        ax1.pcolormesh(y*lfac,z*lfac,in0s.transpose(),cmap='gist_gray',vmax=1.e4,vmin=-1.e4,shading=shading)
        ax1.set_ylabel("z [Mm]")
        ax1.set_title("Vertical velocity")
    else:
        in0 = d.qt["in"].copy()
        in0s = np.roll(in0,[jx//2-jc,kx//2-kc],axis=[0,1])
        ax1.pcolormesh(y*lfac,z*lfac,in0s.transpose(),cmap='gist_gray',vmax=3.2e10,vmin=1.e10,shading=shading)
        ax1.set_ylabel("z [Mm]")
        ax1.set_title("Emergent intensity")

    bx = np.roll(d.qt["bx"],[jx//2-jc,kx//2-kc],axis=[0,1])
    by = np.roll(d.qt["by"],[jx//2-jc,kx//2-kc],axis=[0,1])
    bz = np.roll(d.qt["bz"],[jx//2-jc,kx//2-kc],axis=[0,1])
    ax2.tick_params(labelbottom=False)
    ax2.tick_params(labelleft=False)
    ax2.pcolormesh(y*lfac,z*lfac,bx.transpose(),cmap='gist_gray',vmax=2.5e3,vmin=-2.5e3,shading=shading)
    ax2.set_title(r"LOS magnetic field@$\tau=1$")

    ses = np.roll((d.vc['sep']-d.vc['sem'])/d.vc['serms'],jx//2-jc,axis=1)
    ax3.pcolormesh(y*lfac,(x-rsun)*lfac,ses,vmin=-3.,vmax=3.,cmap='gist_heat',shading=shading)
    tus = np.roll(d.vc["tup"],[jx//2-jc],axis=1)
    ax3.contour(y*lfac,(x-rsun)*lfac,tus,levels=[1.],colors="w")
    ax3.set_ylim((max(xmax-yran,xmin)-rsun)*lfac,(xmax-rsun)*lfac)
    ax3.set_ylabel("x [Mm]")
    ax3.set_xlabel("y [Mm]")
    ax3.set_title(r"$T$")
    
    bb = np.sqrt(d.vc["bxm"]**2 + d.vc["bym"]**2 + d.vc["bzm"]**2)
    bbs = np.roll(bb,[jx//2-jc],axis=1)
    ax4.tick_params(labelleft=False)
    ax4.pcolormesh(y*lfac,(x-rsun)*lfac,bbs,vmax=2.e3,vmin=0.,cmap='gist_heat',shading=shading)
    ax4.contour(y*lfac,(x-rsun)*lfac,tus,levels=[1.],colors="w")
    ax4.set_ylim((max(xmax-yran,xmin)-rsun)*lfac,(xmax-rsun)*lfac)
    ax4.set_xlabel("y [Mm]")
    ax4.set_title(r"$|B|$")

    bbox_props = dict(boxstyle="round,pad=0.3", fc="white", ec="black", lw=2,alpha=0.9)
    ax3.annotate(s="t="+"{:.2f}".format((t-t0)/60/60)+" [hour]"\
                     ,xy=[0.05,0.05],xycoords="figure fraction"\
                     ,fontsize=18,color='black',bbox=bbox_props)

    if(n == n0):
        fig.tight_layout()
    
    plt.pause(0.1)
    plt.savefig(pngdir+"py"+'{0:08d}'.format(n)+".png")

    if(n != nd):
        clf()
