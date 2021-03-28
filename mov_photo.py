import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from tqdm import tqdm
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
pngdir="../figs/"+caseid+"/mov_photo/"
os.makedirs(pngdir,exist_ok=True)

d = R2D2.R2D2_data(datadir)
for key in d.p:
    exec('%s = %s%s%s' % (key, 'd.p["',key,'"]'))
    
try:
    n0
except NameError:
    n0 = 0
if  n0 > d.p["nd_tau"]:
    n0 = d.p["nd_tau"]

print("Maximum time step (nd_tau) = ",nd_tau," time ="\
          ,dtout*float(nd)/3600./24.," [day]")

plt.close('all')

# read initial time
t0 = d.read_time(0,silent=True)

yran = ymax - ymin
xran = min(xmax-xmin,yran)

xsize = 12
ysize = xsize*(yran + xran)/2/yran
fig = plt.figure(num=1,figsize=(xsize,ysize))

grid = GridSpec(2,2,height_ratios=[yran,xran])

te2, tmp = np.meshgrid(te0,y,indexing='ij')

#n0 = 18
#nd_tau = n0

d.read_vc(0,silent=True)

for n in tqdm(range(n0,nd_tau+1)):
#for n in range(0,1):
    #print(n)
    ##############################
    # read time
    t = d.read_time(n,tau=True,silent=True)
        
    ##############################
    # read time
    if xmax > rsun:
        d.read_qq_tau(n,silent=True)

    ##############################
    # read value

    k_slice = 1
    d.read_qq_slice(k_slice,'z',n,silent=True)
    kl = np.argmin(np.abs(z - z_slice[k_slice]))
    ##############################

    shading = "auto"

    lfac = 1.e-8
    
    ax1 = fig.add_subplot(grid[0,0],aspect='equal')
    ax2 = fig.add_subplot(grid[0,1],aspect='equal')
    ax3 = fig.add_subplot(grid[1,0],aspect='equal')
    ax4 = fig.add_subplot(grid[1,1],aspect='equal')
    
    ax1.tick_params(labelbottom=False)
    in0 = d.qt["in"].copy()
    in0s = np.roll(in0,[jx//2-jc,kx//2-kc],axis=[0,1])
    ax1.pcolormesh(y*lfac,z*lfac,in0s.transpose(),cmap='gist_gray',vmax=3.2e10,vmin=1.e10,shading=shading)
    ax1.set_ylabel("Mm")
    ax1.set_title("Emergent intensity")
    ax3.pcolormesh(y*lfac,(x-rsun)*lfac,d.ql['te']+te2,vmin=3000,vmax=20000,cmap='gist_heat',shading=shading)
    ax3.plot(y*lfac,(d.qt['he'][:,kl]-rsun)*lfac,color='w')
    ax3.set_title(r"$T$")
    ax4.plot(y*lfac,(d.qt['he'][:,kl]-rsun)*lfac,color='w')

    bx = np.roll(d.qt["bx"],[jx//2-jc,kx//2-kc],axis=[0,1])
    ax2.tick_params(labelbottom=False)
    ax2.tick_params(labelleft=False)
    ax2.pcolormesh(y*lfac,z*lfac,bx.transpose(),cmap='gist_gray',vmax=0.25e3,vmin=-0.25e3,shading=shading)
    ax2.set_title(r"LOS magnetic field@$\tau=1$")

    ax3.set_ylim((max(xmax-yran,xmin)-rsun)*lfac,(xmax-rsun)*lfac)
    ax3.set_ylabel("Mm")
    ax3.set_xlabel("Mm")
    
    bb = np.sqrt(d.ql["bx"]**2 + d.ql["by"]**2 + d.ql["bz"]**2)
    bbs = np.roll(bb,[jx//2-jc],axis=1)
    ax4.tick_params(labelleft=False)
    ax4.pcolormesh(y*lfac,(x-rsun)*lfac,bbs,vmax=2.e3,vmin=0.,cmap='gist_heat',shading=shading)
    ax4.set_ylim((max(xmax-yran,xmin)-rsun)*lfac,(xmax-rsun)*lfac)
    ax4.set_xlabel("$y$ [Mm]")
    ax4.set_title(r"$|B|$")

    ax3.annotate(text="$t$="+"{:.2f}".format((t-t0)/60/60)+" [hour]"\
                     ,xy=[0.03,0.03],xycoords="figure fraction"\
                     ,color='black')#,bbox=bbox_props)

    if(n == n0):
        fig.tight_layout(pad=0.1)
    
    plt.pause(0.1)
    plt.savefig(pngdir+"py"+'{0:08d}'.format(n)+".png")

    if(n != nd_tau):
        clf()
