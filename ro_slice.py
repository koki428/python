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
pngdir="../figs/"+caseid+"/ro_slice/"
os.makedirs(pngdir,exist_ok=True)

d = R2D2.R2D2_data(datadir)
for key in d.p:
    exec('%s = %s%s%s' % (key, 'd.p["',key,'"]'))
    
try:
    n0
except NameError:
    n0 = 0
if  n0 > nd_tau:
    n0 = nd_tau

print("Maximum time step= ",nd_tau," time ="\
          ,dtout*float(nd)/3600./24.," [day]")

plt.close('all')

plt.rcParams['font.size'] = 16
# read time
t0 = d.read_time(0,silent=True)

yran = ymax - ymin
xran = min(xmax-xmin,yran)

xsize = 7
ysize = xsize*(yran + xran)/yran
fig = plt.figure(num=1,figsize=(xsize,ysize))

grid = GridSpec(2,1,height_ratios=[yran,xran])

te2, tmp = np.meshgrid(te0,y,indexing='ij')

#n0 = 0
#nd_tau = n0 + 1
for n in range(n0,nd_tau+1):
#for n in range(0,1):
    print(n)
    ##############################
    # read time
    t = d.read_time(n,tau=True,silent=True)
        
    ##############################

    shading = "flat"
    #shading = "groroud"

    lfac = 1.e-8
    
    ax1 = fig.add_subplot(grid[0,0],aspect='equal')
    ax2 = fig.add_subplot(grid[1,0],aspect='equal')
    
    vmax = 8.e-8

    ax1.tick_params(labelbottom=False)
    d.read_qq_slice(n,2,'x',silent=True)
    ax1.pcolormesh(y*lfac,z*lfac,d.ql['ro'].T-d.ql['ro'].mean(),cmap='gist_gray',shading=shading,vmax=vmax,vmin=-vmax)
    ax1.set_ylabel("z [Mm]")
    ax1.set_title("Density perturbation")


    d.read_qq_slice(n,1,'y',silent=True)
    
    ax2.pcolormesh(y*lfac,(x-rsun)*lfac,d.ql['ro'],cmap='gist_gray',vmax=vmax,vmin=-vmax)
    ax2.set_ylim((max(xmax-yran,xmin)-rsun)*lfac,(xmax-rsun)*lfac)
    ax2.set_ylabel("x [Mm]")
    ax2.set_xlabel("y [Mm]")
    

    ax2.annotate(s="t="+"{:.2f}".format((t-t0)/60/60)+" [hour]"\
                     ,xy=[0.05,0.03],xycoords="figure fraction"\
                     ,fontsize=18,color='black')#,bbox=bbox_props)

    if(n == n0):
        fig.tight_layout()
    
    plt.pause(0.1)
    plt.savefig(pngdir+"py"+'{0:08d}'.format(n)+".png")

    if(n != nd_tau):
        clf()
