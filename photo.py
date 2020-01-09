
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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

dir="../run/"+caseid+"/data/"
pngdir="../figs/"+caseid+"/photo/"
os.makedirs(pngdir,exist_ok=True)

d = R2D2.R2D2_data(dir)
for key in d.p:
    exec('%s = %s%s%s' % (key, 'd.p["',key,'"]'))

try:
    n0
except NameError:
    n0 = 0
if  n0 > d.p["ni"]:
    n0 = d.p["ni"]

print("Maximum time step= ",ni," time ="\
      ,dtout/ifac*float(nd)/3600./24.," [day]")

plt.close('all')

xsize = 15
ysize = 7.5
fig = plt.figure(num=1,figsize=(xsize,ysize))

t0 = d.read_time(0,tau=True,silent=True)

for n in range(n0,ni+1):
#for n in range(20,21):
    print(n)
    ##############################
    # read time
    t = d.read_time(n,tau=True,silent=True)

    ##############################
    # read time
    d.read_qq_tau(n,silent=True)

    shading = "flat"
    #shading = "groroud"
    ax1 = fig.add_subplot(121,aspect="equal")
    ax2 = fig.add_subplot(122,aspect="equal")
    
    in0 = np.roll(d.qt["in"],[jx//2-jc,kx//2-kc],axis=(0,1))
    bx0 = np.roll(d.qt["bx"],[jx//2-jc,kx//2-kc],axis=(0,1))

    lfac = 1.e-8
    ax1.pcolormesh(y*lfac,z*lfac,in0.transpose(),cmap='gist_gray',vmax=3.0e10,vmin=0.2e10,shading=shading)
    ax2.pcolormesh(y*lfac,z*lfac,bx0.transpose(),cmap='gist_gray',vmax=5.e2,vmin=-5.e2,shading=shading)
    ax1.set_xlabel("Mm")
    ax1.set_ylabel("Mm")

    ax2.set_xlabel("Mm")
    ax2.tick_params(labelleft=False,left=False)

    bbox_props = dict(boxstyle="round,pad=0.3", fc="white", ec="black", lw=2,alpha=0.9)
    ax1.annotate(s="t="+"{:.2f}".format((t-t0)/3600.)+" [hr]"\
                     ,xy=[0.03,0.03],xycoords="axes fraction"\
                     ,fontsize=18,color='black',bbox=bbox_props)

    if n == n0:
        fig.tight_layout()


    plt.savefig(pngdir+"py"+'{0:08d}'.format(n)+".png")
    plt.pause(0.1)

    if(n != ni):
        clf()
