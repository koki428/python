import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import R2D2
import config as c
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
pngdir="../figs/"+caseid+"/png/"
os.makedirs(pngdir,exist_ok=True)

R2D2.read_init(dir,"3d")
for key in c.p:
    exec('%s = %s%s%s' % (key, 'c.p["',key,'"]'))

try:
    n0
except NameError:
    n0 = 0
if  n0 > c.p["ni"]:
    n0 = c.p["ni"]

print("Maximum time step= ",nd," time ="\
      ,dtout*float(nd)/3600./24.," [day]")

plt.close('all')

xsize = 12
ysize = 6
fig = plt.figure(num=1,figsize=(xsize,ysize))

# read time
f = open(dir+"time/t.dac."+'{0:08d}'.format(0),"rb")
t0 = np.fromfile(f,endian+'d',1)
f.close()    
t0 = np.reshape(t0,(1),order="F")

plt.rcParams["font.size"] = 15
#for n in range(193,194):
t_flag = 0
for n in range(n0,ni+1):
#for n in range(n0,200):
    print(n)
    ##############################
    # read time
    f = open(dir+"time/ti.dac."+'{0:08d}'.format(n),"rb")
    t = np.fromfile(f,endian+'d',1)
    f.close()    
    t = np.reshape(t,(1),order="F")

    ##############################
    # read time
    qq_in = R2D2.read_tau_one(dir,n)

    shading = "flat"
    #shading = "groroud"
    ax1 = fig.add_subplot(121,aspect="equal")
    ax2 = fig.add_subplot(122,aspect="equal")
    
    in0 = np.roll(qq_in["in"],[jx//2-jc,kx//2-kc],axis=(0,1))
    bx0 = np.roll(qq_in["bx"],[jx//2-jc,kx//2-kc],axis=(0,1))

    lfac = 1.e-8
    ax1.pcolormesh(y*lfac,z*lfac,in0.transpose(),cmap='gist_gray',vmax=2.5e10,vmin=0.2e10,shading=shading)
    ax2.pcolormesh(y*lfac,z*lfac,bx0.transpose(),cmap='gist_gray',vmax=3e3,vmin=-3.e3,shading=shading)
    ax1.set_xlabel("Mm")
    ax1.set_ylabel("Mm")

    ax2.set_xlabel("Mm")
    ax2.tick_params(labelleft="off",left="off")

    bbox_props = dict(boxstyle="round,pad=0.3", fc="white", ec="black", lw=2,alpha=0.9)
    ax1.annotate(s="t="+"{:.2f}".format((t[0]-t0[0])/3600.)+" [hr]"\
                     ,xy=[0.03,0.03],xycoords="axes fraction"\
                     ,fontsize=18,color='black',bbox=bbox_props)

    if(t_flag == 0):
        t_flag = 1
        fig.tight_layout()


    plt.savefig(pngdir+"py"+'{0:08d}'.format(n)+".png")
    plt.pause(0.1)

    if(n != ni):
        clf()
