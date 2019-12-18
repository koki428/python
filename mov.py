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
if  n0 > c.p["nd"]:
    n0 = c.p["nd"]

print("Maximum time step= ",nd," time ="\
          ,dtout*float(nd)/3600./24.," [day]")

tmp, te2 = np.meshgrid(y,te0)
plt.close('all')

vsize = 8

zfac0 = 1
xfac0 = zfac0*(xmax - xmin)/(zmax - zmin)
marginfac_vbot0 = 0.1
marginfac_vint0 = 0.06
marginfac_vtop0 = 0.06
vsize0 = zfac0 + xfac0 + marginfac_vbot0 + marginfac_vint0 + marginfac_vtop0
xfac = xfac0/vsize0
zfac = zfac0/vsize0
marginfac_vbot = marginfac_vbot0/vsize0
marginfac_vint = marginfac_vint0/vsize0
marginfac_vtop = marginfac_vtop0/vsize0

marginlen_vbot = vsize*marginfac_vbot
marginlen_vint = vsize*marginfac_vint
xlen = vsize*xfac
zlen = vsize*zfac
ylen = zlen*(ymax - ymin)/(zmax - zmin)
marginlen_hbot = ylen*0.15
marginlen_hint = ylen*0.02
marginlen_htop = ylen*0.02

hsize = 2*ylen + marginlen_hbot + marginlen_hint + marginlen_htop
yfac = ylen/hsize

h0 =  marginlen_hbot/hsize
h1 = (marginlen_hbot + marginlen_hint + ylen)/hsize

v0 =  marginlen_vbot/vsize
v1 = (marginlen_vbot + marginlen_vint + xlen)/vsize

fig = plt.figure(num=1,figsize=(hsize,vsize))

# read time
t0 = R2D2.read_time(dir,0)

plt.rcParams["font.size"] = 15

#n0 = 7
#nd = n0

for n in range(n0,nd+1):
#for n in range(0,1):
    print(n)
    ##############################
    # read time
    t = R2D2.read_time(dir,n)
        
    ##############################
    # read time
    qq_in = R2D2.read_tau_one(dir,n*int(ifac))

    ##############################
    # read value

    vc = R2D2.read_vc(dir,n)
    ##############################

    shading = "flat"
    #shading = "groroud"

    lfac = 1.e-8
    
    ax1 = fig.add_axes([h0,v1,yfac,zfac])
    ax2 = fig.add_axes([h1,v1,yfac,zfac])
    ax3 = fig.add_axes([h0,v0,yfac,xfac])
    ax4 = fig.add_axes([h1,v0,yfac,xfac])

    ax1.tick_params(labelbottom=False)
    in0 = qq_in["in"].copy()
    in0s = np.roll(in0,[jx//2-jc,kx//2-kc],axis=[0,1])
    ax1.pcolormesh(y*lfac,z*lfac,in0s.transpose(),cmap='gist_gray',vmax=3.2e10,vmin=1.e10,shading=shading)
    ax1.set_ylabel("z [Mm]")
    ax1.set_title("Emergent intensity")

    bx = np.roll(qq_in["bx"],[jx//2-jc,kx//2-kc],axis=[0,1])
    by = np.roll(qq_in["by"],[jx//2-jc,kx//2-kc],axis=[0,1])
    bz = np.roll(qq_in["bz"],[jx//2-jc,kx//2-kc],axis=[0,1])
    ax2.tick_params(labelbottom=False)
    ax2.tick_params(labelleft=False)
    ax2.pcolormesh(y*lfac,z*lfac,bx.transpose(),cmap='gist_gray',vmax=2.5e3,vmin=-2.5e3,shading=shading)
    ax2.set_title(r"LOS magnetic field@$\tau=1$")

    #ses = np.roll(vc['tep']+te2,jx//2-jc,axis=1)
    ses = np.roll((vc['sep']-vc['sem'])/vc['serms'],jx//2-jc,axis=1)
    #ax3.pcolormesh(y*lfac,(x-rsun)*lfac,ses,vmin=3000.,vmax=18000.,cmap='gist_heat',shading=shading)
    ax3.pcolormesh(y*lfac,(x-rsun)*lfac,ses,vmin=-3.,vmax=3.,cmap='gist_heat',shading=shading)
    tus = np.roll(vc["tup"],[jx//2-jc],axis=1)
    ax3.contour(y*lfac,(x-rsun)*lfac,tus,levels=[1.],colors="w")
    ax3.set_ylabel("x [Mm]")
    ax3.set_xlabel("y [Mm]")
    ax3.set_title(r"$T$")
    
    bb = np.sqrt(vc["bxm"]**2 + vc["bym"]**2 + vc["bzm"]**2)
    bbs = np.roll(bb,[jx//2-jc],axis=1)
    ax4.tick_params(labelleft=False)
    ax4.pcolormesh(y*lfac,(x-rsun)*lfac,bbs,vmax=1.e3,vmin=0.,cmap='gist_heat',shading=shading)
    ax4.contour(y*lfac,(x-rsun)*lfac,tus,levels=[1.],colors="w")
    ax4.set_xlabel("y [Mm]")
    ax4.set_title(r"$|B|$")

    bbox_props = dict(boxstyle="round,pad=0.3", fc="white", ec="black", lw=2,alpha=0.9)
    ax3.annotate(s="t="+"{:.2f}".format((t-t0)/60/60)+" [hour]"\
                     ,xy=[0.02,0.02],xycoords="figure fraction"\
                     ,fontsize=18,color='black',bbox=bbox_props)
        
    plt.pause(0.1)
    plt.savefig(pngdir+"py"+'{0:08d}'.format(n)+".png")

    if(n != nd):
        clf()
