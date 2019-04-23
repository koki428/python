import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import read
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
figdir="../figs/"
casedir="../figs/"+caseid
pngdir="../figs/"+caseid+"/png/"
if not os.path.exists(figdir):
    os.mkdir(figdir)
if not os.path.exists(casedir):
    os.mkdir(casedir)
if not os.path.exists(pngdir):
    os.mkdir(pngdir)

read.read_init(dir,"3d")
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

plt.close('all')

xsize = 10
ysize = 10
fig = plt.figure(num=1,figsize=(xsize,ysize))

# read time
f = open(dir+"time/t.dac."+'{0:08d}'.format(0),"rb")
t0 = np.fromfile(f,endian+'d',1)
f.close()    
t0 = np.reshape(t0,(1),order="F")

plt.rcParams["font.size"] = 15

#n0 = 30
#nd = n0
    
#for n in range(40,41):

#jc = jx//2
#kc = kx//2
for n in range(n0,nd+1):
#for n in range(0,1):
#for n in range(n0,200):
    print(n)
    ##############################
    # read time
    f = open(dir+"time/t.dac."+'{0:08d}'.format(n),"rb")
    t = np.fromfile(f,endian+'d',1)
    f.close()    
    t = np.reshape(t,(1),order="F")

    ##############################
    # read time
    qq_in = read.read_tau_one(dir,n*int(ifac))

    ##############################
    # read value
    f = open(dir+"remap/vla.dac."+'{0:08d}'.format(n),"rb")
    vl0 = np.fromfile(f,c.p["endian"]+'f',m2da*ix*jx)
    f.close()

    vl = np.reshape(vl0,(ix,jx,m2da),order="F")

    vc = {"a":0}
    for m in range(c.p["m2da"]):
        vc[c.p["cl"][m]] = vl[:,:,m]
    ##############################

    shading = "flat"
    #shading = "groroud"

    lfac = 1.e-8
    
    ax1 = fig.add_subplot(221,aspect="equal")
    ax2 = fig.add_subplot(222,aspect="equal")
    ax3 = fig.add_subplot(223,aspect="equal")
    ax4 = fig.add_subplot(224,aspect="equal")

    ax1.tick_params(labelbottom=False)
    in0 = qq_in["in"].copy()
    in0s = np.roll(in0,[jx//2-jc,kx//2-kc],axis=[0,1])
    ax1.pcolormesh(y*lfac,z*lfac,in0s.transpose(),cmap='gist_gray',vmax=3.2e10,vmin=1.e10,shading=shading)
    ax1.set_ylabel("z [$R_\odot$]")
    ax1.set_title("Emergent intensity")

    bx = np.roll(qq_in["bx"],[jx//2-jc,kx//2-kc],axis=[0,1])
    by = np.roll(qq_in["by"],[jx//2-jc,kx//2-kc],axis=[0,1])
    bz = np.roll(qq_in["bz"],[jx//2-jc,kx//2-kc],axis=[0,1])
    #bx = qq_in[8,:,:]
    ax2.tick_params(labelbottom=False)#,bottom="off")
    ax2.tick_params(labelleft=False)#,left="off")
    ax2.pcolormesh(y*lfac,z*lfac,bx.transpose(),cmap='gist_gray',vmax=2.5e3,vmin=-2.5e3,shading=shading)
    ax2.set_title(r"LOS magnetic field@$\tau=1$")

    ses = np.roll((vc["sep"]-vc["sem"])/vc["serms"],jx//2-jc,axis=1)
    tus = np.roll(vc["tup"],[jx//2-jc],axis=1)
    ax3.pcolormesh(y*lfac,(x-rsun)*lfac,ses,vmax=3.,vmin=-3.,cmap='gist_heat',shading=shading)
    ax3.contour(y*lfac,(x-rsun)*lfac,tus,levels=[1.],colors="w")
    ax3.set_ylim((xmax-rsun-(ymax-ymin))*lfac,(xmax-rsun)*lfac)
    #ax3.pcolormesh((vc["sep"]-vc["sem"])/vc["serms"],vmax=3.,vmin=-3.,cmap='gist_heat',shading=shading)
    ax3.set_xlabel("y [$R_\odot$]")
    ax3.set_title(r"$(s-\langle s\rangle)/s_{rms}$")

    bb = np.sqrt(vc["bxp"]**2 + vc["byp"]**2 + vc["bzp"]**2)

    bbs = np.roll(bb,[jx//2-jc],axis=1)
    ax4.tick_params(labelleft=False)
    ax4.pcolormesh(y*lfac,(x-rsun)*lfac,bbs,vmax=1.e3,vmin=0.,cmap='gist_heat',shading=shading)
    ax4.contour(y*lfac,(x-rsun)*lfac,tus,levels=[1.],colors="w")
    ax4.set_ylim((xmax-rsun-(ymax-ymin))*lfac,(xmax-rsun)*lfac)
    ax4.set_xlabel("y [$R_\odot$]")
    ax4.set_ylabel("x [$R_\odot$]")
    ax4.set_title(r"$\langle B^2\rangle$")
    #ax4.pcolormesh(bb,vmax=1.e4,vmin=0.,cmap='gist_heat',shading=shading)

    bbox_props = dict(boxstyle="round,pad=0.3", fc="white", ec="black", lw=2,alpha=0.9)
    ax3.annotate(s="t="+"{:.2f}".format((t[0]-t0[0])/3600.)+" [hr]"\
                     ,xy=[0.02,0.02],xycoords="figure fraction"\
                     ,fontsize=18,color='black',bbox=bbox_props)

    if(n == n0):
        plt.tight_layout(pad=0.1)

    plt.pause(0.1)
    plt.savefig(pngdir+"py"+'{0:08d}'.format(n)+".png")

    if(n != nd):
        clf()
