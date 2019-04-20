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

read.read_init(dir,"debug")
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
ysize = 5
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

    shading = "flat"
    shading = "groroud"
    
    ax1 = fig.add_subplot(121,aspect="equal")
    ax2 = fig.add_subplot(122,aspect="equal")

    ax1.tick_params(labelbottom=False)#,bottom="off")
    #ax1.tick_params(labelleft="off",left="off")
    in0 = qq_in["in"].copy()
    in0s = np.roll(in0,[jx//2-jc,kx//2-kc],axis=[0,1])
    ax1.pcolormesh(y/rsun,z/rsun,in0s.transpose(),cmap='gist_gray',vmax=3.2e10,vmin=1.e10,shading=shading)
    ax1.set_ylabel("z [$R_\odot$]")
    ax1.set_title("Emergent intensity")

    bx = np.roll(qq_in["bx"],[jx//2-jc,kx//2-kc],axis=[0,1])
    by = np.roll(qq_in["by"],[jx//2-jc,kx//2-kc],axis=[0,1])
    bz = np.roll(qq_in["bz"],[jx//2-jc,kx//2-kc],axis=[0,1])
    #bx = qq_in[8,:,:]
    ax2.tick_params(labelbottom=False)#,bottom="off")
    ax2.tick_params(labelleft=False)#,left="off")
    ax2.pcolormesh(y/rsun,z/rsun,bx.transpose(),cmap='gist_gray',vmax=2.5e3,vmin=-2.5e3,shading=shading)
    ax2.set_title(r"LOS magnetic field@$\tau=1$")


    if(n == n0):
        plt.tight_layout(pad=0.1)

    plt.pause(0.1)

    if(n != nd):
        clf()
