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

datadir="../run/"+caseid+"/data/"
casedir="../figs/"+caseid
os.makedirs(casedir,exist_ok=True)

pngdir="../figs/"+caseid+"/divb/"
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

ip = np.argmin(abs(x-rsun))
bx = np.zeros((5,jx,kx))
by = np.zeros((5,jx,kx))
bz = np.zeros((5,jx,kx))
xl = x[ip-2:ip+3]

plt.clf()
fig = plt.figure(100,figsize=(10,5))

#n0 = 0
#nd = n0

print("Maximum time step= ",nd," time ="\
          ,dtout*float(nd)/3600./24.," [day]")

for n in range(n0,nd+1):
    print(n)
    for i in range(ip-2,ip+3):
        d.read_qq_select(x[i],n,silent=True)
        ii = i - (ip-2)
        bx[ii,:,:] = d.qs['bx']
        by[ii,:,:] = d.qs['by']
        bz[ii,:,:] = d.qs['bz']

        divb = (R2D2.derv.d_x(bx,xl) + R2D2.derv.d_y(by,y) + R2D2.derv.d_z(bz,z))[2,:,:]
        
    ax1 = fig.add_subplot(121,aspect='equal')
    ax2 = fig.add_subplot(122,aspect='equal')

    ax1.pcolormesh(divb,vmin=-1.e-7,vmax=1.e-7)
    ax1.set_title(r'$\nabla\cdot B$')

    d.read_qq_select(x[ip],n,silent=True)
    ax2.pcolormesh(d.qs['ph'],vmin=-5.e5,vmax=5.e5)
    ax2.set_title(r'$\Psi$')

    if(n == n0):
        fig.tight_layout()

    plt.pause(0.1)
    plt.savefig(pngdir+'py'+'{0:08d}'.format(n)+'.png')

    if(n != nd):
        plt.clf()
