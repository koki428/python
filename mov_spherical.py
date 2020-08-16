import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from mpl_toolkits.basemap import Basemap
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

X, Y = np.meshgrid(x,y,indexing='ij')
SINY = sin(Y)
SINYM = sin(Y).sum(axis=1)

plt.rcParams['font.size'] = 16
# read time
t0 = d.read_time(0,silent=True)

yran = ymax - ymin
xran = min(xmax-xmin,yran)

xsize = 12
ysize = 12
fig = plt.figure(num=1,figsize=(xsize,ysize))

#grid = GridSpec(2,2,height_ratios=[yran,xran])

RA, TH = np.meshgrid(x,y,indexing='ij')
XX, YY = RA*cos(TH), RA*sin(TH)

yc = y/pi*180
zc = z/pi*180


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

    ax1 = fig.add_subplot(221,aspect='equal')
    ax2 = fig.add_subplot(222,aspect='equal')
    ax3 = fig.add_subplot(223,aspect='equal')
    ax4 = fig.add_subplot(224,aspect='equal')

    d.read_qq_select(xmax,n,silent=True)
    vx = d.qs['vx']
    bx = d.qs['bx']
    ax1.pcolormesh(yc,zc,vx)
    ax2.pcolormesh(yc,zc,bx)

    sem, tmp   = np.meshgrid((d.vc['sem']*SINY).sum(axis=1)/SINYM,y,indexing='ij')
    serms, tmp = np.meshgrid(sqrt((d.vc['serms']**2*SINY).sum(axis=1)/SINYM),y,indexing='ij')
    bbp = sqrt(d.vc['bxp']**2 + d.vc['byp']**2 + d.vc['bzp']**2)
    
    lfac = 1/rsun
    ax3.pcolormesh(XX.T*lfac,YY.T*lfac,((d.vc['sep']-sem)/serms).T,vmin=-2.,vmax=2.)
    ax4.pcolormesh(XX.T*lfac,YY.T*lfac,bbp.T)
    
    if(n == n0):
        fig.tight_layout()
    
    plt.pause(0.1)
    plt.savefig(pngdir+"py"+'{0:08d}'.format(n)+".png")

    if(n != nd):
        clf()
