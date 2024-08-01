#磁場とエントロピーを横に並べた画像を作成
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import R2D2
import sys
import os

#try:
#    caseid
#except NameError:
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

# read time
t0 = d.read_time(0,silent=True)

yran = ymax - ymin
xran = min(xmax-xmin,yran)

xsize = 12
ysize = xsize*xran/yran/2
fig = plt.figure(num=1,figsize=(10,4))

for n in range(n0,nd+1):
#for n in range(0,1):
    print(n)
    ##############################
    # read time
    t = d.read_time(n,silent=True)
        

    ##############################
    # read value

    d.read_qq_2d(n,silent=True)
    ##############################

    #shading = "flat"
    shading = "gouraud"

    lfac = 1.e-8
    plt.rcParams['font.size']=14
    
    ax1 = fig.add_subplot(111,aspect='equal')
    # ax2 = fig.add_subplot(122,aspect='equal')

    sem = d.q2['se'].mean(axis=1)
    sem2, tmp = np.meshgrid(sem,y,indexing='ij')
    serms = np.sqrt(((d.q2['se'] - sem2)**2).mean(axis=1))
    serms2, tmp = np.meshgrid(serms,y,indexing='ij')
    
    
    ax1.pcolormesh(y*1.e-8,(x-rsun)*1.e-8,(d.q2['se']-sem2)/serms2,vmin=-4,vmax=4)
    ax1.contour(y*1.e-8,(x-rsun)*1.e-8,d.q2['tu'],levels=[1.])
    ax1.set_xlabel("y [Mm]")
    ax1.set_ylabel("x [Mm]")
    ax1.annotate("t = "+str(n*8)+" [hours]", xy=(0.03,0.90), xycoords="figure fraction",fontsize=18, color='black')

    # bz = np.sqrt(d.q2['bz']**2)

    # ax2.pcolormesh(y*1.e-8,(x-rsun)*1.e-8,bz,cmap='gist_stern',shading=shading)
    # ax2.contour(y*1.e-8,(x-rsun)*1.e-8,d.q2['tu'],levels=[1.])
    
    # ax1.annotate(s="t="+"{:.2f}".format((t-t0)/60/60)+" [hour]"\
    #                  ,xy=[0.05,0.03],xycoords="figure fraction"\
    #                  ,fontsize=18,color='black')#,bbox=bbox_props)

    if(n == n0):
        fig.tight_layout(pad=0.5)
    
    plt.pause(0.1)
    plt.savefig(pngdir+"py_se"+'{0:08d}'.format(n)+".png")

    if(n != nd):
        clf()
