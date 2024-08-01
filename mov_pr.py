#圧力の画像を作成

import numpy as np
import matplotlib.pyplot as plt
import R2D2
import sys
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable

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

fig = plt.figure(figsize=(12,6))

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
    ro=d.q2['ro']
    se=d.q2['se']
    DprDro,tmp=np.meshgrid(dprdro,y,indexing='ij')
    DprDse,tmp=np.meshgrid(dprdse,y,indexing='ij')
    
    ax1 = fig.add_subplot(111,aspect='equal')
    
    pr=DprDro*ro+DprDse*se
    p=ax1.pcolormesh(y*1.e-8,(x-rsun)*1.e-8,pr)
    divider=make_axes_locatable(ax1)
    ax_cb=divider.new_horizontal(size="5%", pad=0.1)
    fig.add_axes(ax_cb)
    fig.colorbar(p,cax=ax_cb,label='$[g/cm/s^2]$')
    s="n = "+str(n)+", t = "+str(round(n/24,2))+" (days)"
    ax1.text(-15,0,s,fontsize=20)

    #ax1.contour(y*1.e-8,(x-rsun)*1.e-8,d.q2['tu'],levels=[1.])

    
    
    #ax1.annotate(s="t="+"{:.2f}".format((t-t0)/60/60)+" [hour]"\
    #                ,xy=[0.05,0.03],xycoords="figure fraction"\
    #                ,fontsize=18,color='black')#,bbox=bbox_props)
    
    if(n == n0):
        fig.tight_layout(pad=0.1)
    
    plt.pause(0.5)
    plt.savefig(pngdir+"py_pr"+'{0:08d}'.format(n)+".png")

    if(n != nd):
        clf()
