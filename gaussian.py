##ガウシアンフィルターをかけるプログラム
import numpy as np
import matplotlib.pyplot as plt
import R2D2
import sys
import os
from scipy.ndimage import gaussian_filter
from mpl_toolkits.axes_grid1 import make_axes_locatable

print("input caseid id")
caseid=0
caseid=input()
caseid='d'+caseid.zfill(3)

datadir="../run/"+caseid+"/data/"
pngdir="../figs/"+caseid+"/png/"
os.makedirs(pngdir,exist_ok=True)

d=R2D2.R2D2_data(datadir)
for key in d.p:
    exec('%s=%s%s%s' % (key, 'd.p["',key,'"]'))

print("Maximum time step= ",nd," time ="\
          ,dtout*float(nd)/3600./24.," [day]")

plt.close('all')

print('input max step number')
nd=input()
nd=int(nd)

t0 = d.read_time(0,silent=True)
fig= plt.figure(figsize=(16, 6))
pr=np.zeros((ix,jx,nd+1))

gg=50
gf=3
for n in range(0,nd+1):
    t=d.read_time(n,silent=True)
    d.read_qq_2d(n,silent=True)
    ro=d.q2['ro']
    se=d.q2['se']
    DprDro,tmp=np.meshgrid(dprdro,y,indexing='ij')
    DprDse,tmp=np.meshgrid(dprdse,y,indexing='ij')
    pr[:,:,n]=DprDro*ro+DprDse*se

pr_gauss=np.zeros((ix,jx,nd+1))
gauss=np.zeros((ix,jx,gg+1))
for n in range(0,nd+1):
    #print(n)
    gauss[:,:,0:gg]=gaussian_filter(pr[:,:,n:n+gg],sigma=[0,0,gf])
    pr_gauss[:,:,n+int(gg/2)]=gauss[:,:,int(gg/2)]
    if (n+gg+1 > nd):
        break
for n in range(0,int(gg/2)):
    pr_gauss[:,:,n]=pr[:,:,n]
    pr_gauss[:,:,nd-n]=pr[:,:,nd-n]

# 結果のプロット
for n in range(0,nd+1):
    ax1=fig.add_subplot(1,2,1)
    p=ax1.pcolormesh(y*1.e-8,(x-rsun)*1.e-8,pr[:,:,n])
    divider=make_axes_locatable(ax1)
    ax_cb=divider.new_horizontal(size="5%", pad=0.1)
    fig.add_axes(ax_cb)
    fig.colorbar(p,cax=ax_cb,label='$[g/cm/s^2]$')
    ax1.set_title('Original')
    s="n = "+str(n)+", t = "+str(round(n/24,2))+" (days)"
    ax1.text(-15,-10,s,fontsize=20)

    ax2=fig.add_subplot(1,2,2)
    p=ax2.pcolormesh(y*1.e-8,(x-rsun)*1.e-8,pr_gauss[:,:,n])
    divider=make_axes_locatable(ax2)
    ax_cb=divider.new_horizontal(size="5%", pad=0.1)
    fig.add_axes(ax_cb)
    fig.colorbar(p,cax=ax_cb,label='$[g/cm/s^2]$')
    ax2.set_title('Gaussian Filter')

    if (n == 0):
        fig.tight_layout(pad=1)

    plt.pause(0.1)
    plt.savefig(pngdir+"gaussian_sigma_"+str(gf)+"_"+str(gg+1)+'{0:08d}'.format(n)+".png")

    if (n != nd):
        clf()


