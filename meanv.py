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

d = R2D2.R2D2_data(datadir)
for key in d.p:
    exec('%s = %s%s%s' % (key, 'd.p["',key,'"]'))

try:
    n0
except NameError:
    n0 = 0
if  n0 > d.p["nd"]:
    n0 = d.p["nd"]

n0 = 0
t = np.zeros(nd-n0+1)
ekm = np.zeros(nd-n0+1)
anm = np.zeros(nd-n0+1)
bxmt = np.zeros((ix,jx,nd-n0+1))
bymt = np.zeros((ix,jx,nd-n0+1))
bzmt = np.zeros((ix,jx,nd-n0+1))

RR, TH = np.meshgrid(x,y,indexing='ij')
ro2, tmp = np.meshgrid(ro0,y,indexing='ij')

for n in range(n0,nd+1):
    print(n)

    d.read_time(n)
    t[n-n0] = d.t


    d.read_vc(n,silent=True)
    ekm[n-n0] = (RR**2*sin(TH)*d.vc['vzm']**2).mean()/(RR**2*sin(TH)).mean()
    anm[n-n0] = ((RR*sin(TH)*d.vc['vzm'])*RR**2*sin(TH)*ro2).sum()
    bxmt[:,:,n-n0] = d.vc['bxm']
    bymt[:,:,n-n0] = d.vc['bym']
    bzmt[:,:,n-n0] = d.vc['bzm']

plt.close('all')
plt.clf()
plt.figure(100,figsize=(8,4))
plt.pcolormesh(t/86400/365,90-y/np.pi*180,bzmt[0,:,:],vmin=-8.e3,vmax=8.e3,cmap='bwr',shading='auto',rasterized=True)
plt.xlabel('t [year]')
plt.ylabel('latitude [degree]')
plt.xlim(0,30)
plt.tight_layout()
plt.savefig(caseid+'.pdf')
