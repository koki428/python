import numpy as np
import math
import sympy as sym
import matplotlib.pyplot as plt
import R2D2
import sys
import os

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

try:
    n0
except NameError:
    n0=0
if n0>d.p["nd"]:
    n0=d.p["nd"]

print("Maximum time step=",nd,"time="\
        ,dtout*float(nd)/3600./24.,"[day]")

plt.close("all")

t0=d.read_time(0,silent=True)

dx=(xmax-xmin)/ix
dy=(ymax-ymin)/jx

print('input step number')
n=input()
n=int(n)

t=d.read_time(n,silent=True)
d.read_qq_2d(n,silent=True)

bx=d.q2['bx']
by=d.q2['by']
bz=d.q2['bz']

#psi cal
psi_x=np.zeros((ix,jx))
psi_y=np.zeros((ix,jx))
psi=np.zeros((ix,jx))

psi_x[0,0]=0.0
j=0
for i in range(1,ix):
    psi_x[i,j]=psi_x[i-1,j]+by[i,j]*dx
for j in range(1,jx):
    i=0
    psi_x[i,j]=psi_x[i,j-1]-bx[i,j]*dy
    for i in range(1,ix):
        psi_x[i,j]=psi_x[i-1,j]+by[i,j]*dx

psi_y[0,0]=0.0
i=0
for j in range(1,jx):
    psi_y[i,j]=psi_y[i,j-1]-bx[i,j]*dy
for i in range(1,ix):
    j=0
    psi_y[i,j]=psi_y[i-1,j]+by[i,j]*dx
    for j in range(1,jx):
        psi_y[i,j]=psi_y[i,j-1]-bx[i,j]*dy

#for i in range(0,ix):
#    for j in range(0,jx):
#        if (j==0) and (i==0):
#            psi_y4[i,j]=0.0
#        elif (j==0) and (i>=1):
#            psi_y4[i,j]=psi_y4[i-1,j]+by[i,j]*dx
#        elif (j>=1):
#            psi_y4[i,j]=psi_y4[i,j-1]-bx[i,j]*dy

for j in range(0,jx):
    for i in range(0,ix):
        psi[i,j]=(psi_x[i,j]+psi_y[i,j])*0.5

#print(np.max(psi))

#main tube b cal
#tube=0
#for j in range(0,jx):
#    for i in range(0,ix):
#       if (psi[i,j] >= np.max(psi)*0.2):
#            tube=tube+bz[i,j]
#print(tube)



fig=plt.figure(figsize=(10,12))
lfac=1.e-8

ax4=fig.add_subplot(111,aspect='equal')
ax4=ax4.pcolormesh(y*lfac,(x-rsun)*lfac,psi,vmin=0,cmap='gist_stern')
plt.axis('off')

plt.savefig(pngdir+'psi_simple_'+'{0:08d}'.format(n)+'.png')
#ax44=plt.contour(y*lfac,(x-rsun)*lfac,psi,colors='g',levels=[0.5e+12,1.e+12,1.5e+12,2.e+12,2.5e+12,3.e+12,3.5e+12])
#ax44.clabel(fmt='%1.1f',fontsize=6,colors='white')
#fig.colorbar(ax4)

