import cv2
import numpy as np
import R2D2
import matplotlib.pyplot as plt
import math
from decimal import Decimal, ROUND_HALF_UP, ROUND_HALF_EVEN

def psi_compare(xmax,xmin,ymax,ymin):
    #print("input caseid")
    #caseid=0
    #caseid=input()
    #caseid="d"+caseid.zfill(3)

    #datadir="../run/"+caseid+"/data/"
    #pngdir="../figs/"+caseid+"/png/"
    #os.makedirs(pngdir,exist_ok=True)

    #d=R2D2.R2D2_data(datadir)
    #for key in d.p:
    #    exec('%s=%s%s%s' % (key, 'd.p["',key,'"]'))

    #try:
    #    n0
    #except NameError:
    #    n0=0
    #if n0>d.p["nd"]:
    #    n0=d.p["nd"]

   # print("Maximum time step=",nd,"time="\
   #         ,dtout*float(nd)/3600./24.,"[day]")

    #plt.close("all")

    t0=d.read_time(0,silent=True)

    dx=(xmax-xmin)/jx
    dy=(ymax-ymin)/ix

    print('input step number')
    n=input()
    n=int(n)

    t=d.read_time(n,silent=True)
    d.read_qq_2d(n,silent=True)

    bx=d.q2['bx']
    by=d.q2['by']
    #bz=d.q2['bz']

    #psi cal
    #psi_x=np.zeros((ix,jx))
    #psi_y=np.zeros((ix,jx))
    psi=np.zeros((ix,jx))

    #psi_x[0,0]=0.0
    #j=0
    #for i in range(1,ix):
    #    psi_x[i,j]=psi_x[i-1,j]+by[i,j]*dx
    #for j in range(1,jx):
    #    i=0
    #    psi_x[i,j]=psi_x[i,j-1]-bx[i,j]*dy
    #    for i in range(1,ix):
    #        psi_x[i,j]=psi_x[i-1,j]+by[i,j]*dx

    #psi_y[0,0]=0.0
    #i=0
    #for j in range(1,jx):
    #    psi_y[i,j]=psi_y[i,j-1]-bx[i,j]*dy
    #for i in range(1,ix):
    #    j=0
    #    psi_y[i,j]=psi_y[i-1,j]+by[i,j]*dx
    #    for j in range(1,jx):
    #        psi_y[i,j]=psi_y[i,j-1]-bx[i,j]*dy

    psi[0,0]=0.0
    i=0
    for j in range(1,jx):
        psi[i,j]=psi[i,j-1]-bx[i,j]*dy
    for i in range(1,ix):
        j=0
        psi[i,j]=psi[i-1,j]+by[i,j]*dx
        for j in range(1,jx):
            psi[i,j]=psi[i,j-1]-bx[i,j]*dy

    #for j in range(0,jx):
    #    for i in range(0,ix):
    #        psi[i,j]=(psi_x[i,j]+psi_y[i,j])*0.5


    # 2値化を行う
    #threshold=0.1e+10
    psi_thr=np.zeros((ix,jx))
    psi_max=np.max(psi)
    psi_min=100.0
    print('psi_max =',np.max(psi))
    while True:
        center=(psi_max+psi_min)*0.5
        for j in range(0,jx):
            for i in range(0,ix):
                if (psi[i,j] > center):
                    psi_thr[i,j]=110.0
                else:
                    psi_thr[i,j]=0.0

        ret,binary_psi = cv2.threshold(psi_thr, 100, 255, cv2.THRESH_BINARY)

        psi_label=binary_psi.astype('uint8')
        nlabels,labellmages,stats,centroids=cv2.connectedComponentsWithStats(psi_label)
        if (nlabels-1 > 1):
            psi_min=center
        else:
            psi_max=center
        if (abs(psi_max-psi_min) < 1.e+8):
            break
    print('border =',center)
    print('area = ',stats[1,4]*dx*dy)
    print('r_eff =',np.sqrt(stats[1,4]*dx*dy/math.pi))
    #print('abs = ',abs(psi_max-psi_min))

    fig=plt.figure(figsize=(10,12))
    lfac=1.e-8

    ax=fig.add_subplot(111,aspect='equal')
    ax=ax.pcolormesh(y*lfac,(x-rsun)*lfac,binary_psi,cmap='gray')
    plt.savefig(pngdir+'psi_binary_'+'{0:08d}'.format(n)+'.png')

    #磁束の計算
    t=d.read_time(0,silent=True)

    d.read_qq_2d(0,silent=True)
    #bz=abs(d.q2['bz'])
    bx=d.q2['bx']
    by=d.q2['by']
    bz=d.q2['bz']
    initial_flux=0.0

    for j in range(0,jx):
        for i in range(0,ix):
            initial_flux=initial_flux+bz[i,j]*dx*dy
    print('step 0 flux = ',initial_flux)

    t=d.read_time(n,silent=True)

    d.read_qq_2d(n,silent=True)
    #bz=abs(d.q2['bz'])
    bz=d.q2['bz']
    after_flux=0.0

    for j in range(0,jx):
        for i in range(0,ix):
            if (psi[i,j] > center):
                after_flux=after_flux+bz[i,j]*dx*dy
    print('step',n,'flux =',after_flux)

    #保持率の計算
    ratio=after_flux/initial_flux*100
    print('ratio =',ratio)

    return ratio


print("input caseid")
caseid=0
caseid=input()
caseid="d"+caseid.zfill(3)

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
ratio1=psi_compare(xmax,xmin,ymax,ymin)

print("input caseid")
caseid=0
caseid=input()
caseid="d"+caseid.zfill(3)

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
ratio2=psi_compare(xmax,xmin,ymax,ymin)

print("input caseid")
caseid=0
caseid=input()
caseid="d"+caseid.zfill(3)

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
ratio3=psi_compare(xmax,xmin,ymax,ymin)

print("input caseid")
caseid=0
caseid=input()
caseid="d"+caseid.zfill(3)

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
ratio4=psi_compare(xmax,xmin,ymax,ymin)

print("input caseid")
caseid=0
caseid=input()
caseid="d"+caseid.zfill(3)

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
ratio5=psi_compare(xmax,xmin,ymax,ymin)


#make plot
fig=plt.figure(figsize=(10,12))
#x=[256*512,512*1024,1024*2048,2048*4096]
#x=[0.10,0.15,0.25,0.40]
#y=[ratio1,ratio2,ratio3,ratio4]
x=[0.05,0.10,0.15,0.25,0.40]
y=[ratio1,ratio2,ratio3,ratio4,ratio5]


yy=np.zeros(len(y))
for m in range(0,len(y)):
    yy[m]=Decimal(str(y[m])).quantize(Decimal('0.1'),ROUND_HALF_UP)
plt.xlim(0,0.45)
plt.ylim(0,100)
#plt.xscale("log")
#plt.xlabel("grid number",fontsize=20)
plt.xlabel("λ",fontsize=20)
plt.ylabel("flux ratio (%)",fontsize=20)
#plt.text(0.05,ratio1,ratio1)
plt.plot(x,y,marker='D',linestyle='None')
for m in range(0,len(y)):
    plt.text(x[m],y[m],yy[m],fontsize=16)
