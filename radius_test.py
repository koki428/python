import cv2
import numpy as np
import R2D2
import math
import matplotlib.pyplot as plt

#def psi_compare(xmax,xmin,ymax,ymin):
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

t0=d.read_time(0,silent=True)

dx=(xmax-xmin)/jx
dy=(ymax-ymin)/ix

print('input max step number')
nd=input()
nd=int(nd)


t=d.read_time(nd,silent=True)
d.read_qq_2d(nd,silent=True)

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

#for j in range(0,jx):
#    for i in range(0,ix):
#        if (i==0) and (j==0):
#            psi_x[i,j]=0
#        elif (j >= 1) and (i==0):
#            psi_x[i,j]=psi_x[i,j-1]-bx[i,j]*dy
#        elif (i >= 1):
#            psi_x[i,j]=psi_x[i-1,j]+by[i,j]*dx
psi_y[0,0]=0.0
i=0
for j in range(1,jx):
    psi_y[i,j]=psi_y[i,j-1]-bx[i,j]*dy
for i in range(1,ix):
    j=0
    psi_y[i,j]=psi_y[i-1,j]+by[i,j]*dx
    for j in range(1,jx):
        psi_y[i,j]=psi_y[i,j-1]-bx[i,j]*dy

#for j in range(0,jx):
#    for i in range(0,ix):
#        psi[i,j]=(psi_x[i,j]+psi_y[i,j])*0.5


psi[0,0]=0.0
i=0
for j in range(1,jx):
    psi[i,j]=psi[i,j-1]-bx[i,j]*dy
for i in range(1,ix):
    j=0
    psi[i,j]=psi[i-1,j]+by[i,j]*dx
    for j in range(1,jx):
        psi[i,j]=psi[i,j-1]-bx[i,j]*dy

# 2値化を行う
#threshold=0.1e+10
psi_thr=np.zeros((ix,jx))
psi_max=np.max(psi)
psi_min=100
#print('psi_max =',np.max(psi))
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
    if (abs(psi_max-psi_min) < 1.e+8) and (nlabels-1==1):
        break

r_eff=np.zeros(nd+1)
lambda_ave=np.zeros(nd+1)
"""
#n=0(初期状態)の計算
n=n0
print(n)
t=d.read_time(n,silent=True)
d.read_qq_2d(n,silent=True)

bx=d.q2['bx']
by=d.q2['by']
bz=d.q2['bz']

area=0
for j in range(0,jx):
    for i in range(0,ix):
        if (bz[i,j] > 0):
            area=area+dx*dy
r0=np.sqrt(area/math.pi)
r_eff[0]=r0
print('r_eff =',r_eff[n])

# 2値化を行う
#threshold=0.1e+10
bz_thr=np.zeros((ix,jx))
#psi_max=np.max(psi)
#psi_min=100.0
#print('psi_max =',np.max(psi))
#while True:
#    center=(psi_max+psi_min)*0.5
for j in range(0,jx):
    for i in range(0,ix):
        if (bz[i,j] > 0):
            bz_thr[i,j]=110.0
        else:
            bz_thr[i,j]=0.0

ret,binary_psi = cv2.threshold(psi_thr, 100, 255, cv2.THRESH_BINARY)

psi_label=binary_psi.astype('uint8')
nlabels,labellmages,stats,centroids=cv2.connectedComponentsWithStats(psi_label)

#<λ>の計算
rr=0
count=0
for j in range(0,jx):
    for i in range(0,ix):
        if (bz[i,j] > 0):
            r=np.sqrt((i*dx-centroids[1,0]*dx)**2+(j*dy-centroids[1,1]*dy)**2)
            bs=np.sqrt((bx[i,j])**2+(by[i,j])**2)
            rr=rr+bs/(r*abs(bz[i,j]))
            count=count+1
lambda_ave[n]=r_eff[n]*rr/count
print('lambda average =',lambda_ave[n])

fig=plt.figure(figsize=(8,10))
lfac=1.e-8
ax=fig.add_subplot(111,aspect='equal')
ax=ax.pcolormesh(y*lfac,(x-rsun)*lfac,binary_psi,cmap='gray')
if (n==n0):
    fig.tight_layout(pad=0.1)
plt.pause(0.1)
plt.savefig(pngdir+'psi_binary_'+'{0:08d}'.format(n)+'.png')
if(n != nd-1):
    clf()
"""
fig=plt.figure(figsize=(8,10))
lfac=1.e-8
#for n in range(n0,nd+1):
    #print(n)
ax3=fig.add_subplot(1,2,1,aspect='equal')
ax3=ax3.pcolormesh(y*lfac,(x-rsun)*lfac,binary_psi,cmap='gray')
ax4=fig.add_subplot(1,2,2,aspect='equal')
ax4=ax4.pcolormesh(y*lfac,(x-rsun)*lfac,psi,cmap='gist_stern')
ax44=plt.contour(y*lfac,(x-rsun)*lfac,psi,colors='g',levels=[center])
fig.colorbar(ax4)
if (nlabels-1 > 1):
    area_max=0
    for m in range(1,nlabels):
        if (stats[m,4] > area_max):
            area_max=stats[m,4]
            mm=m
    stats[1,4]=area_max
r_eff[nd]=np.sqrt(stats[1,4]*dx*dy/math.pi)
print('r_eff = ',r_eff[nd])
print('rsun*0.005*2 =',rsun*0.005*2)
print('max =',np.max(psi))
print('閾値 =',center)

"""
t=d.read_time(0,silent=True)
d.read_qq_2d(0,silent=True)

bx=d.q2['bx']
by=d.q2['by']
bz=d.q2['bz']

#bb=np.sqrt((bx)**2+(by)**2+(bz)**2)

#psi cal
psi_x=np.zeros((ix,jx))
psi_y=np.zeros((ix,jx))
psi=np.zeros((ix,jx))

psi_x[0,0]=0.0
j=0
for i in range(1,ix):
    psi_x[i,j]=psi_x[i-1,j]+by[i,j]*dx
    #psi_x[i,j]=psi_x[i-1,j]+by[j,i]*dx
for j in range(1,jx):
    i=0
    psi_x[i,j]=psi_x[i,j-1]-bx[i,j]*dy
    #psi_x[i,j]=psi_x[i,j-1]-bx[j,i]*dy
    for i in range(1,ix):
        psi_x[i,j]=psi_x[i-1,j]+by[i,j]*dx
        #psi_x[i,j]=psi_x[i-1,j]+by[j,i]*dx

psi_y[0,0]=0.0
i=0
for j in range(1,jx):
    psi_y[i,j]=psi_y[i,j-1]-bx[i,j]*dy
    #psi_y[i,j]=psi_y[i,j-1]-bx[j,i]*dy
for i in range(1,ix):
    j=0
    psi_y[i,j]=psi_y[i-1,j]+by[i,j]*dx
    #psi_y[i,j]=psi_y[i-1,j]+by[j,i]*dx
    for j in range(1,jx):
        psi_y[i,j]=psi_y[i,j-1]-bx[i,j]*dy
        #psi_y[i,j]=psi_y[i,j-1]-bx[j,i]*dy

for j in range(0,jx):
    for i in range(0,ix):
        psi[i,j]=(psi_x[i,j]+psi_y[i,j])*0.5


#2値化を行う
#threshold=0.1e+10
psi_thr=np.zeros((ix,jx))
#psi_max=np.max(psi)
#psi_min=100.0
#print('psi_max =',np.max(psi))
#while True:
#    center=(psi_max+psi_min)*0.5
for j in range(0,jx):
    for i in range(0,ix):
        if (psi[i,j] > center):
            psi_thr[i,j]=110.0
        else:
            psi_thr[i,j]=0.0

ret,binary_psi = cv2.threshold(psi_thr, 100, 255, cv2.THRESH_BINARY)

psi_label=binary_psi.astype('uint8')
nlabels,labellmages,stats,centroids=cv2.connectedComponentsWithStats(psi_label)
    #    if (nlabels-1 > 1):
    #        psi_min=center
    #    else:
    #        psi_max=center
    #    if (abs(psi_max-psi_min) < 0.5) and (nlabels-1==1):
    #        break
    #print('center =',center)
ax=fig.add_subplot(2,2,3,aspect='equal')
ax=ax.pcolormesh(y*lfac,(x-rsun)*lfac,binary_psi,cmap='gray')

ax2=fig.add_subplot(2,2,4,aspect='equal')
ax2=ax2.pcolormesh(y*lfac,(x-rsun)*lfac,psi,cmap='gist_stern')
ax22=plt.contour(y*lfac,(x-rsun)*lfac,psi,colors='g',levels=[center])
fig.colorbar(ax2)   
   #if (n==n0):
    #    fig.tight_layout(pad=0.1)
    #plt.pause(0.1)
    #plt.savefig(pngdir+'psi_binary_'+'{0:08d}'.format(n)+'.png')
    #if(n != nd-1):
    #    clf()

    #print('stats = ',stats[1,:])
    #print('中心 =',centroids[1,:])
    #print('中心 =',centroids[1,1])

    #print('area =',stats[1,4]*dx*dy)

    #有効半径の計算
    #area=0
    #if (n==0):
    #    for j in range(0,jx):
    #        for i in range(0,ix):
    #            if (bz[i,j] > 0):
    #                area=area+dx*dy
    #print('area =',area)
    #r_eff2=np.sqrt(area/math.pi)
    #print('r_eff2 =',r_eff2)

    #r0=rsun*0.005*2
if (nlabels-1 > 1):
    area_max=0
    for m in range(1,nlabels):
        if (stats[m,4] > area_max):
            area_max=stats[m,4]
            mm=m
    stats[1,4]=area_max
r_eff[nd]=np.sqrt(stats[1,4]*dx*dy/math.pi)
print('r_eff = ',r_eff[nd])
print('rsun*0.005*2 =',rsun*0.005*2)
print('max =',np.max(psi))
print('閾値 =',center)
    #<λ>の計算
    #rr=0
    #count=0
    #if (nlabels-1 > 1):
    #    centroids[1,:]=centroids[mm,:]
    #for j in range(0,jx):
    #    for i in range(0,ix):
    #        if (psi[i,j] > center):
    #            r=np.sqrt((i*dx-centroids[1,0]*dx)**2+(j*dy-centroids[1,1]*dy)**2)
    #            #r=np.sqrt((j*dy-centroids[1,0]*dy)**2+(i*dx-centroids[1,1]*dx)**2)
    #            bs=np.sqrt((bx[i,j])**2+(by[i,j])**2)
    #            #bs=np.sqrt((bx[j,i])**2+(by[j,i])**2)
    #            rr=rr+bs/(r*abs(bz[i,j]))
    #            #rr=rr+bs/(r*abs(bz[j,i]))
    #            count=count+1
    #lambda_ave[n]=r_eff[n]*rr/count
    #print('r_eff =',r_eff[n]/r_eff[0])
    #print('lambda average =',lambda_ave[n]/lambda_ave[0])

#for m in range(0,10):
#    r_eff[m]=r0
#    lambda_ave[m]=0.15

#plt.clf()
#plt.plot(r_eff/(rsun*0.005*2),lambda_ave/0.25,marker='D',linestyle='None')
#plt.xlim(0.5,2.0)
    #print('lambda average =',lambda_ave)

    #面積の計算
    #area=0
    #for j in range(0,jx):
    #    for i in range(0,ix):
    #        if (psi[i,j] > center):
    #            area=area+dx*dy
    #print('area = ',area)
    #r_eff=np.sqrt(area/math.pi)
    #print('r_eff = ',r_eff)

    #print('psi = ',center)
    #print('abs = ',abs(psi_max-psi_min))

    #fig=plt.figure(figsize=(10,12))
    #lfac=1.e-8

    #ax=fig.add_subplot(111,aspect='equal')
    #ax=ax.pcolormesh(y*lfac,(x-rsun)*lfac,binary_psi,cmap='gray')
"""
"""
    #磁束の計算
    t=d.read_time(0,silent=True)

   d.read_qq_2d(0,silent=True)
    bz=abs(d.q2['bz'])
    initial_flux=0.0

    for j in range(0,jx):
        for i in range(0,ix):
            initial_flux=initial_flux+bz[i,j]*dx*dy
    print('step 0 flux = ',initial_flux)

    t=d.read_time(n,silent=True)

    d.read_qq_2d(n,silent=True)
    bz=abs(d.q2['bz'])
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
#x=[256*512,512*1024,1024*2048,2048*4096]
#y=[ratio1,ratio2,ratio3,ratio4]
x=[0.05,0.10,0.15,0.25,0.40]
y=[ratio1,ratio2,ratio3,ratio4,ratio5]

plt.ylim(0,100)
#plt.xscale("log")
#plt.xlabel("grid number")
plt.xlabel("λ")
plt.ylabel("flux ratio (%)")
plt.plot(x,y,marker='D',linestyle='None')
"""
