import cv2
import numpy as np
import R2D2
import math
import matplotlib.pyplot as plt

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

print('input start step number')
n0=input()
n0=int(n0)

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
psi_min=100.0
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

pp=1.2
center=center*pp

r_eff=np.zeros(nd+1)
lambda_ave=np.zeros(nd+1)
fig=plt.figure(figsize=(7,9))
lfac=1.e-8
pre=np.zeros(nd+1)
xp=np.zeros(nd+1)
R=np.zeros(nd+1)

for n in range(n0,nd+1):
    print(n)
    t=d.read_time(n,silent=True)
    d.read_qq_2d(n,silent=True)

    bx=d.q2['bx']
    by=d.q2['by']
    bz=d.q2['bz']


    #psi cal
    psi_x=np.zeros((ix,jx))
    psi_y=np.zeros((ix,jx))
    psi=np.zeros((ix,jx))

    psi[0,0]=0.0
    i=0
    for j in range(1,jx):
        psi[i,j]=psi[i,j-1]-bx[i,j]*dy
    for i in range(1,ix):
        j=0
        psi[i,j]=psi[i-1,j]+by[i,j]*dx
        for j in range(1,jx):
            psi[i,j]=psi[i,j-1]-bx[i,j]*dy
    """
    psi_thr=np.zeros((ix,jx))
    psi_max=np.max(psi)
    psi_min=100.0
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
    """
   #2値化を行う
    psi_thr=np.zeros((ix,jx))
    for j in range(0,jx):
        for i in range(0,ix):
            if (psi[i,j] > center):
                psi_thr[i,j]=110.0
            else:
                psi_thr[i,j]=0.0

    ret,binary_psi = cv2.threshold(psi_thr, 100, 255, cv2.THRESH_BINARY)

    psi_label=binary_psi.astype('uint8')
    nlabels,labellmages,stats,centroids=cv2.connectedComponentsWithStats(psi_label)
    
    ax=fig.add_subplot(1,2,1,aspect='equal')
    ax=ax.pcolormesh(y*lfac,(x-rsun)*lfac,binary_psi,cmap='gray')
    ax2=fig.add_subplot(1,2,2,aspect='equal')
    ax2=ax2.pcolormesh(y*lfac,(x-rsun)*lfac,psi,vmin=0,cmap='gist_stern')
    ax22=plt.contour(y*lfac,(x-rsun)*lfac,psi,colors='g',levels=[center])
    if (n==n0):
        fig.tight_layout(pad=0.1)
    plt.pause(0.1)
    plt.savefig(pngdir+'psi_binary_contour_'+'{0:08d}'.format(n)+'.png')
    if(n != nd-1):
        clf()

    if (nlabels-1 > 1):
        area_max=0
        for m in range(1,nlabels):
            if (stats[m,4] > area_max):
                area_max=stats[m,4]
                mm=m
        stats[1,4]=area_max
    r_eff[n]=np.sqrt(stats[1,4]*dx*dy/math.pi)
    print('r_eff = ',r_eff[n])

    idx=np.unravel_index(np.argmax(bz),bz.shape)
    cxx=int(idx[0])
    cyy=int(idx[1])

    pr=d.q2['pr']
    Pr0,tmp=np.meshgrid(pr0,y,indexing='ij')
    pre[n]=Pr0[cxx,cyy]+pr[cxx,cyy]

    xp[n]=pre[n]/pre[0]
    
    gam=5/3
    R[n]=r_eff[0]*xp[n]**(-1/(2*gam))

    #<λ>の計算
    rr=0
    count=0
    if (nlabels-1 > 1):
        centroids[1,:]=centroids[mm,:]
    #print('(cx,cy) = (',centroids[1,0],',',centroids[1,1],')')
    cx=round(centroids[1,0],0)
    cy=round(centroids[1,1],0)

    #idx=np.unravel_index(np.argmax(bz),bz.shape)
    #cy=int(idx[0])
    #cx=int(idx[1])
    for j in range(0,jx):
        for i in range(0,ix):
            if (psi[i,j] > center):
                #print('(i,j) = (',i,',',j,')')
                #r=np.sqrt((j*dx-centroids[1,0]*dx)**2+(i*dy-centroids[1,1]*dy)**2)
                dis_x=(j-cx)*dx
                dis_y=(i-cy)*dy
                r=np.sqrt(dis_x**2+dis_y**2)
                #bs=np.sqrt((bx[i,j])**2+(by[i,j])**2)
                if (r != 0):
                    bs=-dis_y/r*by[i,j]+dis_x/r*bx[i,j]
                    rr=rr+bs/(r*bz[i,j])
                    count=count+1
    lambda_ave[n]=r_eff[n]*rr/count
    print('lambda average =',lambda_ave[n])
    print('r_eff/R0 =',r_eff[n]/(rsun*0.005))
    #print('lambda average =',lambda_ave[n])
    print('lam/lam0 =',lambda_ave[n]/0.25)

print(np.max(abs(r_eff/(rsun*0.005)-lambda_ave/0.25)))

#for n  in range(0,13):
#    r_eff[n]=r_eff[13]
#    lambda_ave[n]=lambda_ave[13]

plt.clf()
plt.axes().set_aspect('equal')

plt.plot(r_eff/(rsun*0.005),lambda_ave/0.25,marker='D',linestyle='None',label=r"$\left<\lambda\right>=r_{eff}\left<\frac{B_{\theta}}{r^{'}B_{z}}\right>$")

x=np.linspace(0,10,10)
y=x
"""
rr_eff=np.zeros(nd+1-13)
llambda_ave=np.zeros(nd+1-13)
for n in range(0,nd+1-13):
    rr_eff[n]=r_eff[n+13]
    llambda_ave[n]=lambda_ave[n+13]
ab=np.polyfit(rr_eff/(rsun*0.005),llambda_ave/0.25,1)
print(ab)
yy=np.poly1d(ab)(rr_eff/(rsun*0.005))
"""
plt.plot(x,y,label=r"$\frac{\lambda}{\lambda_{0}}=\frac{R}{R_{0}}$")
#plt.plot(rr_eff/(rsun*0.005),yy)

#plt.plot(r_eff/(rsun*0.005),R/(rsun*0.005))
#plt.xlim([np.min(r_eff)/r_eff[0],np.max(r_eff)/r_eff[0]])
plt.xlim(1.0,3.0)
#plt.ylim([np.min(lambda_ave)/lambda_ave[0],np.max(lambda_ave)/lambda_ave[0]])
plt.ylim(1.0,3.0)
plt.xlabel('radius (R$_{0}$)',fontsize=20)
plt.ylabel('<λ> (λ$_{0}$)',fontsize=20)
plt.legend(fontsize=20)
