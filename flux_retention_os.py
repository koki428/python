#磁束の保持率の計算
import numpy as np
import math
import sympy as sym
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import R2D2
import cv2
import sys
import os
import pandas as pd
from decimal import Decimal, ROUND_HALF_UP, ROUND_HALF_EVEN

"""
print('compare case count')
count=0
count=input()
count=int(count)
"""

# print('roll number')
# roll=input()
# roll=int(roll)

# print('case number')
# count=0
# count=input()
# count=int(count)

# print('convection number')
# con_num=0
# con_num=input()
# con_num=int(con_num)

#caseid=np.zeros(count)
# typ=np.zeros(count)
# lam=np.zeros(count)
# flux=np.zeros(count)
# lam1=np.zeros(con_num)
# lam2=np.zeros(count-con_num)
# ratio1=np.zeros(con_num)
# ratio2=np.zeros(count-con_num)

# for m in range(0,count):
    # print('input caseid id',m+1)
    # globals()['caseid'+str(m)]=0
    # globals()['caseid'+str(m)]=input()
    # globals()['caseid'+str(m)]='d'+globals()['caseid'+str(m)].zfill(3)

    # print('convection = 1, no convection = 2')
    # typ[m]=input()
    # typ[m]=int(typ[m])

    # print('λ')
    # lam[m]=input()


# for m in range(0,count):
"""
print("input caseid id",m+1)
caseid=0
caseid=input()
caseid='d'+caseid.zfill(3)

print('convection = 1, no convection = 2')
typ=0
typ=input()
typ=int(typ)
"""
#print(caseid,'λ')
# if (typ[m] == 1):
#     mm=m-np.sum(lam2 > 0)
#     #lam1[mm]=0
#     lam1[mm]=lam[m]
#     #tt=typ[m]

# if (typ[m] == 2):
#     mm=m-np.sum(lam1 > 0)
#     #lam2[mm]=0
#     lam2[mm]=lam[m]
#     #tt=typ[m]

# datadir="../run/"+globals()['caseid'+str(m)]+"/data/"
# pngdir="../figs/"+globals()['caseid'+str(m)]+"/png/"

# d=R2D2.R2D2_data(datadir)
# for key in d.p:
#     exec('%s=%s%s%s' % (key, 'd.p["',key,'"]'))

# try:
#     n0
# except NameError:
#     n0=0
# if n0>d.p["nd"]:
#     n0=d.p["nd"]

# plt.close("all")

# t0=d.read_time(0,silent=True)

# dx=(xmax-xmin)/ix
# dy=(ymax-ymin)/jx

"""    
print('input max step number')
nd=input()
nd=int(nd)
"""

##t＝0の計算
t=d.read_time(0,silent=True)
d.read_qq_2d(0,silent=True)
# bx=d.q2['bx']
# by=d.q2['by']
bz=d.q2['bz']

#平行移動
# bx=np.roll(bx,roll,axis=1)
# by=np.roll(by,roll,axis=1)
bz=np.roll(bz,roll,axis=1)

#t=0での磁束の計算
initial_flux=0.0
# for j in range(0,jx):
#     for i in range(0,ix):
#         initial_flux=initial_flux+bz[i,j]*dx*dy
bz=pd.DataFrame(bz)
initial_flux=bz.sum().sum()*dx*dy

##t=nの計算
# found=False
# for n in range(n0,nd+1):
    #print(n)
t=d.read_time(achirved_step,silent=True)
d.read_qq_2d(achirved_step,silent=True)
# bx=d.q2['bx']
# by=d.q2['by']
bz=d.q2['bz']

#平行移動
# bx=np.roll(bx,roll,axis=1)
# by=np.roll(by,roll,axis=1)
bz=np.roll(bz,roll,axis=1)

#t=achirved_stepでのpsiの計算
# psi=np.zeros((ix,jx))
# psi[0,0]=0.0
# i=0
# for j in range(1,jx-1):
#     psi[i,j]=psi[i,j-1]-(bx[i,j-1]*dy+4*bx[i,j]*dy+bx[i,j+1]*dy)/3
# for i in range(1,ix):
#     j=0
#     psi[i,j]=psi[i-1,j]+by[i,j]*dx
#     for j in range(1,jx-1):
#         psi[i,j]=psi[i,j-1]-(bx[i,j-1]*dy+4*bx[i,j]*dy+bx[i,j+1]*dy)/3

# #2値化
# psi_thr=np.zeros((ix,jx))
# psi_max=np.max(psi)
# psi_min=100.0
# while True:
#     center=(psi_max+psi_min)*0.5
#     for j in range(0,jx):
#         for i in range(0,ix):
#             if (psi[i,j] > center):
#                 psi_thr[i,j]=110.0
#             else:
#                 psi_thr[i,j]=0.0

#     ret,binary_psi = cv2.threshold(psi_thr, 100, 255, cv2.THRESH_BINARY)

#     psi_label=binary_psi.astype('uint8')
#     nlabels,labellmages,stats,centroids=cv2.connectedComponentsWithStats(psi_label)
#     if (nlabels-1 > 1):
#         psi_min=center
#     else:
#         psi_max=center
#     if (abs(psi_max-psi_min) < 1.e+1):
#         break

#中心の高さの計算
# cx=centroids[1,1]
# hight=(xmin-rsun+cx*dx)*1.e-8

# rr=np.sqrt(stats[1,4]*dx*dy/math.pi)*1.e-8

#上部の到達したか判定
# if np.any(psi_thr[950,:] > 0):
# print('step',n)
bz=pd.DataFrame(bz)
psi=pd.DataFrame(psi)
flux=bz[psi > center].sum().sum()*dx*dy

##保持率の計算
retention=flux/initial_flux*100
print(retention)
# found=True

# if found:
#     break

# if (n==nd):
#     globals()['ratio'+str(int(typ[m]))][mm]=0.0
#     print(globals()['caseid'+str(m)],'ratio =',globals()['ratio'+str(int(typ[m]))][mm])


"""
#t=nでの磁束の計算
for j in range(0,jx):
    for i in range(0,ix):
        if (psi[i,j] > center):
            flux[n]=flux[n]+bz[i,j]*dx*dy

##保持率の計算
ratio[n]=flux[n]/initial_flux*100
print('caseid',caseid,'ratio =',ratio[n])
"""
##plot
#x=[0.10,0.15,0.20,0.25,0.30,0.35,0.40]
#x=[1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7]
# x1=lam1
# x2=lam2
# y1=np.zeros(len(ratio1))
# y2=np.zeros(len(ratio2))
# yy1=np.zeros(len(y1))
# yy2=np.zeros(len(y2))

# for n in range(0,len(y1)):
#     y1[n]=globals()['ratio1'][n]
#     yy1[n]=Decimal(str(y1[n])).quantize(Decimal('0.1'),ROUND_HALF_UP)
# for n in range(0,len(y2)):
#     y2[n]=globals()['ratio2'][n]
#     yy2[n]=Decimal(str(y2[n])).quantize(Decimal('0.1'),ROUND_HALF_UP)

# fig=plt.figure(figsize=(10,6))
# ax=fig.add_subplot(1,1,1)
# ax.plot(x1,y1,marker='D',linestyle='None',label='convection')
# ax.plot(x2,y2,marker='D',linestyle='None',label='no convection')
# plt.xlim(0,0.45)
# plt.ylim(0,100)
# plt.ylabel('flux retained (%)',fontsize=22)
# plt.xlabel('$\lambda$',fontsize=22)
# #plt.xlabel('x $10^5$ [G]',fontsize=22)
# for n in range(0,len(y1)):
#     plt.text(x1[n],y1[n],yy1[n],fontsize=16)
# for n in range(0,len(y2)):
#     plt.text(x2[n],y2[n],yy2[n],fontsize=16)
# plt.legend()