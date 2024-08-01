import numpy as np
import math
import matplotlib.pyplot as plt
import R2D2
import sys
import os

b_border_256x512=2000
b_border_512x1024=2500
b_border_1024x2048=2500
#b_border_2024x4096=

#1grid scale cal
def x_grid_scale_cal(xmax,xmin,ix):
    dx=(xmax-xmin)/ix
    return dx

def y_grid_scale_cal(ymax,ymin,jx):
    dy=(ymax-ymin)/jx
    return dy

#main tube flux cal
def initial_maintube_cal(bb,ro,jx,ix,dx,dy):
    btube=1.0e+5
    flux_initial=0
    btotal=0
    mass=0

    for l in range(0,jx):
        for k in range(0,ix):
            if(bb[k,l]>btube*np.exp(-4)):
                btotal=btotal+bb[k,l]*dx*dy
                mass=mass+ro[k,l]*dx*dy

    #print('step 0 b_total = ',btotal)
    ratio_initial=mass/btotal
    #print('step 0 flux = ',flux_initial)

    return ratio_initial

def after_maintube_cal(bb,ro,b_border,jx,ix,dx,dy):
    t=0
    frag=0
    flux=0
    tube=[0]*100
    flux_max=0
    under=0
    mass=[0]*100

    for l in range(0,jx):
        for k in range(0,ix):
            if (bb[k,l] > b_border):
                for m in range(0,ix):
                    if (bb[m,l-1] < b_border) and (bb[k-1,l] < b_border) and (frag != 1):
                        t=t+1
                        frag=1
                tube[t]=tube[t]+bb[k,l]*dx*dy
                mass[t]=mass[t]+ro[k,l]*dx*dy

            else:
                under=under+1
                if (under>ix-1):
                    frag=0
        under=0

    tmax=t
    #print('tmax = ',tmax)

    for t in range(1,tmax+1):
        if (tube[t]>flux_max):
            flux_max=tube[t]
            tt=t
    #print('flux max =',flux_max)
    ratio_after=mass[tt]/flux_max
    return ratio_after

#ratio cal
def ratio_cal(ratio_initial,ratio_after):
    ratio=ratio_after/ratio_initial*100
    return ratio


#256x512cal
print("input 256x512 caseid id (3digit)")
caseid=0
caseid=input()
caseid="d"+caseid.zfill(3)

datadir="../run/"+caseid+"/data/"

d=R2D2.R2D2_data(datadir)
for key in d.p:
    exec('%s = %s%s%s' % (key, 'd.p["',key,'"]'))

try:
    n0
except NameError:
    n0=0

if n0>d.p["nd"]:
    n0=d.p["nd"]

#print("Maximum time step = ",nd," time ="\
#        ,dtout*float(nd)/3600./24.," [day]")

plt.close('all')

#read time
t0=d.read_time(0,silent=True)

dx=x_grid_scale_cal(xmax,xmin,ix)
dy=y_grid_scale_cal(ymax,ymin,jx)

#initial step cal
t=d.read_time(0,silent=True)

#read value
d.read_qq_2d(0,silent=True)

bb=abs(d.q2['bz'])
ro=abs(d.q2['ro'])

#main tube flux cal
ratio_initial_256x512=initial_maintube_cal(bb,ro,jx,ix,dx,dy)
print('step 0 mass/flux 256x512 = ',ratio_initial_256x512)

#60 step cal
t=d.read_time(60,silent=True)

#read value
d.read_qq_2d(60,silent=True)

bb=abs(d.q2['bz'])
ro=abs(d.q2['ro'])

#main tube flux cal
ratio_after_256x512=after_maintube_cal(bb,ro,b_border_256x512,jx,ix,dx,dy)
print('step 60 mass/flux 256x512 = ',ratio_after_256x512)

#ratio
ratio_256x512=ratio_cal(ratio_initial_256x512,ratio_after_256x512)
print('ratio 256x512 = ',ratio_256x512)

#512x1024cal
print("input 512x1024 caseid id (3digit)")
caseid=0
caseid=input()
caseid="d"+caseid.zfill(3)

datadir="../run/"+caseid+"/data/"

d=R2D2.R2D2_data(datadir)
for key in d.p:
    exec('%s = %s%s%s' % (key, 'd.p["',key,'"]'))

try:
    n0
except NameError:
    n0=0

if n0>d.p["nd"]:
    n0=d.p["nd"]

#print("Maximum time step = ",nd," time ="\
#        ,dtout*float(nd)/3600./24.," [day]")

plt.close('all')

#read time
t0=d.read_time(0,silent=True)

dx=x_grid_scale_cal(xmax,xmin,ix)
dy=y_grid_scale_cal(ymax,ymin,jx)

#initial step cal
t=d.read_time(0,silent=True)

#read value
d.read_qq_2d(0,silent=True)

bb=abs(d.q2['bz'])
ro=abs(d.q2['ro'])

#main tube flux cal
ratio_initial_512x1024=initial_maintube_cal(bb,ro,jx,ix,dx,dy)
print('step 0 mass/flux 512x1024 = ',ratio_initial_512x1024)

#60 step cal
t=d.read_time(60,silent=True)

#read value
d.read_qq_2d(60,silent=True)

bb=abs(d.q2['bz'])
ro=abs(d.q2['ro'])

#main tube flux cal
ratio_after_512x1024=after_maintube_cal(bb,ro,b_border_512x1024,jx,ix,dx,dy)
print('step 60 mass/flux 512x1024 = ',ratio_after_512x1024)

#ratio
ratio_512x1024=ratio_cal(ratio_initial_512x1024,ratio_after_512x1024)
print('ratio 512x1024 = ',ratio_512x1024)

#1024x2048cal
print("input 1024x2048 caseid id (3digit)")
caseid=0
caseid=input()
caseid="d"+caseid.zfill(3)

datadir="../run/"+caseid+"/data/"

d=R2D2.R2D2_data(datadir)
for key in d.p:
    exec('%s = %s%s%s' % (key, 'd.p["',key,'"]'))

try:
    n0
except NameError:
    n0=0

if n0>d.p["nd"]:
    n0=d.p["nd"]

#print("Maximum time step = ",nd," time ="\
#        ,dtout*float(nd)/3600./24.," [day]")

plt.close('all')

#read time
t0=d.read_time(0,silent=True)

dx=x_grid_scale_cal(xmax,xmin,ix)
dy=y_grid_scale_cal(ymax,ymin,jx)

#initial step cal
t=d.read_time(0,silent=True)

#read value
d.read_qq_2d(0,silent=True)

bb=abs(d.q2['bz'])
ro=abs(d.q2['ro'])

#main tube flux cal
ratio_initial_1024x2048=initial_maintube_cal(bb,ro,jx,ix,dx,dy)
print('step 0 mass/flux 1024x2048 = ',ratio_initial_1024x2048)

#60 step cal
t=d.read_time(60,silent=True)

#read value
d.read_qq_2d(60,silent=True)

bb=abs(d.q2['bz'])
ro=abs(d.q2['ro'])

#main tube flux cal
ratio_after_1024x2048=after_maintube_cal(bb,ro,b_border_1024x2048,jx,ix,dx,dy)
print('step 60 mass/flux 1024x2048 = ',ratio_after_1024x2048)

#ratio
ratio_1024x2048=ratio_cal(ratio_initial_1024x2048,ratio_after_1024x2048)
print('ratio 1024x2048 = ',ratio_1024x2048)

#2048x4096cal
'''
print("input 2048x4096 caseid id (3digit)")
caseid=0
caseid=input()
caseid="d"+caseid.zfill(3)

datadir="../run/"+caseid+"/data/"

d=R2D2.R2D2_data(datadir)
for key in d.p:
    exec('%s = %s%s%s' % (key, 'd.p["',key,'"]'))

try:
    n0
except NameError:
    n0=0

if n0>d.p["nd"]:
    n0=d.p["nd"]

print("Maximum time step = ",nd," time ="\
        ,dtout*float(nd)/3600./24.," [day]")

plt.close('all')

#read time
t0=d.read_time(0,silent=True)

dx=x_grid_scale_cal(xmax,xmin,ix)
dy=y_grid_scale_cal(ymax,ymin,jx)

#initial step cal
t=d.read_time(0,silent=True)

#read value
d.read_qq_2d(0,silent=True)

bb=abs(d.q2['bz'])
ro=abs(d.q2['ro'])

#main tube flux cal
ratio_initial_2048x4096=initial_maintube_cal(bb,ro,jx,ix,dx,dy)
print('step 0 mass/flux 2048x4096 = ',ratio_initial_2048x4096)

#60 step cal
t=d.read_time(60,silent=True)

#read value
d.read_qq_2d(60,silent=True)

bb=abs(d.q2['bz'])
ro=abs(d.q2['ro'])

#main tube flux cal
ratio_after_2048x4096=after_maintube_cal(bb,ro,b_border_2048x4096,jx,ix,dx,dy)
print('step 60 mass/flux 2048x4096 = ',ratio_after_2048x4096)

#ratio
ratio_2048x4096=ratio_cal(ratio_initial_2048x4096,ratio_after_2048x4096)
print('ratio 2048x4096 = ',ratio_2048x4096)
'''

#make plot
x=[256*512,512*1024,1024*2048]
y=[ratio_256x512,ratio_512x1024,ratio_1024x2048]

plt.ylim(0,100)

plt.xscale("log")

plt.xlabel("grid number")
plt.ylabel("chenge M/Φ (%)")

plt.plot(x,y,marker="D")
