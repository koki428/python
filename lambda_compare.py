import numpy as np
import math
import matplotlib.pyplot as plt
import R2D2
import sys
import os

#b_border_lam_000=2000
b_border_lam_005=5000
b_border_lam_010=2500
b_border_lam_015=2500
b_border_lam_025=2500
b_border_lam_040=2500

#1grid scale cal
def x_grid_scale_cal(xmax,xmin,ix):
    dx=(xmax-xmin)/ix
    return dx

def y_grid_scale_cal(ymax,ymin,jx):
    dy=(ymax-ymin)/jx
    return dy

#main tube flux cal
def initial_maintube_cal(bb,jx,ix,dx,dy):
    btube=1.0e+5
    flux_initial=0
    btotal=0

    for l in range(0,jx):
        for k in range(0,ix):
            if(bb[k,l]>btube*np.exp(-4)):
                btotal=btotal+bb[k,l]*dx*dy

    #print('step 0 b_total = ',btotal)
    flux_initial=btotal
    #print('step 0 flux = ',flux_initial)

    return flux_initial

def after_maintube_cal(bb,b_border,jx,ix,dx,dy):
    t=0
    frag=0
    flux=0
    tube=[0]*1000
    flux_max=0
    under=0

    for l in range(0,jx):
        for k in range(0,ix):
            if (bb[k,l] > b_border):
                for m in range(0,ix):
                    if (bb[m,l-1] < b_border) and (bb[k-1,l] < b_border) and (frag != 1):
                        t=t+1
                        frag=1
                tube[t]=tube[t]+bb[k,l]*dx*dy
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
    #print('flux max =',flux_max)

    return flux_max

#ratio cal
def ratio_cal(flux_initial,flux_max):
    ratio=flux_max/flux_initial*100
    return ratio

#####################################################################################################
#lambda=0.05
print("input lambda=0.05 caseid id (3digit)")
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

#main tube flux cal
flux_initial_lam_005=initial_maintube_cal(bb,jx,ix,dx,dy)
print('step 0 lam=0.05 = ',flux_initial_lam_005)

#60 step cal
t=d.read_time(60,silent=True)

#read value
d.read_qq_2d(60,silent=True)

bb=abs(d.q2['bz'])

#main tube flux cal
flux_max_lam_005=after_maintube_cal(bb,b_border_lam_005,jx,ix,dx,dy)
print('step 60 lam=0.05 = ',flux_max_lam_005)

#ratio
ratio_lam_005=ratio_cal(flux_initial_lam_005,flux_max_lam_005)
print('ratio lam=0.05 = ',ratio_lam_005)
###########################################################################

#lambda=0.10
print("input lambda=0.10 caseid id (3digit)")
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

#main tube flux cal
flux_initial_lam_010=initial_maintube_cal(bb,jx,ix,dx,dy)
print('step 0 lam=0.10 = ',flux_initial_lam_010)

#60 step cal
t=d.read_time(60,silent=True)

#read value
d.read_qq_2d(60,silent=True)

bb=abs(d.q2['bz'])

#main tube flux cal
flux_max_lam_010=after_maintube_cal(bb,b_border_lam_010,jx,ix,dx,dy)
print('step 60 lam=0.10  = ',flux_max_lam_010)

#ratio
ratio_lam_010=ratio_cal(flux_initial_lam_010,flux_max_lam_010)
print('ratio lam=0.10 = ',ratio_lam_010)
##################################################################################

#lambda=0.15 cal
print("input lambda=0.15 caseid id (3digit)")
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

#main tube flux cal
flux_initial_lam_015=initial_maintube_cal(bb,jx,ix,dx,dy)
print('step 0 lam=0.15 = ',flux_initial_lam_015)

#60 step cal
t=d.read_time(60,silent=True)

#read value
d.read_qq_2d(60,silent=True)

bb=abs(d.q2['bz'])

#main tube flux cal
flux_max_lam_015=after_maintube_cal(bb,b_border_lam_015,jx,ix,dx,dy)
print('step 60 lam=0.15 = ',flux_max_lam_015)

#ratio
ratio_lam_015=ratio_cal(flux_initial_lam_015,flux_max_lam_015)
print('ratio lam=0.15 = ',ratio_lam_015)
###############################################################################

#lambda=0.25cal
print("input lambda=0.25 caseid id (3digit)")
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

#main tube flux cal
flux_initial_lam_025=initial_maintube_cal(bb,jx,ix,dx,dy)
print('step 0 lam=0.25 = ',flux_initial_lam_025)

#60 step cal
t=d.read_time(60,silent=True)

#read value
d.read_qq_2d(60,silent=True)

bb=abs(d.q2['bz'])

#main tube flux cal
flux_max_lam_025=after_maintube_cal(bb,b_border_lam_025,jx,ix,dx,dy)
print('step 60 lam=0.25 = ',flux_max_lam_025)

#ratio
ratio_lam_025=ratio_cal(flux_initial_lam_025,flux_max_lam_025)
print('ratio lam=0.25 = ',ratio_lam_025)
####################################################################################

#lambda=0.40cal
print("input lambda=0.40 caseid id (3digit)")
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

#main tube flux cal
flux_initial_lam_040=initial_maintube_cal(bb,jx,ix,dx,dy)
print('step 0 lam=0.40 = ',flux_initial_lam_040)

#60 step cal
t=d.read_time(60,silent=True)

#read value
d.read_qq_2d(60,silent=True)

bb=abs(d.q2['bz'])

#main tube flux cal
flux_max_lam_040=after_maintube_cal(bb,b_border_lam_040,jx,ix,dx,dy)
print('step 60 lam=0.40 = ',flux_max_lam_040)

#ratio
ratio_lam_040=ratio_cal(flux_initial_lam_040,flux_max_lam_040)
print('ratio lam=0.40 = ',ratio_lam_040)
####################################################################################

#make plot
x=[0.05,0.10,0.15,0.25,0.40]
y=[ratio_lam_005,ratio_lam_010,ratio_lam_015,ratio_lam_025,ratio_lam_040]

plt.xlim(0,0.45)
plt.ylim(0,100)

plt.xlabel("Twist parameter")
plt.ylabel("flux ratio (%)")

plt.plot(x,y,marker="D",linestyle='None')
