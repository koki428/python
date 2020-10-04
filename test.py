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

plt.clf()
plt.close('all')

#d.read_qq_slice(5,'x',0)
d.read_qq_slice(5,'x',nd_tau)

fig = plt.figure(1)
ax = fig.add_subplot(111,projection='mollweide')

#vmin = -2.e-1
#vmax =  2.e-1

vmin = -2.e3
vmax =  2.e3

var = 'vx'

ax.pcolormesh(zog_yy[jxg_yy//2:jxg_yy,:],yog_yy[jxg_yy//2:jxg_yy,:]-0.5*pi,d.ql_yan[var][jxg_yy//2:jxg_yy,:]
              ,vmin=vmin,vmax=vmax,shading='auto')

ax.pcolormesh(zog_yy[0:jxg_yy//2,:],yog_yy[0:jxg_yy//2,:]-0.5*pi,d.ql_yan[var][0:jxg_yy//2,:]              
              ,vmin=vmin,vmax=vmax,shading='auto')
ax.pcolormesh(zzg_yy,yyg_yy-0.5*pi,d.ql_yin[var],vmin=vmin,vmax=vmax,shading='auto')
ax.set_xticklabels('')
ax.set_yticklabels('')

