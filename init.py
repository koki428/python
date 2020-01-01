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

R2D2.init(datadir)
for key in R2D2.p:
    exec('%s = %s%s%s' % (key, 'R2D2.p["',key,'"]'))

try:
    n0
except NameError:
    n0 = 0
if  n0 > R2D2.p["nd"]:
    n0 = R2D2.p["nd"]

print('### calculation domain ###')
print('xmax - rsun = ', '{:6.2f}'.format((xmax - rsun)*1.e-8),'[Mm], xmin - rsun = ', '{:.2f}'.format((xmin - rsun)*1.e-8),'[Mm]')
print('ymax        = ', '{:6.2f}'.format(ymax*1.e-8)       ,'[Mm], ymin        = ', '{:.2f}'.format(ymin*1.e-8),'[Mm]' )
print('zmax        = ', '{:6.2f}'.format(zmax*1.e-8)       ,'[Mm], zmin        = ', '{:.2f}'.format(zmin*1.e-8),'[Mm]' )

print('')
print('### number of grid ###')
print('(ix,jx,kx)=(',ix,',',jx,',',kx,')')

print('')
print('### calculation time ###')
print('time step (nd) =',nd)
t = R2D2.read_time(nd)
print('time =','{:.2f}'.format(t/3600),' [hour]')
