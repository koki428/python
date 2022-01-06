import numpy as np
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

ReadFlag = False
try:
    d
except NameError:
    ReadFlag = True
else:
    if d.p['datadir'] != datadir:
        ReadFlag = True

if ReadFlag:
    d = R2D2.R2D2_data(datadir)
    for key in d.p:
        exec('%s = %s%s%s' % (key, 'd.p["',key,'"]'))

try:
    n0
except NameError:
    n0 = 0
if  n0 > d.p["nd"]:
    n0 = d.p["nd"]

if geometry == 'Cartesian':
    print('### calculation domain ###')
    print('xmax - rsun = ', '{:6.2f}'.format((xmax - rsun)*1.e-8),'[Mm], xmin - rsun = ', '{:6.2f}'.format((xmin - rsun)*1.e-8),'[Mm]')
    print('ymax        = ', '{:6.2f}'.format(ymax*1.e-8)       ,'[Mm], ymin        = ', '{:6.2f}'.format(ymin*1.e-8),'[Mm]' )
    print('zmax        = ', '{:6.2f}'.format(zmax*1.e-8)       ,'[Mm], zmin        = ', '{:6.2f}'.format(zmin*1.e-8),'[Mm]' )

    print('')
    print('### number of grid ###')
    print('(ix,jx,kx)=(',ix,',',jx,',',kx,')')

    print('')
    print('### calculation time ###')
    print('time step (nd) =',nd)
    t = d.read_time(nd)
    print('time =','{:.2f}'.format(t/3600),' [hour]')

if geometry == 'Spherical':
    pi2rad = 180/pi
    print('### calculation domain ###')
    print('xmax - rsun = ', '{:6.2f}'.format((xmax - rsun)*1.e-8),'[Mm],  xmin - rsun = ', '{:6.2f}'.format((xmin - rsun)*1.e-8),'[Mm]')
    print('ymax        = ', '{:6.2f}'.format(ymax*pi2rad)        ,'[rad], ymin        = ', '{:6.2f}'.format(ymin*pi2rad),'[rad]' )
    print('zmax        = ', '{:6.2f}'.format(zmax*pi2rad)        ,'[rad], zmin        = ', '{:6.2f}'.format(zmin*pi2rad),'[rad]' )

    print('')
    print('### number of grid ###')
    print('(ix,jx,kx)=(',ix,',',jx,',',kx,')')

    print('')
    print('### calculation time ###')
    print('time step (nd) =',nd)
    t = d.read_time(nd)
    print('time =','{:.2f}'.format(t/3600),' [hour]')

if geometry == 'YinYang':
    pi2rad = 180/pi
    print('### calculation domain ###')
    print('xmax - rsun = ', '{:6.2f}'.format((xmax - rsun)*1.e-8),'[Mm],  xmin - rsun = ', '{:6.2f}'.format((xmin - rsun)*1.e-8),'[Mm]')
    print('Yin-Yang grid is used to cover the whole sphere')

    print('')
    print('### number of grid ###')
    print('(ix,jx,kx)=(',ix,',',jx,',',kx,')')

    print('')
    print('### calculation time ###')
    print('time step (nd) =',nd)
    t = d.read_time(nd)
    print('time =','{:.2f}'.format(t/3600),' [hour]')
    
