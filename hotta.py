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

try:
    d
except NameError:
    ReadFlag = True
else:
    if d.p['datadir'] != datadir:
        ReadFlag = True

d = R2D2.R2D2_data(datadir)
for key in d.p:
    exec('%s = %s%s%s' % (key, 'd.p["',key,'"]'))

DprDro,tmp=np.meshgrid(dprdro,y,indexing='ij')
DprDse,tmp=np.meshgrid(dprdse,y,indexing='ij')

pr = np.zeros((ix,jx,nd+1))
vx = np.zeros((ix,jx,nd+1))
ro = np.zeros((ix,jx,nd+1))
for n in range(0,nd+1):
    print(n)
    d.read_qq_2d(n,silent=True)
    pr[:,:,n] = DprDro*d.q2['ro']+DprDse*d.q2['se']
    vx[:,:,n] = d.q2['vx']
    ro[:,:,n] = d.q2['ro']
    