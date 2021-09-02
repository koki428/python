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

n0 = 0
t = np.zeros(nd-n0+1)
ekm = np.zeros(nd-n0+1)
ekt = np.zeros(nd-n0+1)
emm = np.zeros(nd-n0+1)
emt = np.zeros(nd-n0+1)
anm = np.zeros(nd-n0+1)
bzmt = np.zeros((ix,jx,nd-n0+1))

RR, TH = np.meshgrid(x,y,indexing='ij')
ro2, tmp = np.meshgrid(ro0,y,indexing='ij')

pii8 = 1/pi/8

for n in range(n0,nd+1):
    print(n)

    d.read_time(n)
    t[n-n0] = d.t

    d.read_vc(n,silent=True)
    ekm[n-n0] = (RR**2*sin(TH)*0.5*ro2*(d.vc['vxm']**2   + d.vc['vym']**2   + d.vc['vzm']**2)).mean()/(RR**2*sin(TH)).mean()
    ekt[n-n0] = (RR**2*sin(TH)*0.5*ro2*(d.vc['vxrms']**2 + d.vc['vyrms']**2 + d.vc['vzrms']**2)).mean()/(RR**2*sin(TH)).mean()
    emm[n-n0] = (RR**2*sin(TH)*pii8*(d.vc['bxm']**2   + d.vc['bym']**2   + d.vc['bzm']**2)).mean()/(RR**2*sin(TH)).mean()
    emt[n-n0] = (RR**2*sin(TH)*pii8*(d.vc['bxrms']**2 + d.vc['byrms']**2 + d.vc['bzrms']**2)).mean()/(RR**2*sin(TH)).mean()
    anm[n-n0] = ((RR*sin(TH)*d.vc['vzm'])*RR**2*sin(TH)*ro2).sum()
    bzmt[:,:,n-n0] = d.vc['bzm']

plt.close('all')
plt.clf()
fig = plt.figure(100,figsize=(8,4))

ax = fig.add_subplot(111)
t_day = 86400
ax.plot(t/t_day,ekm,label='kin, mean',color=R2D2.blue)
ax.plot(t/t_day,ekt,label='kin, turb',color=R2D2.magenta)
ax.plot(t/t_day,emm,label='mag, mean',color=R2D2.green)
ax.plot(t/t_day,emt,label='mag, turb',color=R2D2.orange)
ax.legend()
plt.xlabel('t [day]')
plt.yscale('log')
plt.ylim(1.e4,1.e9)
plt.tight_layout()
plt.savefig(caseid+'.pdf')
