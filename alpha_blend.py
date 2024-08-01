#磁場と熱対流のx方向の合成
import numpy as np
import matplotlib.pyplot as plt
import cv2
import R2D2

print("input caseid id (3 digit)")
caseid = 0
caseid = input()
caseid = "d"+caseid.zfill(3)

datadir="../run/"+caseid+"/data/"
pngdir="../figs/"+caseid+"/png/"

d = R2D2.R2D2_data(datadir)
for key in d.p:
    exec('%s = %s%s%s' % (key, 'd.p["',key,'"]'))

try:
    n0
except NameError:
    n0 = 0
if  n0 > d.p["nd"]:
    n0 = d.p["nd"]

print("Maximum time step= ",nd," time ="\
          ,dtout*float(nd)/3600./24.," [day]")

print('input max step number')
nd=input()
nd=int(nd)

for n in range(n0,nd):
    print('step',n)
    bz_img=cv2.imread(pngdir+'py_bz'+'{0:08d}'.format(n)+'.png')
    vx_img=cv2.imread('../figs/d004/png/convection_vector_nobar'+'{0:08d}'.format(n)+'.png')

    alpha=0.5
    blended=cv2.addWeighted(bz_img,alpha,vx_img,1-alpha,0)

    plt.imshow(cv2.cvtColor(blended,cv2.COLOR_BGR2RGB))
    plt.axis('off')

    if(n == n0):
        plt.tight_layout(pad=0.5)
    
    plt.pause(0.3)
    plt.savefig(pngdir+"py_bz_convection_blend"+'{0:08d}'.format(n)+".png")

    if(n != nd):
        clf()
