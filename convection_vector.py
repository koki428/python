#熱対流のベクトル表示

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

pngdir='../figs/d004/png/'
datadir="../run/d004/data/"
d=R2D2.R2D2_data(datadir)
for key in d.p:
    exec('%s=%s%s%s' % (key, 'd.p["',key,'"]'))

plt.close("all")

print('input max step number')
nd=input()
nd=int(nd)

#fig=plt.figure(figsize=(15,5))
fig=plt.figure(figsize=(12,6))

n0=0
for n in range(n0,nd):
    print(n)
    t=d.read_time(300+n,silent=True)
    d.read_qq_2d(300+n,silent=True)

    vx=d.q2['vx']
    vy=d.q2['vy']
    """
    vx_sub=vx.reshape(128,8,64,8)
    vy_sub=vy.reshape(128,8,64,8)

    vx_2=vx_sub.mean(axis=(1,3))
    vy_2=vy_sub.mean(axis=(1,3))

    x_2=np.linspace(xmin,xmax,128)
    y_2=np.linspace(ymin,ymax,64)
    """
    ax1=fig.add_subplot(1,1,1,aspect="equal")
    p=ax1.streamplot(y*1.e-8,(x-rsun)*1.e-8,vy,vx,density=6,color=np.sqrt(vx**2+vy**2)*1.e-5,cmap='jet')
    #fig.colorbar(p.lines,label='[km/s]')
    #ax1.set_title('convection vector')
    #ax1.set_xlabel('Y [Mm]')
    #ax1.set_ylabel('X [Mm]')

    if (n == n0):
        fig.tight_layout(pad=0.5)

    plt.pause(0.1)
    plt.savefig(pngdir+'convection_vector_nobar_d6'+'{0:08d}'.format(n)+'.png')

    if (n != nd-1):
        clf()



    

