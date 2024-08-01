#各時間での画像から動画を作成
from PIL import Image
import os
import numpy as np
import cv2
import R2D2
import sys

CLIP_FPS=10.0

# print("input caseid id (3 digit)")
caseid = 0
caseid = sys.argv[1]
caseid = "d"+caseid.zfill(3)
print(caseid)
    
datadir="../run/"+caseid+"/data/"
pngdir="../figs/"+caseid+"/png/"
movdir="../figs/"+caseid+"/mov/"
os.makedirs(movdir,exist_ok=True)

png_name='py_bz_se'
mov_name='py_bz_se'

filepath_mp4=movdir+caseid+mov_name+'.mp4'
filepath_gif=movdir+caseid+mov_name+'.gif'


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

# print('input end step number')
# nd=input()
# nd=int(nd)

picture=[]

#f=open(data_dir+'num.dat','r')
#num=f.read().split()
#n=int(num[0])

img_mp4=cv2.imread(pngdir+png_name+'{0:08d}'.format(n0)+'.png')

w=img_mp4.shape[1]
h=img_mp4.shape[0]
codec=cv2.VideoWriter_fourcc(*'mp4v')
video=cv2.VideoWriter(filepath_mp4,codec,CLIP_FPS,(w,h))

for idx in range(n0,nd+1):
    filename=pngdir+png_name+'{0:08d}'.format(idx)+'.png'
    img_mp4=cv2.imread(filename)
    img_gif=Image.open(filename)

    video.write(img_mp4)
    picture.append(img_gif)

video.release()
picture[0].save(filepath_gif,save_all=True,append_images=picture[1:],optimize=True,duration=500,loop=0)
