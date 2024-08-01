import cv2
import numpy as np
import R2D2
import math
import matplotlib.pyplot as plt

plt.close("all")

test=np.zeros((100,100))
for i in range(35,65):
    for j in range(35,65):
        test[i,j]=1000
for i in range(0,10):
    for j in range(0,10):
        test[i,j]=500
        #print(test[i,j])
#print(test)
# 2値化を行う
test_thr=np.zeros((100,100))
#psi_max=np.max(psi)
#psi_min=100
#while True:
    #center=(psi_max+psi_min)*0.5
center=800
for j in range(0,100):
    for i in range(0,100):
        if (test[i,j] > center):
            test_thr[i,j]=110.0
        else:
            test_thr[i,j]=0.0

ret,binary_test = cv2.threshold(test_thr, 50, 255, cv2.THRESH_BINARY)

test_label=binary_test.astype('uint8')
nlabels,labellmages,stats,centroids=cv2.connectedComponentsWithStats(test_label)
    #print("label =",nlabels)
    #if (nlabels-1 > 1):
        #psi_min=center
    #else:
    #    psi_max=center
    #if (abs(psi_max-psi_min) < 1.e+0) and (nlabels-1==1):
    #    break

x=np.linspace(0,100,100)
y=np.linspace(0,100,100)

fig=plt.figure(figsize=(8,10))
ax1=fig.add_subplot(1,2,1,aspect='equal')
ax1=ax1.pcolormesh(x,y,test,cmap='gist_stern')
ax2=fig.add_subplot(1,2,2,aspect='equal')
ax2=ax2.pcolormesh(x,y,binary_test,cmap='gist_stern')
