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

dir="../run/"+caseid+"/data/"
casedir="../figs/"+caseid
os.makedirs(casedir,exist_ok=True)

R2D2.read_init(dir,"3d")
for key in R2D2.p:
    exec('%s = %s%s%s' % (key, 'R2D2.p["',key,'"]'))

try:
    n0
except NameError:
    n0 = 0
if  n0 > R2D2.p["nd"]:
    n0 = R2D2.p["nd"]

