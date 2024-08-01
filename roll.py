import numpy as np
import math
import sympy as sym
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import R2D2
import cv2
import sys
import os

a=np.arange(12).reshape(3,4)

print(a)

b=np.roll(a,2,axis=1)
print(b)

a=np.roll(a,2,axis=1)
print(a)