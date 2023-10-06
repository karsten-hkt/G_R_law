#!/usr/bin/env python
#
#
import os
import numpy as np
import sys
import pokada
import matplotlib.pyplot as plt
#
#
if True:
    #
    x = np.linspace(10,50,num=50)
    y = np.linspace(-20,20,num=50)
    #
    pts = [30,0,3,90,-1.5]
    #
    mx,my = np.meshgrid(x,y)
    #
    dE,dN,dU = pokada.dc3dpts(np.reshape(np.array(pts),[1,5]),mx,my,mx*0)
    #
    plt.imshow(dU,extent=[np.min(mx.ravel()),np.max(mx.ravel()),\
                       np.min(my.ravel()),np.max(my.ravel())])
    plt.quiver(mx,my,dE,dN)
    plt.show()
