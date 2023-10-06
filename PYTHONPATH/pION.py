#!/usr/bin/env python
#
# A new module implementing range-spectrum splitting method to reduce 
# effects of ionospheric contributions 
#
# by Wanpeng Feng, @CCRS/NRCan, 2017-08-25
#
#
from scipy import ndimage
import cv2
import numpy as np
import pSAR
# from scipy.interpolate import interp1d
from scipy.interpolate import griddata
#
#
def ionCorrect(phs,ion,cc=None,theilsen=True,scale=2):
    #
    X,Y = np.meshgrid(np.r_[0:phs.shape[1]],np.r_[0:phs.shape[0]])
    corphs = np.copy(phs)
    #
    # Orbital errors can also be removed with a plane 
    #
    for i in range(3):
       tmpPhs,sim = pSAR.ts.orbitalCor(corphs,X,Y,dem=ion,thirddep='NULL',\
               istheilsen=theilsen,scale=scale,npara=0,cc=cc)
       #
       corphs[corphs!=0] = corphs[corphs!=0] - sim[corphs!=0]
       
    return corphs,ion
#
def Gaussian_kernel(Sx, Sy, sig_x,sig_y):
    #
    # adapted from ISCE-2.1.0
    #
    if np.mod(Sx,2) == 0:
        Sx = Sx + 1

    if np.mod(Sy,2) ==0:
            Sy = Sy + 1
    #
    x,y = np.meshgrid(np.arange(Sx),np.arange(Sy))
    x = x + 1
    y = y + 1
    x0 = (Sx+1)/2
    y0 = (Sy+1)/2
    fx = ((x-x0)**2.)/(2.*sig_x**2.)
    fy = ((y-y0)**2.)/(2.*sig_y**2.)
    k = np.exp(-1.0*(fx+fy))
    a = 1./np.sum(k)
    k = a*k
    #
    return k
    #
def fillzerov2(data):
    #
    X,Y = np.meshgrid(np.r_[0:data.shape[1]],np.r_[0:data.shape[0]])
    points = np.vstack((X[data!=0],Y[data!=0])).T
    values = data[data!=0]
    grid_z1 = griddata(points, values, (X, Y), method='linear')
    return grid_z1
    
def fillzero(data):
    #
    #
    cdata = np.copy(data)
    cdata[cdata==0] = np.nan
    ind = ndimage.distance_transform_edt(np.isnan(cdata),
                                    return_distances=False,
                                    return_indices=True)
    return cdata[tuple(ind)]
    #
def lowpassfilter(data, Sx, Sy, sig_x, sig_y, iters=2,tols=2.0,fill=True):
    #
    # adapted from ISCE-2.1.0
    #
    #
    kernel = Gaussian_kernel(Sx, Sy, sig_x, sig_y) #(800, 800, 15.0, 100.0)
    if fill:
       cdata = fillzero(data)
    else:
       cdata = np.copy(data)
       cdata[cdata==0] = np.min(np.abs(cdata[cdata!=0])) / 1000.
    #
    cdata[cdata==0.] = np.min(np.abs(cdata[cdata!=0]))/1000.
    #
    for i in range(iters):
      #
      filtereddata = cv2.filter2D(cdata,-1,kernel)
      diffda = np.abs((filtereddata - cdata)/cdata)
      # 
      cdata[diffda > tols] = filtereddata[diffda>tols]/10.
      #
    #
    return filtereddata

def unwrapping_errors(phs_L,phs_H,f0,bw,ION,nonION):
    #
    # unwrapping error checking
    # 
    d = np.round(1./(2*np.pi) * \
            (phs_H - phs_L - 2 * bw/(3.*f0) * nonION + 2*bw/(3.*f0) * ION ))
    m = np.round(1./(4*np.pi) * (phs_H + phs_L - 2 * nonION - 2 * ION)-d/2.)
    return d,m
#
def ion_refine(fl,fh,f0,bw,phs_low,phs_hig,ion,ion_non):
    #
    errD,errM = unwrapping_errors(phs_low,phs_hig,f0,bw,ion,ion_non)
    phs_low = phs_low - 2*np.pi*errM
    phs_hig = phs_hig - 2*np.pi*(errM+errD)
    #
    ion = (fl * fh) / (f0 * (fh**2 - fl**2)) * (fh * phs_low - fl * phs_hig)
    ion_non = f0 / (fh**2 - fl**2) * (fh * phs_hig - fl * phs_low)
    #
    return phs_low,phs_hig,ion,ion_non

