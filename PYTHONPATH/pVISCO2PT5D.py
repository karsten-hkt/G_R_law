#!/usr/bin/env python
#
#
import pSAR
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
#
import matplotlib.pyplot as plt
import pSIMP
#
def fpara2topline(fpara,refll=[0,0]):
    #
    outxy = []
    for i in range(fpara.shape[0]):
      #
      cfpara = pSIMP.simp_fpara2whole(fpara[i,:],0)
      poly,deps = pSIMP.simp_corners(cfpara)
      poly,deps = poly[0:3,:],deps[0:3]
      xy = poly[deps == 0,:]
      #
      llxy = pSAR.util.enu2llh_m(xy[:,0]*1000,xy[:,1]*1000,xy[:,0]*0,\
                refll[1],refll[0],0)
      outxy.append(llxy)
    #
    return outxy
#
def fpara2llxy(fpara,refll=[0,0]):
    # fpara in km
    #
    outllxy = []
    for i in range(fpara.shape[0]):
        # below returns coordinates in km
        xy,deps  = pSIMP.simp_corners(fpara[i,:])
        llxy = pSAR.util.enu2llh_m(xy[:,0]*1000,xy[:,1]*1000,xy[:,0]*0,\
                refll[1],refll[0],0)   
        #
        outllxy.append(llxy)
        #
    return outllxy
        #
def llfpara2kmfpara(fpara,refll=[0,0]):
    #
    kmfpara = np.copy(fpara)
    for i in range(fpara.shape[0]):
        #
        enu = pSAR.util.llh2enu(fpara[i,1],fpara[i,0],0,refll[1],refll[0],0)
        kmfpara[i,0] = enu[0]/1000.
        kmfpara[i,1] = enu[1]/1000.
    #
    return kmfpara
#
def vssource2fpara(vssource):
    #
    fpara = np.zeros([vssource.shape[0],10])
    for i in range(vssource.shape[0]):
        #
        cf = vssource[i,:]
        #
        width = (cf[0]-cf[1])/np.sin(cf[2]*np.pi/180)
        #
        strike_slip = np.cos(cf[7]*np.pi/180) * cf[8]
        dip_slip    = np.cos(cf[7]*np.pi/180) * cf[8]
        #
        fpara[i,:] = [cf[4],cf[3],cf[6],cf[2],cf[1],width,cf[5],strike_slip,dip_slip,0]
    #
    return fpara
#
def source2faults(insource):
    #
    fpara = []
    with open(insource,'r') as fid:
        #
        for cline in fid:
            #
           if '# finite fault with # segments' in cline:
               num = fid.readline().split('\n')[0]
               num = int(num)
               # jump two lines 
               fid.readline()
               fid.readline()
               #
               for i in range(num):
                   #
                   c1 = fid.readline().split('\n')[0]
                   c2 = fid.readline().split('\n')[0]
                   # 
                   c1 = [float(tmp) for tmp in c1.split()]
                   c2 = [float(tmp) for tmp in c2.split()]
                   #
                   paras = c1+c2
                   fpara.append(paras)
           #
        #
    #
    return np.array(fpara) 
    #
def visco_source2corners(insource,plot=False):
    # insource = 'source-spherical.param'
    #
    vsfpara = source2faults(insource)
    #
    # vsfpara, [maxdep,mindep,dip,lat,lon,length,strike,rake,slip]
    # fpara, [lon,lat,strike,dip,dep(top),width,length,strike_slip,dip_slip,open_slip]
    #
    fpara = vssource2fpara(vsfpara)
    #
    #
    refll = [fpara[:,0].mean(),fpara[:,1].mean()]
    #
    # kfpara, fault parameters in km
    # reference at lower edge (right) lr
    # 
    kfpara_v = llfpara2kmfpara(fpara,refll=refll)
    #
    # switch reference from right-bottom to top-center
    #
    kfpara   = pSIMP.sim_fparaconv(kfpara_v,99,0)
    #
    # llxy, 4 corners of faults in lonlat
    # topllxy, top-line with depth of 0
    #
    llxy = fpara2llxy(kfpara,refll=refll)
    topllxy = fpara2topline(kfpara,refll=refll)
    #
    #
    if plot:
      for i in range(len(llxy)):
          xy = llxy[i]
          ctop = topllxy[i]
          plt.plot(xy[:,0],xy[:,1],'-b',linewidth=1.5)
          #
          plt.plot(ctop[:,0],ctop[:,1],'-r',linewidth=3.5)
      #
      plta = plt.gca()
      plt.show()
      #
    #
    return llxy,topllxy
