#hkt change in 230908
#based on the file from Wangpeng

import sys

import pandas as pd

sys.path.append('../src/PYTHONPATH')
import pSAR
import pSIMP
import numpy as np
import area_cut
import pandas as pd
import matplotlib.pyplot as plt


profile_len = 300 #km
nprofile = 1
buffer_width = 30000 #m
faults_file_line = '../data/cross_section_line.inp'
faults_file = '../data/Prof_test_0.inp'
#using Prof_test_0 to generate a buffer parallel the fault strike
aftershock_file = '../data/data_processing_7_8.csv'
faults = np.loadtxt(faults_file)
faults_line = np.loadtxt(faults_file_line)
aftershocks  = pd.read_csv(aftershock_file)
aftershocks  = np.array(aftershocks)
faults[faults[:,0]>180,0] = faults[faults[:,0]>180,0] - 360 #if the data using 360 need to be tranlated
#the faults need to be changed into the biggest one and the smallest one
ind = [0,-1]
ft = faults[ind,:]
#using the utm projection
'''
zone utm zone number
nl   utm zone letter
dip  dip angle of fault
dep  burial depth
'''
zone = None
nl = None
cfpara,zone,nl = pSIMP.simp_2pts2fault(ft[:,0],ft[:,1],dip=85.,wid=10,dep=0.,zone=zone,nl=nl)
#print(cfpara.shape)
xm = (ft[0,0]+ft[1,0])/2#np.linspace(ft[0,0],ft[1,0],num=nprofile)
ym = (ft[0,1]+ft[1,1])/2#np.linspace(ft[0,1],ft[1,1],num=nprofile)
fpara = np.zeros([1,10])
outdata = []
for i in range(1):
    ux,uy,a,b = pSAR.utm_conversion.from_latlon(ym,xm,force_zone_number=zone,force_zone_letter=nl)#ym[i],xm[i],force_zone_number=zone,force_zone_letter=nl)
    ux,uy = ux/1000,uy/1000
    #print(ux)
    #print(uy)
    tfpara = [ux,uy,cfpara[0,2]+90,89.999,0,10,profile_len,0,0,0]
    fpara[0,:] = tfpara
    cpolygon = pSIMP.simp_llcorners(fpara[0,:],zone,nl)
    outdata.append(cpolygon[[0,1],:])
    cdata = cpolygon[[0,1],:]
    np.savetxt('../out_file/buffer_generate/Prof_test_%d.inp' % i,cdata,fmt='%f %f',newline='\n')
plt.scatter(aftershocks[:,2],aftershocks[:,3])#,s = aftershocks[:,5]/10)
#plt.plot(ft[:,0],ft[:,1],linewidth=5.0,color='red')
plt.plot(faults_line[:,0],faults_line[:,1],linewidth=5.0,color='red')
plt.plot(xm,ym,'ob',markersize=1)
#build the buffer
for i in range(len(outdata)):
    print('start to select buffer %d' %i)
    cdata = outdata[i]
    ux,uy = pSAR.utm_conversion.lltoutm(cdata[:,1],cdata[:,0],force_zone_number=zone,force_zone_letter=nl)
    prof = np.vstack((ux,uy)).T
    x,y = pSAR.opt.line2region(prof,buf=buffer_width)
    plat,plon = pSAR.utm_conversion.utmtoll(x,y,zone,force_zone_letter=nl)
    plt.plot(cdata[:,0],cdata[:,1],'-g')
    plt.plot(plon,plat,'-y')
    #change to 2 dimential
    out_buffer=np.zeros((plon.size,2))
    out_buffer[:,0]=plon
    out_buffer[:,1]=plat
    out_aftershocks=[]
    # jugde if the aftershocks inside the buffer
    for j in range(aftershocks.shape[0]):
        point = np.zeros(2)
        point[0] = aftershocks[j,2]
        point[1] = aftershocks[j,3]
        if area_cut.is_in_poly(point,out_buffer):
            out_aftershocks.append(aftershocks[j,:])
    np.savetxt('../out_file/b_variation_after_7_8/buffer_test_%d.txt' % i,out_buffer,fmt='%f,%f',newline='\n')
    data = pd.DataFrame(out_aftershocks)
    data.to_csv('../out_file/b_variation_after_7_8/aftershocks_in_buffer_%d.csv' % i,index=0,header=0)
    print('success')
plt.savefig('../out_file/b_variation_after_7_8/buffer_shape.png',bbox_inches='tight')
