#! /usr/bin/env python
#
# A new module to operate SLC data. This is intented to use for analysis of SLC
# frequency and spectrum. The first version was developed on 28/08/2017
# This is part of ginsar. pGAMMA and pSAR are another two modules develoepd 
# in ginsar. numpy, scipy, os and matplotlib are external modules managed by 
# annaconda. 
# Develoepd by Wanpng Feng, @CCRS/NRCan. 
# Please return any comments and bugs to wanpeng.feng@hotmail.com.
#
###############################################################################
import numpy as np
import pGAMMA
import pSAR
from scipy import fftpack
import os
import matplotlib.pyplot as plt
#
###############################################################################
def rangesplit_WIN(freqs,fc,bw):
    #
    flag1 = freqs <= fc+bw/2.
    flag2 = freqs >= fc-bw/2.
    flag  = flag1 & flag2
    outW  = freqs * 0.
    outW[flag] = 1.
    return outW

def rangesplit_SLC(inslc,outslc,ref0,rbw,plot=False,fc=None,bw=None,amp=1.,demod=True):
    '''
    inslc,  SLC in the GAMMA format with .par header
    outslc, output SLC in the GAMMA format
    ref0 -0.5,0.5
    rbw  0 - 1.0
    developed by Wanpeng Feng, @CCRS/NRCan, 20170101
    please return comments to wanpeng.feng@hotmail.com
     
    '''
    slc_par   = inslc+'.par'
    freqinfo = pGAMMA.gamma_slcpar2freqs(slc_par)
    dims  = pGAMMA.gamma_slcpar2dim(slc_par)
    f_s   = freqinfo[1]
    t_s   = np.linspace(0,dims[0]-1,num=dims[0]) / f_s
    freqs  = fftpack.fftfreq(dims[0]) * freqinfo[1]
    #
    if fc is None:
       fc = ref0 * freqinfo[1]
    if bw is None:
       bw = freqinfo[1] * rbw
    
    print(" +++++++++++++++++++++++++++++++++++++ ")
    print(" + Start splitting %s at %f with %f " % (inslc,fc,bw))
    print(" +++++++++++++++++++++++++++++++++++++ ")
    filter_WIN = rangesplit_WIN(freqs,fc,bw)
    #
    if plot:
        plt.plot(freqs,filter_WIN,'-r',label='filtering_window')
        plt.xlabel('Frequency in Hertz [Hz]')
        plt.ylabel('Frequency Domain (Spectrum) Magnitude')
        plt.legend()
        plt.show()
    #
    fid = open(inslc,'rb')
    fidout = open(outslc,'wb')
    #
    for i in range(dims[1]):
        if i % 5000 == 0:
           fac = i / float(dims[1]) * 100.
           print(" %-10.3f%% done..." % fac)
        #
        cdata = pSAR.br.br_complex_fid(fid,i+1,dims[0],bigendian=True)
        #
        # filtering 
        #
        cdata = cdata * amp
        x = fftpack.fft(cdata)
        x = x * filter_WIN 
        cdata = fftpack.ifft(x)
        #
        # demodualting or shifting data to a new center
        #
        if demod:
          cdata = demodulate(cdata,fc,t_s)
        pSAR.br.br_complexarray2fid(cdata,fidout,big=True)   
    #
    fid.close()
    fidout.close()
    #
    # Updating par to <outslc>.par
    pGAMMA.gamma_updateF0(inslc+'.par',outslc+'.par',fc,bw=bw)
    #
    #
    if os.path.exists(outslc):
        return fc,bw,True
    else:
        return fc,bw,False
#
def demodulate_SLC(inslc,outslc,threshold=None,fmt='>ff',big=True,f0=None):
    '''
    '''
    slc_par   = inslc+'.par'
    #
    if f0 is None:
       estFactor = gSLC_rngspec(inslc,bigendian=big,threshold=threshold,step=200)
       f0  = np.mean(estFactor[:,2])
    # 
    dims  = pGAMMA.gamma_slcpar2dim(slc_par)
    freqs = pGAMMA.gamma_slcpar2freqs(slc_par)
    f_s   = freqs[1]
    t_s   = np.linspace(0,dims[0]-1,num=dims[0]) / f_s
    #
    fid = open(inslc,'rb')
    fidout = open(outslc,'wb')
    for i in range(dims[1]):
        if i % 5000 == 0:
           fac = i / float(dims[1]) * 100.
           print(" %-10.5f%% done..." % fac)
        #
        cdata = pSAR.br.br_complex_fid(fid,i+1,dims[0],bigendian=True)
        outdata = demodulate(cdata,f0,t_s)
        # pSAR.br.br_complex2fid(outdata,fidout,fmt=fmt)
        pSAR.br.br_complexarray2fid(outdata,fidout,big=True)
    #
    fid.close()
    fidout.close()
    #
    if os.path.exists(outslc):
        return f0,True
    else:
        return f0,False

def demodulate(data,f,ts):
    '''
    To center subband SLC at f, here
    1) f can be fL or fH
    2) ts is configured by
       ts = pixel_index / f_s
       f_s adc_sampling_rate
    
    demoulating kernal was proposed by (Heresh Fattahi, 2017)
    Coded by Wanpeng Feng, @NRCan, 2017-09-06
    '''
    return data * np.exp(-1j * 2 * np.pi * f * ts)
    #
def spectral(cdata,f_s,threshold=None):
    '''
    in_data  complex series in default
    f_s      sampling rate
    '''
    #
    x = fftpack.fft(cdata)
    freqs = fftpack.fftfreq(len(x)) * f_s
    #
    # magnitude of frequency spectra, np.abs(x)
    if threshold is None:
        threshold = np.max(np.abs(x)) * 0.05
    #
    x1     = freqs[np.abs(x) > threshold]
    return x, np.max(x1) - np.min(x1),np.mean([x1])
#
###############################################################################
def gSLC_rngspec(in_slc,bigendian=True,threshold=None,step=200,isprint=False):
    '''
    
    '''
    slc_par  = in_slc+'.par'
    dims     = pGAMMA.gamma_slcpar2dim(slc_par)
    freqs    = pGAMMA.gamma_slcpar2freqs(slc_par)
    f_s      = freqs[1]   # adc_sampling rate
    cen_lines = range(100,dims[1]-101,step)
    outdata   = np.zeros([len(cen_lines),3])
    outdata[:,0] = cen_lines
    #
    fid = open(in_slc,'br')
    for i,cline in enumerate(cen_lines):
      if isprint:
         print(" NO: %d of %s" % (cline,in_slc))
      outdata[i,1],outdata[i,2] = gSLC_bandwidth(fid,sample_lines=50,\
             cen_line=cline,width=dims[0],\
             file_length=dims[1],f_s=f_s,bigendian=True,threshold=threshold)
    fid.close()
    return outdata
    
    
#
def gSLC_bandwidth(in_slc,sample_lines=100,cen_line=None,width=None,\
                   file_length=None,f_s=None,bigendian=True,threshold=None):
    '''
    <gSLC> will be used for a SLC in GAMMA format. <gSLC>.par should be available
    in default with all data information. 
    by FWP, 2017-08-28
    '''
    if width == None:
      slc_par  = in_slc+'.par'
      dims     = pGAMMA.gamma_slcpar2dim(slc_par)
      freqs    = pGAMMA.gamma_slcpar2freqs(slc_par)
      f_s      = freqs[1]   # adc_sampling rate
      width = dims[0]
      file_length = dims[1]
      fid = open(in_slc,'br')
    else:
      fid = in_slc
    #
    if cen_line == None:
       cen_line = int(dims[1]/2)
    #
    l_s = int(cen_line - sample_lines/2)     # starting line
    l_e = int(cen_line + sample_lines/2)     # ending line
    if l_s < 0:
        l_s = 0
    if l_e >= file_length:
        l_e = file_length - 1
    #
    freq_wids = np.zeros([int(l_e - l_s + 1),2])
    counter = 0
    #
    for cl in range(l_s,l_e+1):
        #
        # return a line of complex data from SLC
        cdata = pSAR.br.br_complex_fid(fid,cl,width,bigendian=bigendian)
        x,x_wid,x_mean = spectral(cdata,f_s,threshold=threshold)
        freq_wids[counter,0] = x_wid
        freq_wids[counter,1] = x_mean
        counter += 1
    #
    if width == None:
       fid.close()
    #
    return np.mean(freq_wids[:,0]),np.mean(freq_wids[:,1])
