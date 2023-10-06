#!/usr/bin/env python
#
# a module to drive GAMMA to process SAR data in Python
# Developed by Wanpeng Feng, @NRCan, 2016-04-15
#
from __future__ import division
#
import re
import sys
import os
import numpy as np
import subprocess
import datetime
import glob
import xml.dom.minidom
import shutil
#
# Below are my personal python-based modules
#
import pDATA
import pALOS
import pSAR
# 
# Make sure to run this smoothly in HPC
# Modified by Wanpeng Feng, @NRCan, 2017-07-11
#
try:
    from shapely.geometry import Polygon
    from shapely.ops import cascaded_union
except:
    1
    #
def gamma_asarN1toSLC(N1,outslc):
    # Create SLC 
    # from a Level 1.1 N1 ENVSAT data
    # using par_ASAR 
    goSLC='par_ASAR %s %s' % (N1,outslc)
    flag,info,errs = gamma_run(goSLC)
    return flag,info,errs
    
###############################################################################
def gamma_slconroi(inslcpar,roi):
    #
    numROI = [float(croi) for croi in roi.split(',')]
    if numROI[0] == numROI[1] or numROI[2] == numROI[3]:
        return True
    else:
        # in default in lon,lat
        polygonROI = pSAR.roipac.ext2polygon(numROI)
        #
        # corners of SLC in lat,lon
        pts = gamma_slc_corners(inslcpar)
        polygonSLC = np.zeros([5,2])
        polygonSLC[0:4,0] = pts[:,1]
        polygonSLC[0:4,1] = pts[:,0]
        polygonSLC[4,:]   = polygonSLC[0,:]
        #
        polyROI           = Polygon(polygonROI)
        polySLC           = Polygon(polygonSLC)
        flag = polyROI.intersects(polySLC)
        return flag
#
###############################################################################
#
def gamma_cfg(incfg):
    '''
    Read gInSAR processing configure files into a dict variable
    keywords show the physical meanings.
    #
    Adding some information by Wanpeng Feng, @CCRS/NRCan, 2018-01-26
    
    '''
    #
    info = {}
    with open(incfg,'r') as fid:
        for cline in fid:
            #
            cline = pDATA.bytestoutf8(cline)
            if cline[0] != "*":
              cline = cline.split('\n')[0].split()
              if len(cline) > 1:
                 info[cline[0]] = cline[1:]
              else:
                 info[cline[0]] = ['']
    return info
#
###############################################################################
#            
def gamma_cfg2MLN(incfg):
    #
    info = gamma_cfg(incfg)
    rlk  = info['ISAR_RLOOKS'][0]
    azlk = info['ISAR_ZLOOKS'][0]
    #
    if len(rlk) == 0 or len(azlk) == 0:
       postrng = float(info['ISAR_RNG_POS'][0])
       postazi = float(info['ISAR_AZI_POS'][0])
       slcs = glob.glob('SLC/*/*.slc')
       if len(slcs) == 0:
           print(" Error: you need to work in a gInSAR working folder.")
           sys.exit(-1)
       parinfo = gamma_slcpar(slcs[0]+'.par')
       rng_pixel_size = float(parinfo['range_pixel_spacing'][0])
       azi_pixel_size = float(parinfo['azimuth_pixel_spacing'][0]) 
       rlk = int(postrng/rng_pixel_size)
       azlk= int(postazi/azi_pixel_size)
       #
    return rlk,azlk
#
###############################################################################
#
def gamma_updates1rsc(inrsc,mpar,spar):
    #
    info,ext = pSAR.roipac.rsc_read(inrsc)
    sensor1 = gamma_s1parTOsensor(mpar)
    sensor2 = gamma_s1parTOsensor(spar)
    time1 = gamma_s1parTOtime(mpar)
    time2 = gamma_s1parTOtime(spar)
    #
    info['MASTER'] = time1.replace('-','')[0:8]
    info['SLAVE'] = time2.replace('-','')[0:8]
    info['TEMPB'] = pSAR.ts.diff2dates(info['MASTER'],info['SLAVE'])
    info['MTIME'] = time1
    info['STIME'] = time2
    info['POL_MASTER'] = gamma_s1parTOpol(mpar)
    info['POL_SLAVE']  = gamma_s1parTOpol(spar)
    info['SENSOR_M'] = sensor1
    info['SENSOR_S'] = sensor2
    pSAR.roipac.info_to_rsc(info,inrsc)
    return True
#
###############################################################################
def gamma_s1parTOpol(inslcpar):
    info = gamma_slcpar(inslcpar)
    tmp = info['sensor']
    if len(tmp)>3:
        return tmp[3]
    else:
        return 'NONE'
#    
def gamma_s1parTOsensor(inslcpar):
    info = gamma_slcpar(inslcpar)
    tmp = info['title'][0][0:3]
    return tmp
#
###############################################################################
#
def gamma_s1parTOtime(inslcpar):
    #
    info = gamma_slcpar(inslcpar)
    tmp = info['date']
    mm = tmp[1]
    dd = tmp[2]
    hh = tmp[3]
    M = tmp[4]
    if len(mm) < 2:
        mm = '0%s' % tmp[1]
    if len(dd) < 2:
        dd = '0%s' % tmp[2]
    if len(hh) < 2:
        hh = '0%s' % tmp[3]
    if len(M) < 2:
        M = '0%s' % tmp[4]   
    times = '%s-%s-%sT%s:%s:%sZ' % (tmp[0],mm,dd,hh,M,tmp[5])
    return times
#
def gamma_slcpar2gmt(inslcpar,outgmt):
    #
    # return an ascii file which saves four corners of SLC data 
    #
    pts = gamma_slc_corners(inslcpar)
    pts = np.vstack((pts,pts[0,:]))
    np.savetxt(outgmt,pts,fmt='%f %f',newline='\n')
    #
    if os.path.exists(outgmt):
        return True
    else:
        return False
#
###############################################################################
#
def gamma_slcpar2kml(inslcpar,outkml):
    #
    # Return 4 corners from a SLC_par
    #
    pts = gamma_slc_corners(inslcpar)
    pts = np.vstack((pts,pts[0,:]))
    #
    # generate a kml structure from pDATA
    #
    kmlstr_s,kmlstr_e,polystr_s,polystr_e = pDATA.kml_poly()
    fidkml = open(outkml,'w')
    fidkml.write(kmlstr_s % outkml)
    #
    parinfo = gamma_slcpar(inslcpar)
    keys = parinfo.keys()
    descri = ["%s:%s" % (key,parinfo[key]) for key in keys]
    descri = '\n'.join(descri)
    fidkml.write(polystr_s % (os.path.basename(inslcpar),descri))
    #
    for ind in range(pts.shape[0]):
        outloc = ("                %f,%f,0\n" % (pts[ind,1],pts[ind,0]))
        fidkml.write(outloc)
    #
    fidkml.write(polystr_e)
    fidkml.write(kmlstr_e)
    fidkml.close()
    return True
#
###############################################################################
#
def gamma_SLC_intf(slc1,slc2,rlk,azlk,outint,offpar=None):
#
    slc1par = slc1+'.par'
    slc2par = slc2+'.par'
    log     = outint+'.log'
    if offpar is None:
       #
       offpar = outint+'.off'
       coffpar = gamma_create_offset(slc1par,slc2par,offpar,algorithm=1,\
                            rlks=rlk,azlks=azlk,iflg=0)
       gamma_run(coffpar)
    #
    # print(slc1,slc2,slc1par,slc2par,offpar,outint,str(rlk),str(azlk))
    intf = 'SLC_intf %s %s %s %s %s %s %s %s 0 - 0 0 0 0' % \
          (slc1,slc2,slc1par,slc2par,offpar,outint,str(rlk),str(azlk))
    if not os.path.exists(outint):
        fidlog = open(log,'a')
        gamma_run(intf,log=fidlog)
        fidlog.close()
    return True
#
###############################################################################
#
def gamma_freqshift(slcpar1,slcpar2):
    #
    freqs_1 = gamma_slcpar2freqs(slcpar1)
    freqs_2 = gamma_slcpar2freqs(slcpar2)
    return freqs_2[0]-freqs_1[0],freqs_2[2]
#
###############################################################################
#
def gamma_bpf(inslc,outslc,fa,bw,beta=1.0,fir_len=128):
    '''
    A function was used in ginsar_IONest.py to correct for ionospheric signals
    by Wanpeng Feng, @CCRS/NRCan, 2017-09-11
    '''
    par = inslc+'.par'
    info = gamma_slcpar(par)
    width = int(info['range_samples'][0])
    #freqs = gamma_slcpar2freqs(par)
    #
    #
    bpf = ' bpf %s %s %d %f %f 0.0 1.0 0 0 - - 0 0 %f %d' % \
          (inslc,outslc,width,fa,bw,beta,fir_len)
    if not os.path.exists(outslc):
       print(' + %s' % bpf)
       os.system(bpf)
    #
    outpar = outslc+'.par'
    if not os.path.exists(outpar):
       #
       # shutil.copy(par,outpar)
       info = gamma_slcpar(par)
       #tfc = float(info['radar_frequency'][0])
       #tbw = float(info['chirp_bandwidth'][0])
       shutil.copy(par,outpar)
       #
       # radar_frequency and chirp_bandwidth will be re-estimated based on data
       #
       #cbw,rf0 = pSLC.gSLC_bandwidth(outslc)
       #
    #
    #
    return True

def gamma_slcpar2dir(inpar):
    #
    info = gamma_slcpar(inpar)
    azi  = float(info["azimuth_angle"][0])
    if azi > 0:
        return 'Right'
    else:
        return 'Left'
#
###############################################################################
def gamma_updateF0(inpar,outpar,fshift,bw=None):
    #
    freqs = gamma_slcpar2freqs(inpar)
    f0 = freqs[0] + fshift
    set_value = 'set_value %s %s %s "%e Hz"' % \
                 (inpar,outpar,"radar_frequency",f0)
    print(' pGAMMA: %s' % set_value)
    os.system(set_value)
    #
    if bw is not None:
       set_value = 'set_value %s %s %s "%e Hz"' % \
                 (outpar,outpar,"chirp_bandwidth",bw)
       print(' pGAMMA: %s' % set_value) 
       os.system(set_value)
       #
    return True
#       
def gamma_set_value(inpar,outpar,keyword,new_value):
    #
    '''
    To replace a value in a GAMMA par file.
    '''
    set_value=('set_value %s %s %s "%s"' % (inpar,outpar,keyword,new_value))
    flag,info,errs = gamma_run(set_value)
    #
    return flag
#
###############################################################################
#
def gamma_slc_ovr(inslc,outslc,ovr=2):
    '''
    To resample a SLC by an integer sample value
    '''
    inslc_par  = inslc+'.par'
    outslc_par = outslc+'.par'
    SLC_ovr = 'SLC_ovr %s %s %s %s %d' % (inslc,inslc_par,outslc,\
                                          outslc_par,ovr)
    #
    flag,prints,errs = gamma_run(SLC_ovr)
    if flag == 0:
        return True
    else:
        return False
###############################################################################   
def alosCEOS2SLCpar(ledfile,imagefile,pol='HH'):
    #
    '''
    To focuse an ALOS CEOS raw file
    '''
    info = pALOS.led2info(ledfile)
    outname = info['yyyymmdd']+'_'+pol
    msppar  = outname+'.msp.par'
    procpar = outname+'.proc.par'
    slc_par = outname+'.slc.par'
    cant    = outname+'.ant'
    log     = outname+'.log'
    #
    gamma_home = os.environ["GAMMA_HOME"]
    ant = gamma_home+'/MSP/sensors/palsar_ant_20061024.dat'
    ##############################################################
    #
    if not os.path.exists(slc_par):
      #
      fidlog  = open(log,'w')     
      PALSAR_proc = 'PALSAR_proc %s %s %s %s -' % (ledfile,msppar,procpar,\
                                    imagefile)
      flag,prints,errs = gamma_run(PALSAR_proc,log=fidlog)
      PALSAR_antpat = 'PALSAR_antpat %s %s %s %s' % (msppar,procpar,ant,cant)
      flag,prints,errs = gamma_run(PALSAR_antpat,log=fidlog)
      #
      par_MSP= 'par_MSP %s %s %s' % (msppar,procpar,slc_par)
      flag,prints,errs = gamma_run(par_MSP,log=fidlog)
      fidlog.close()
    if os.path.exists(slc_par):
        return True
    else:
        return False
###############################################################################
def alosCEOS2SLC_WS(ledfile,imagefile,pol='HH'):
    #
    '''
    To focuse a CEOS raw file
    '''
    info = pALOS.led2info(ledfile)
    for ni in range(5):
      outname = info['yyyymmdd']+'_'+pol+'W'+str(int(ni+1))
      msppar  = outname+'.msp.par'
      procpar = outname+'.proc.par'
      slc     = outname+'.slc'
      slc_par = outname+'.slc.par'
      raw     = outname+'.raw'
      cant    = outname+'.ant'
      azsp    = outname+'.azsp'
      rspec   = outname+'.rspec'
      rc      = outname+'.rc'
      autof   = outname+'.autof'
      log     = outname+'.log'
      mlcc    = outname+'.mlcc'
      dop     = outname+'.dop'
      #
      gamma_home = os.environ["GAMMA_HOME"]
      ant = gamma_home+'/MSP/sensors/palsar_ant_20061024.dat'
      ##############################################################
      #
      if not os.path.exists(slc_par) or \
         not os.path.exists(slc):
        #
        fidlog  = open(log,'w')     
        PALSAR_proc = 'PALSAR_proc_WB %s %s %s %s %d %s - 0' % (ledfile,msppar,procpar,\
                                      imagefile,ni+1,raw)
        flag,prints,errs = gamma_run(PALSAR_proc,log=fidlog)
        PALSAR_antpat = 'PALSAR_antpat %s %s %s %s' % (msppar,procpar,ant,cant)
        flag,prints,errs = gamma_run(PALSAR_antpat,log=fidlog)
        if flag != 0:
           print(prints)
        dop_mlcc = 'dop_mlcc %s %s %s %s' % (msppar,procpar,raw,mlcc)
        flag,prints,errs = gamma_run(dop_mlcc,log=fidlog)
        #
        doppler = 'doppler %s %s %s %s' % (msppar,procpar,raw,dop)
        flag,prints,errs = gamma_run(doppler,log=fidlog)
        #
        azsp_IQ = 'azsp_IQ %s %s %s %s' % (msppar,procpar,raw,azsp)
        flag,prints,errs = gamma_run(azsp_IQ,log=fidlog)
        if flag != 0:
           print(prints)
        #   
        rspec_IQ = 'rspec_IQ %s %s %s %s' % (msppar,procpar,raw,rspec)
        flag,prints,errs = gamma_run(rspec_IQ,log=fidlog)
        if flag != 0:
           print(prints)
        #   
        pre_rc = 'pre_rc %s %s %s %s' % (msppar,procpar,raw,rc)
        flag,prints,errs = gamma_run(pre_rc,log=fidlog)
        if flag != 0:
           print(prints)
        #  
        autof  = 'autof %s %s %s %s ' % (msppar,procpar,rc,autof)
        flag,prints,errs = gamma_run(autof,log=fidlog)
        flag,prints,errs = gamma_run(autof,log=fidlog)
        azproc = 'az_proc %s %s %s %s %d 0' % (msppar,procpar,rc,slc,8192)
        flag,prints,errs = gamma_run(azproc,log=fidlog)
        if flag != 0:
           print(prints)
        par_MSP= 'par_MSP %s %s %s' % (msppar,procpar,slc_par)
        flag,prints,errs = gamma_run(par_MSP,log=fidlog)
        fidlog.close()
        #
        #
    return True 
        
def alosCEOS2SLC(ledfile,imagefile,pol='HH'):
    '''
    To focuse a CEOS raw file
    '''
    info = pALOS.led2info(ledfile)
    outname = info['yyyymmdd']+'_'+pol
    msppar  = outname+'.msp.par'
    procpar = outname+'.proc.par'
    slc     = outname+'.slc'
    slc_par = outname+'.slc.par'
    raw     = outname+'.raw'
    cant    = outname+'.ant'
    azsp    = outname+'.azsp'
    rspec   = outname+'.rspec'
    rc      = outname+'.rc'
    autof   = outname+'.autof'
    log     = outname+'.log'
    mlcc    = outname+'.mlcc'
    dop     = outname+'.dop'
    #
    gamma_home = os.environ["GAMMA_HOME"]
    ant = gamma_home+'/MSP/sensors/palsar_ant_20061024.dat'
    ##############################################################
    #
    if not os.path.exists(slc_par) or \
       not os.path.exists(slc):
      #
      fidlog  = open(log,'w')     
      PALSAR_proc = 'PALSAR_proc %s %s %s %s %s' % (ledfile,msppar,procpar,\
                                    imagefile,raw)
      flag,prints,errs = gamma_run(PALSAR_proc,log=fidlog)
      PALSAR_antpat = 'PALSAR_antpat %s %s %s %s' % (msppar,procpar,ant,cant)
      flag,prints,errs = gamma_run(PALSAR_antpat,log=fidlog)
      if flag != 0:
         print(prints)
      dop_mlcc = 'dop_mlcc %s %s %s %s' % (msppar,procpar,raw,mlcc)
      flag,prints,errs = gamma_run(dop_mlcc,log=fidlog)
      #
      doppler = 'doppler %s %s %s %s' % (msppar,procpar,raw,dop)
      flag,prints,errs = gamma_run(doppler,log=fidlog)
      #
      azsp_IQ = 'azsp_IQ %s %s %s %s' % (msppar,procpar,raw,azsp)
      flag,prints,errs = gamma_run(azsp_IQ,log=fidlog)
      if flag != 0:
         print(prints)
      #   
      rspec_IQ = 'rspec_IQ %s %s %s %s' % (msppar,procpar,raw,rspec)
      flag,prints,errs = gamma_run(rspec_IQ,log=fidlog)
      if flag != 0:
         print(prints)
      #   
      pre_rc = 'pre_rc %s %s %s %s' % (msppar,procpar,raw,rc)
      flag,prints,errs = gamma_run(pre_rc,log=fidlog)
      if flag != 0:
         print(prints)
      #  
      autof  = 'autof %s %s %s %s ' % (msppar,procpar,rc,autof)
      flag,prints,errs = gamma_run(autof,log=fidlog)
      flag,prints,errs = gamma_run(autof,log=fidlog)
      azproc = 'az_proc %s %s %s %s %d 0' % (msppar,procpar,rc,slc,8192)
      flag,prints,errs = gamma_run(azproc,log=fidlog)
      if flag != 0:
         print(prints)
      par_MSP= 'par_MSP %s %s %s' % (msppar,procpar,slc_par)
      flag,prints,errs = gamma_run(par_MSP,log=fidlog)
      fidlog.close()
      #
      parinfo = gamma_slcpar(slc_par)
      chirpwidth = float(parinfo["chirp_bandwidth"][0])
      #
      # 2.8e7 for fbs, chirp_bandwidth for a single-pol PalSAR data
      # 
      fbf_chirpwidth = 2.8e7
      ovr = round(fbf_chirpwidth/chirpwidth)
      #
      if ovr > 1:
         #
         gamma_slc_ovr(slc,slc+'.ovr',ovr=ovr)
         new_chirp = chirpwidth * ovr
         gamma_set_value(slc+'.ovr.par',slc+'.ovr.par','chirp_bandwidth',\
                         ("%s Hz" % new_chirp))
         shutil.move(slc+'.ovr',    slc)
         shutil.move(slc+'.ovr.par', slc_par)
      #
    if os.path.exists(slc) and os.path.exists(slc_par):
        return slc
    else:
        return False
#
###############################################################################
#    
def gamma_checkslcexists(slc):
    '''
    To check if <SLC> and and <SLC>.par exist...
    True,  yes 
    False, no
    '''
    flag1 = os.path.exists(slc)
    flag2 = os.path.exists(slc+'.par')
    if (not flag1 or not flag2):
        return False
    else:
        return True
#
###############################################################################
#
def gamma_s1coregquality(in_quality_sta):
    #
    offsets = []
    with open(in_quality_sta,'r') as fid:
        #counter = 0
        for cline in fid:
            cline = pSAR.util.bytestoutf8(cline)
            if "azimuth_pixel_offset" in cline:
                cline = cline.split('\n')[0]
                cline = cline.split()
                offsets.append(cline[1])
    #
    if len(offsets)>0:
        return offsets[-1]
    else:
        return None
#
def gamma_rmspace(instring):
    return instring.replace(" ","")
#
def gamma_bursttab(b1,b2,burst_tab):
    #
    fid = open(burst_tab,'w')
    fid.write("%d %d\n" % (b1,b2))
    fid.close()
    if os.path.exists(burst_tab):
        return True
    else:
        return False
#        
def gamma_multi_look(in_slc,out_mli,rlook,azlook):
    '''
    To generate a multilooking image from a SLC
    
    '''
    in_slcpar = in_slc+'.par'
    out_mlipar = out_mli+'.par'
    sTr='multi_look %s %s %s %s %s %s 0' % (in_slc,in_slcpar,out_mli,\
                                            out_mlipar,str(rlook),str(azlook))
    flag,prints,errs = gamma_run(sTr)
    return flag
###############################################################################
#
def gamma_tops2tab(in_slc,in_slc_par,in_toppar,outtab):
    #
    fid = open(outtab,'w')
    fid.write("%s %s %s\n" % (os.path.abspath(in_slc),\
                              os.path.abspath(in_slc_par),\
                              os.path.abspath(in_toppar)))
    fid.close()
    if os.path.exists(outtab):
        return True
    else:
        return False
        
###############################################################################        
def gamma_s1extractburst(in_slc,in_slc_par,in_toppar,outdir,b1,b2):
    # 
    # to extract burst from a TOPS slc
    #
    indate  = os.path.basename(in_slc)[0:8]
    pol = ''
    if len(os.path.basename(in_slc).split('.')[0].split('_')) >= 3:
       pol = os.path.basename(in_slc).split('.')[0].split('_')[2]
       
    #
    in_tab  = os.path.join(outdir,indate+'_in')
    out_tab = os.path.join(outdir,indate+'_out')
    #
    iw = os.path.basename(in_slc).split('.')[0].split('_')[1]
    flag_in = gamma_tops2tab(in_slc,in_slc_par,in_toppar,in_tab)
    #
    SLC_dir = os.path.join(outdir,'SLC')
    #
    if not os.path.exists(SLC_dir):
        os.makedirs(SLC_dir)
    #
    SLC_dir_date = os.path.join(SLC_dir,indate)
    if not os.path.exists(SLC_dir_date):
        os.makedirs(SLC_dir_date)
    #
    out_slc = os.path.join(SLC_dir_date,indate+"_"+iw+'_'+pol+'.slc')
    out_slc_par = out_slc+'.par'
    out_toppar  = os.path.join(SLC_dir_date,indate+'_'+iw+"_"+pol)+'.tops_par'
    #
    flag_out= gamma_tops2tab(out_slc,out_slc_par,out_toppar,out_tab)
    #
    if flag_in and flag_out:
       #
       out_burst_tab = os.path.join(outdir,indate+"_bursts")
       flag_bursts   = gamma_bursttab(b1,b2,out_burst_tab)
       if flag_bursts and not os.path.exists(out_slc):
          #
          SLC_copy_TOPS = ("SLC_copy_S1_TOPS %s %s %s" % \
                           (in_tab,out_tab,out_burst_tab)) 
          os.system(SLC_copy_TOPS)
    if os.path.exists(out_slc):
       return True
    else:
       return False
#
###############################################################################
#       
def ginsarsta2info(in_sta):
    '''
    Return metadata for interferograms generated by gInSAR
    ---
    Updated by Wanpeng Feng, @CCRS/NRCan, 2017-04-10
    
    '''
    # updated by Wanpeng Feng, @RNCan, 2017-04-10
    #
    sta_dir     = os.path.dirname(os.path.abspath(in_sta))
    sta_slc_par = glob.glob(sta_dir+'/*.slc.par')
    #
    if len(sta_slc_par)>0:
       #
       slc_info = gamma_slcpar(sta_slc_par[0])
       azidegree= float(slc_info["azimuth_angle"][0])
       if azidegree == 90.:
          pointing = "right"
       else:
          pointing = "left"  
    else:
       pointing = "NULL"
    ############################################################
    #
    fid = open(in_sta,'r')
    #
    sensor = "SAR"
    inc    = 0.
    pbase  = 0.
    heading = 0.
    wavelength = 0.
    unit   = 'NULL'
    flags  = {}
    flags["SENSOR"]       = 0
    flags["INC"]          = 0
    flags["PBASE"]        = 0
    flags["WAVELENGTH"]   = 0
    flags["UNIT"]         = 0
    flags["HEADING"]      = 0
    flags["Ann_Pointing"] = pointing
    #
    ###########################################################################
    #
    for cline in fid:
        cline = cline.split('\n')[0]
        if ("Sensor" in cline and flags["SENSOR"] == 0 ):
            sensor = cline.split(':')[1]
            #
            flags["SENSOR"] = 1
            #
        if ("Incidence" in cline and flags["INC"]==0 ):
            inc = float(cline.split(':')[1])
            flags["INC"] = 1
            #
        if ("Heading" in cline and flags["HEADING"]==0 ):
            heading = float(cline.split(':')[1])
            flags["HEADING"] = 1
            #
        if ("pBase" in cline and flags["PBASE"]==0):
            pbase=float(cline.split(':')[1][:-1])
            flags["PBASE"] = 1
        if ("Wavelength" in cline and flags["WAVELENGTH"]==0):
            wavelength = float(cline.split(':')[1])
            flags["WAVELENGTH"] = 1
        if ("Unit" in cline and flags["UNIT"]==0):
            unit = cline.split(':')[1]
            flags["UNIT"] = 1
    #
    return gamma_rmspace(sensor),inc,heading,pbase,\
           wavelength,gamma_rmspace(unit),pointing
###############################################################################        
    
def gamma_data2geotiff(dempar,ingammadata,outdata,intype=2,nodata=0.0):
    #
    if not os.path.exists(dempar):
        print("%s cannot be found. Check it before to rerun." % dempar)
        return False
    if not os.path.exists(ingammadata):
        print("%s cannot be found. Check it before to rerun." % ingammadata)
        return False
    #
    gamma_command_STR=("data2geotiff %s %s %d %s %f " % \
                       (dempar,ingammadata,intype,outdata,nodata))
    # print(gamma_command_STR)
    flag,info,error = gamma_run(gamma_command_STR)
    #
    if flag != 0:
       #print(info)
       return False
    else:
       return True
###############################################################################
#    
class strStruc(object):
    def __init__(self,cID,slc,slc_par):
        self.ID=cID
        self.SLC=slc
        self.SLC_PAR=slc_par
        slc_date=os.path.basename(self.SLC)
        self.DATE=slc_date[0:8]
#
class strOffsets(object):
    def __init__(self,no,mno,sno,mdate,sdate,offsets):
        self.ID=no
        self.MID=mno
        self.SID=sno
        self.MDATE=mdate
        self.SDATE=sdate
        self.OFFSETS=offsets
#
###############################################################################
def gamma_s1par2mission(s1_slc_par):
    #
    mission = "S1A"
    fid = open(s1_slc_par,'r')
    for cline in fid:
        if "title:" in cline:
            tmp = cline.split()[1].split('-')
            mission = tmp[0].upper()
    return mission            
#
def gamma_s1par2time(s1_slc_par):
    #
    # return time of S1-TOPS SLC data
    #
    avg_st = 0.
    st1    = 0.
    st2    = 0.
    #
    fid = open(s1_slc_par,'r')
    for cline in fid:
        if "title:" in cline:
            tmp = cline.split()[1].split('-')
            st1 = pSAR.ts.dates2jd(tmp[4])
            st2 = pSAR.ts.dates2jd(tmp[5])
            avg_st = (st1 + st2) / 2.
    return avg_st,st1,st2
###############################################################################    
def gamma_diff_par(in_diff_par):
    #
    if os.path.exists(in_diff_par):   
       fid = open(in_diff_par,'r')
       for cline in fid:
           if "range_samp_1" in cline:
               cline = cline.split('\n')[0]
               tmp   = cline.split()
               width = int(tmp[1])
           if "az_samp_1" in cline:
               cline = cline.split('\n')[0]
               tmp   = cline.split()
               length = int(tmp[1])
    else:
        width,length = None,None
    #
    return width,length
    #
###############################################################################
#
def gamma_off_par(in_off_par):
    #
    if os.path.exists(in_off_par):   
       fid = open(in_off_par,'r')
       for cline in fid:
           if "offset_estimation_range_samples" in cline:
               cline = cline.split('\n')[0]
               tmp   = cline.split()
               width = int(tmp[1])
           if "offset_estimation_azimuth_samples" in cline:
               cline = cline.split('\n')[0]
               tmp   = cline.split()
               length = int(tmp[1])
    else:
        width,length = None,None
    #
    return width,length
#
############################################################################### 
#
def gamma_asar_orb(slcpar,vor,nstat=64):
    #
    doris_vec = 'DORIS_vec %s %s %d' % (slcpar,vor,nstat)
    flag,info,errs = gamma_run(doris_vec)
    return flag,info,errs
#     
###############################################################################
#
def gamma_dempar2roi(in_dem_par):
    #
    try:
      fid = open(in_dem_par,'r')
      for cline in fid:
          cline =cline.split('\r')[0]
          tmp = cline.split()
          if "width" in cline:
              width = int(tmp[1])
          if "nlines" in cline:
              file_length = int(tmp[1])
          if "corner_lat" in cline:
              clat= np.float64(tmp[1])
          if "corner_lon" in cline:
              clon= np.float64(tmp[1])
          if "post_lat" in cline:
              post_lat= np.float64(tmp[1])
          if "post_lon" in cline:
              post_lon= np.float64(tmp[1])
      #
      maxlon = (width-1) * post_lon + clon
      minlat = (file_length - 1) * post_lat + clat
      roilon = np.array([clon,clon,maxlon,maxlon,clon],dtype='float64')
      roilat = np.array([clat,minlat,minlat,clat,clat],dtype='float64')
      roi = np.vstack((roilon,roilat)).T 
      
      return roi
    except:
      return None
############################################################################### 
def s1_poeorb(in_poeorb_dir="/home/wafeng/soft/InSAR/ORB/Sentinel_1/aux_poeorb/",mission=None):
    #
    if mission is None:
        searchStr = ""
    else:
        searchStr = mission+"*"
    stafiles = glob.glob(in_poeorb_dir+"/"+searchStr+"*_POEORB_*.EOF")
    timeinfo = np.zeros([len(stafiles),3])
    for index in range(len(stafiles)):
        #
        cfile = os.path.basename(stafiles[index]).split('_')
        times1 = pSAR.ts.dates2jd(cfile[6][1:])
        times2 = pSAR.ts.dates2jd(cfile[7][:-4])
        times0 = pSAR.ts.dates2jd(cfile[5])
        timeinfo[index,:] = [times1,times2,times0]
    #
    return stafiles, timeinfo
###############################################################################    
def s1_resorb(in_resorb_dir="/home/wafeng/soft/InSAR/ORB/Sentinel_1/aux_resorb/",mission=None):
    #
    #
    if mission is None:
        searchStr = ""
    else:
        searchStr = mission+"*"
    #
    stafiles = glob.glob(in_resorb_dir+"/"+searchStr+"*_RESORB_*.EOF")
    # print('FWP:',stafiles)
    timeinfo = np.zeros([len(stafiles),3])
    for index in range(len(stafiles)):
        #
        cfile = os.path.basename(stafiles[index]).split('_')
        times1 = pSAR.ts.dates2jd(cfile[6][1:])
        times2 = pSAR.ts.dates2jd(cfile[7][:-4])
        times0 = pSAR.ts.dates2jd(cfile[5])
        timeinfo[index,:] = [times1,times2,times0]
    #
    return stafiles, timeinfo
#      
###############################################################################  
#        
def rs2_xml2orb(inxml):
    #
    orb = None
    DOMTree    = xml.dom.minidom.parse(inxml)
    collection = DOMTree.documentElement
    tags       = collection.getElementsByTagName("orbitInformation")
    if len(tags) > 0:
       conts = tags[0].getElementsByTagName("orbitDataFile")
       if len(conts) > 0:
          orb = conts[0].childNodes[0].data      
    return orb
###############################################################################    
def gamma_cfginfo2file(info,outcfg):
    #
    azi_post = 50
    rng_post = 50
    dem_ovr  = 2
    npoly    = 3
    fftsize  = 64
    alpha    = 0.85
    cc_win   = 7
    pbase    = 500.
    cc_thres = 0.15
    rlooks = ''
    azlooks = ''
    tempbase = 24
    master = ''
    geo_dir = ''
    dem_path = ''
    dintdir = ''
    #
    if "ISAR_AZI_POS" in info:
        azi_post = int(info['ISAR_AZI_POS'][0])
    if "ISAR_RNG_POS" in info:
        rng_post = int(info['ISAR_RNG_POS'][0])
    if "DEM_OVR" in info:
        #print(info['DEM_OVR'][0])
        dem_ovr = int(float(info['DEM_OVR'][0]))
    if "ISAR_NPOLY" in info:
        npoly = int(info['ISAR_NPOLY'][0])    
    if "ISAR_FILTER_FFTSIZE" in info:
        fftsize = int(info['ISAR_FILTER_FFTSIZE'][0])  
    if "ISAR_FILTER_ALPHA" in info:
        alpha = float(info['ISAR_FILTER_ALPHA'][0])     
    if 'ISAR_CC_WIN' in info:
        cc_win = int(info['ISAR_CC_WIN'][0])
    if 'ISAR_MAX_SPA_BASELINE' in info:
        pbase = float(info['ISAR_MAX_SPA_BASELINE'][0])    
    if 'ISAR_UNW_CC_THRESH' in info:
        cc_thres = float(info['ISAR_UNW_CC_THRESH'][0]) 
    if "ISAR_ZLOOKS" in info:
        azlooks = (info['ISAR_ZLOOKS'][0])  
    if "ISAR_RLOOKS" in info:
        azlooks = (info['ISAR_ZLOOKS'][0])  
    if "ISAR_MAX_TEM_BASELINE" in info:
        tempbase = int(info['ISAR_MAX_TEM_BASELINE'][0])  
    if "ISAR_MASTER" in info:
        master = (info['ISAR_MASTER'][0])    
    if "GEOCODING_DIR" in info:
        geo_dir = (info['GEOCODING_DIR'][0])   
    if "DEM_PATH" in info:
        dem_path = (info['DEM_PATH'][0])   
    if "DINT_DIR" in info:
        dintdir = (info['DINT_DIR'][0])     
    #
    ginsar_cfg(outcfg,fftsize=fftsize,cc_win=cc_win,alpha=alpha,\
               dem_path=dem_path,master=master,\
               s1_dir='NULL',omp_num=6,cc_thres=cc_thres,\
               tempbase=tempbase,pbase=pbase,\
               min_lat=-90,max_lat=90,rlooks=rlooks,azlooks=azlooks,\
               npoly=npoly,njob=4,wfrac=0.15,\
               usr_inf_list='',start_step="baseline",\
               rng_post=rng_post,azi_post=azi_post,\
               end_step="gamma2out",event_date='',\
               geo_dir=geo_dir,dintdir=dintdir,dem_ovr=dem_ovr)
    return True
#
###############################################################################
def ginsar_cfg(outcfg,fftsize=64,cc_win=7,alpha=0.85,dem_path='NULL',master="",\
            s1_dir='NULL',omp_num=6,cc_thres=0.125,tempbase=24,pbase=500.,\
            min_lat=-90,max_lat=90,rlooks='',azlooks='',\
            npoly=4,njob=4,wfrac=0.05,\
            usr_inf_list='',start_step="baseline",\
            rng_post=50,azi_post=50,roi='0.,0.,0.,0.',\
            end_step="gamma2out",event_date='',seasonal=0,\
            geo_dir='',dintdir='',dem_ovr=2,refine=0):
    #
    intime = datetime.datetime.now()
    fid = open(outcfg,'w')
    #
    fid.write('*********************************************\n')
    fid.write('*Processing Control Parameters for InSAR processing by using gInSAR \n')
    fid.write('*Starting at %s\n' % intime)
    fid.write('*Copyrights @ NRCan\n')
    fid.write('*********************************************\n')
    fid.write('%-25s %d\n' % ("OMP_NUM_THREADS",      omp_num))
    fid.write('%-25s %d\n' % ("PARALLEL_WORKERS",     njob))
    fid.write('%-25s %s\n' % ("RAW_S1A_DATA_DIR",     s1_dir))
    fid.write('%-25s %s\n' % ("DEM_PATH",             dem_path))
    fid.write('%-25s %f\n' % ("DEM_OVR",              dem_ovr))
    fid.write('%-25s %s\n' % ("GEOCODING_DIR",        geo_dir))
    fid.write('%-25s %s\n' % ("DINT_DIR",             dintdir))
    fid.write('%-25s %s\n' % ("ISAR_START_STEP",      start_step))
    fid.write('%-25s %s\n' % ('ISAR_END_STEP',        end_step))
    fid.write('%-25s %s\n' % ('ISAR_INF_USER_LIST',   usr_inf_list))
    fid.write('%-25s %d\n' % ('ISAR_INF_SEASONAL',    seasonal))
    fid.write('%-25s %s\n' % ('ISAR_MASTER',          master))
    fid.write('%-25s %s\n' % ('ISAR_EVENT_DATE',      event_date))
    fid.write('%-25s %s\n' % ('ISAR_SUBROI',          roi))
    fid.write('%-25s %s\n' % ('ISAR_RLOOKS',          str(rlooks)))
    fid.write('%-25s %s\n' % ('ISAR_ZLOOKS',          str(azlooks)))
    fid.write('%-25s %d\n' % ('ISAR_RNG_POS',         rng_post))
    fid.write('%-25s %d\n' % ('ISAR_AZI_POS',         azi_post))
    fid.write('%-25s %d\n' % ('ISAR_DEMASSIST',       0))
    fid.write('%-25s %d\n' % ('ISAR_REFINE_RSLC',    refine))
    fid.write('%-25s %d\n' % ('ISAR_NPOLY',           npoly))
    fid.write('%-25s %f\n' % ('ISAR_MAX_SPA_BASELINE',pbase))
    fid.write('%-25s %d\n' % ('ISAR_MAX_TEM_BASELINE',tempbase))
    fid.write('%-25s %d\n' % ('ISAR_CC_WIN',          cc_win))
    fid.write('%-25s %f\n' % ('ISAR_UNW_CC_THRESH',   cc_thres))
    fid.write('%-25s %f\n' % ('ISAR_FILTER_ALPHA',    alpha))
    fid.write('%-25s %d\n' % ('ISAR_FILTER_FFTSIZE',  fftsize))
    fid.write('%-25s %s\n' % ('ISAR_FILTER_WFRAC',    str(wfrac)))
    fid.write('*********************************************\n')
    fid.write('*End!\n')
    fid.write('*********************************************\n')
    fid.close()
    #
    if os.path.exists(outcfg):
       return True
    else: 
       return False
#       
###############################################################################
def s1a_cfg(outcfg,fftsize=64,cc_win=7,alpha=0.95,dem_path='NULL',\
            s1_dir='NULL',omp_num=12,cc_thres=0.125,tempbase=24,pbase=400.,\
            min_lat=-90,max_lat=90,rlooks=20,azlooks=4):  
    #
    intime = datetime.datetime.now()
    fid = open(outcfg,'w')
    #
    fid.write('*********************************************\n')
    fid.write('Processing Control Parameters for Sentinel-1 TOPS data processing\n')
    fid.write('Starting at %s\n' % intime)
    fid.write('Copyrights @ NRCan\n')
    fid.write('*********************************************\n')
    fid.write('%-25s %d\n' % ("OMP_NUM_THREADS",omp_num))
    fid.write('%-25s %s\n' % ("RAW_S1A_DATA_DIR",s1_dir))
    fid.write('%-25s %s\n' % ("DEM_PATH",dem_path))       
    fid.write('%-25s %f\n' % ('S1A_MIN_LAT',min_lat))
    fid.write('%-25s %f\n' % ('S1A_MAX_LAT',max_lat))
    fid.write('%-25s %d\n' % ('S1A_RLOOKS', rlooks))
    fid.write('%-25s %d\n' % ('S1A_ZLOOKS', azlooks))
    fid.write('%-25s %f\n' % ('S1A_MAX_SPA_BASELINE',pbase))
    fid.write('%-25s %d\n' % ('S1A_MAX_TEM_BASELINE',tempbase))
    fid.write('%-25s %d\n' % ('S1A_CC_WIN',cc_win))
    fid.write('%-25s %f\n' % ('S1A_UNW_CC_THRESH',cc_thres))
    fid.write('%-25s %f\n' % ('S1A_FILTER_ALPHA',alpha))
    fid.write('%-25s %d\n' % ('S1A_FILTER_FFTSIZE',fftsize))
    fid.write('*********************************************\n')
    fid.write('End!\n')
    fid.write('*********************************************\n')
    fid.close()
    return True    
###############################################################################
def gamma_s1a_roionsubwath(slc_par,slc_top_par,slavepar=None,minlat=-90,maxlat=90): 
    #
    startid = 1
    endid   = None
    #
    numburst  = gamma_s1a2burst(slc_top_par)
    #
    if slavepar is not None:
       # val,flag = gamma_s1a_commonbursts(slc_par,slavepar)
       ids1 = gamma_s1a_chkburst(slc_par,None,minlat=minlat,maxlat=maxlat)
       ids2 = gamma_s1a_chkburst(slavepar,None,minlat=minlat,maxlat=maxlat)
       numburst2 = gamma_s1apar2burst(slavepar)
       val  = [ids1[0],ids2[0],numburst,numburst2,0.001]
       # print(minlat,maxlat,ids1)
       #
       if (ids1[0] < numburst2 and ids2[0] < numburst):
          flag = 1
       else:
          flag = 0
        
       if flag == 1:
          startid = val[0]
          endid   = val[2]
    outid = gamma_s1a_chkburst(slc_par,numburst,minlat=minlat,maxlat=maxlat,\
                                          startid=startid,endid=endid)
    return outid
###############################################################################  
def gamma_s1apar2burst(slc_par):
    #
    slc_par = os.path.abspath(slc_par)
    slc_par_dirname = os.path.dirname(slc_par)
    slc_par_rootname= os.path.basename(slc_par).split('.')[0]
    slc_toppar      = os.path.join(slc_par_dirname,slc_par_rootname+'.tops_par')
    numburst        = gamma_s1a2burst(slc_toppar)
    return numburst
#    
###############################################################################
def gamma_updateRS2orb(slc_par,orbfile,nstat=64):
    #
    logfile = os.path.dirname(slc_par)+'/'+os.path.basename(slc_par).split('.')[0]+'.log'
    cmd = 'RSAT2_vec %s %s %d > %s' % (slc_par,orbfile,nstat,logfile)
    print(' pGAMMA: %s' % cmd)
    os.system(cmd)
    return cmd
#
def gamma_opodvector(slc_par, statevector, nstat=64):
    #
    if (statevector not in ["NULL"] and os.path.exists(statevector)):
       gamma_sTr = ("S1_OPOD_vec %s %s %d" % (slc_par, statevector, nstat))
       subinfo = subprocess.Popen(gamma_sTr, stdout=subprocess.PIPE,\
                                   stderr=subprocess.PIPE,shell=True)
       subinfo.wait()
       outlog = subinfo.communicate()[0]
       return outlog
    else:
       return None
#
###############################################################################
def gamma_s1times2statevector(startime,stoptime,model="RESORB",orbdir=None,mission='S1A'):
    #
    if orbdir is None:
       topdir = os.environ['S1_ORB']
       orbdir = topdir+'/aux_'+model.lower()+'/'
    #   
    st = pSAR.ts.timestr2jd(startime)
    et = pSAR.ts.timestr2jd(stoptime)
    #
    avt = (st+et)/2.
    #
    if model.upper() == "RESORB":
       if orbdir is not None:
          # print(model,orbdir,mission)
          orbs,times = s1_resorb(in_resorb_dir=orbdir,mission=mission)
       else:
          orbs,times = s1_resorb(mission=mission)
    else:
       if orbdir is not None:
           orbs, times = s1_poeorb(in_poeorb_dir=orbdir,mission=mission)
       else:
           orbs, times = s1_poeorb(mission=mission)
    #
    #
    flag1 = times[:,0] <= st
    flag2 = times[:,1] >= et
    flag3 = np.logical_and(times[:,0] <= avt,times[:,1] >= avt)
    flag0 = np.logical_and(flag1,flag2)
    #
    orbs   = np.array(orbs)
    ctimes = times[flag0,2]
    outorb = orbs[flag0]
    if len(outorb)==0:
        outorb = orbs[flag3]
        ctimes = times[flag3,2]
    #
    #
    if len(outorb)==0:
        return "NULL"
    if isinstance(outorb, np.str) and os.path.exists(outorb):
        return outorb
    else:
        #
        outorb = outorb[ctimes==np.max(ctimes)][0]
        return outorb
#
###############################################################################
#
def gamma_s1par2statevector(s1_slc_par,model="RESORB",orbdir=None):
    #
    avt,st,et = gamma_s1par2time(s1_slc_par)
    #
    # updated by Wanpeng FEng, @NRCan, 2017-02-07
    # be careful about the mission of the data.
    # 
    mission   = gamma_s1par2mission(s1_slc_par)
    #
    if model.upper() == "RESORB":
       if orbdir is not None:
          orbs,times = s1_resorb(in_resorb_dir=orbdir,mission=mission)
       else:
          orbs,times = s1_resorb(mission=mission)
    else:
       if orbdir is not None:
           orbs, times = s1_poeorb(in_poeorb_dir=orbdir,mission=mission)
       else:
           orbs, times = s1_poeorb(mission=mission)
    #
    flag1 = times[:,0] <= st
    flag2 = times[:,1] >= et
    flag3 = np.logical_and(times[:,0] <= avt,times[:,1] >= avt)
    flag0 = np.logical_and(flag1,flag2)
    #
    # print(flag0.shape)
    orbs   = np.array(orbs)
    ctimes = times[flag0,2]
    outorb = orbs[flag0]
    if len(outorb)==0:
        outorb = orbs[flag3]
        ctimes = times[flag3,2]
    #
    #
    if len(outorb)==0:
        return "NULL"
    if isinstance(outorb, np.str) and os.path.exists(outorb):
        return outorb
    else:
        #
        outorb = outorb[ctimes==np.max(ctimes)][0]
        return outorb
###############################################################################
def gamma_slcs2tab(indir,out_tab,ext='.slc',fixdate=None):
    '''
    To create a tab file with SLC and SLC_par files,
    which will be widely used in the InSAR processing with the GAMMA soft
    
    '''
    if fixdate is None:
       slcs = glob.glob(indir+'/*'+ext)
       if len(slcs) == 0:
          slcs = glob.glob(indir+'/*/*'+ext)
    else:
       #
       slcs = glob.glob(indir+'/%s*' % fixdate + ext)
       if len(slcs) == 0:
          slcs = indir+'/%s' % fixdate+'/%s*' % fixdate + ext
          slcs = glob.glob(slcs)
    #
    # a bug was fixed by FWP, @CCRS/NRCan, 2017-06-08
    # a new list will be in an alphabetical order
    #
    slcs.sort()
    with open(out_tab,'w') as fid:
        for ind in range(len(slcs)):
            fid.write("%s %s\n" % (slcs[ind],slcs[ind]+'.par'))
    if os.path.exists(out_tab):
        return len(slcs),True
    else:
        return len(slcs),False
#
###############################################################################
#    
def gamma_s1getstatevect(indate,model='RESORB',\
                         s1orb="/home/wafeng/soft/InSAR/ORB/Sentinel_1/"):  
    #
    # Two types of State vectors: POEORB and RESORB
    # 
    predate  = pSAR.ts.newdate(indate,-1)
    postdate = pSAR.ts.newdate(indate,1)
    # print(predate,postdate)
    statfile = glob.glob(s1orb+'/*'+model+'_OPOD_*_V'+str(int(predate))\
                          +"T*_"+str(int(postdate))+"*.EOF")
    if len(statfile)>0:
       return statfile[0]
    else:
       return "NULL"
###############################################################################
def gamma_baseline_est(master_par,slave_par,logid="PIPE"):
    '''
    To estimate perpendicular baseline with base_oribt from GAMMA .pars
    Updated by Wanpeng Feng, @CCRS/NRCan, 2017-05-29
    
    '''
    gamma_sTr=(" base_orbit %s %s -" % (master_par,slave_par))
    print(" pGAMMA: %s" % gamma_sTr)
    #
    subinfo = subprocess.Popen(gamma_sTr, stdout=subprocess.PIPE,\
                                   stderr=subprocess.PIPE,shell=True)
    subinfo.wait()
    outlog = subinfo.communicate()[0]
    outlog = pDATA.bytestoutf8(outlog)
    
    # print(outlog)
    newlog = outlog.split('\n')
    perpb  = 0.
    parab  = 0.
    for cnewlog in newlog:
        #
        if "baseline perpendicular component (m):" in cnewlog:
            tmp   = cnewlog.split(":")[1]
            perpb = float(tmp.split()[0])
        if "baseline parallel component (m):" in cnewlog:
            tmp   = cnewlog.split(":")[1]
            parab = float(tmp.split()[0])
            #
    return perpb,parab
    
###############################################################################
def gamma_par2statenum(slc_par):
    #
    numstate = 0
    fid = open(slc_par,'r')
    for cline in fid:
        cline = cline.split('\n')[0]
        if "number_of_state_vectors" in cline:
            tmp = cline.split()
            numstate=int(float(tmp[1]))
    return numstate
#
###############################################################################
def gamma_s1SLC2toppar(slc):
    slc_dirname = os.path.dirname(slc)
    slc_rootname= os.path.basename(slc).split('.')[0]
    slc_toppar      = os.path.join(slc_dirname,slc_rootname+'.tops_par')
    return slc_toppar
#    
def gamma_s1a_chkburst(slc_par,numburst,roipoly=None,minlat=-90,maxlat=90,\
                       startid=1,endid=None,ovr=0.02):
    #
    slc_par = os.path.abspath(slc_par)
    #
    if numburst is None:
       #
       slc_par_dirname = os.path.dirname(slc_par)
       slc_par_rootname= os.path.basename(slc_par).split('.')[0]
       slc_toppar      = os.path.join(slc_par_dirname,slc_par_rootname+'.tops_par')
       numburst        = gamma_s1a2burst(slc_toppar)
    # geopoly = gamma_s1a_geopoly(slc_par,numburst,ovr=ovr)
    geopoly,numburst = gamma_s1a_burstcors(slc_par,None)
    npoly   = geopoly.shape[0]
    outid   = []
    c       = 0
    if endid is None:
       endid = numburst
    #
    if roipoly is not None:
       refpoly = Polygon(roipoly)
    #
    for ni in range(npoly):
        geopolygon = np.reshape(geopoly[ni,:],(5,2))
        # print(geopolygon)
        cminlat    = np.min(geopolygon[:,0])
        cmaxlat    = np.max(geopolygon[:,0])
        # 
        if (ni+1 >= startid and ni+1 <= endid and roipoly is None):
           #
           if (cminlat >= maxlat or cmaxlat <= minlat):
              c += 1
           else:
              outcid = ni + 1
              outid.append(outcid)
        else:
           lons = geopolygon[:,1]
           lats = geopolygon[:,0]
           #
           cpoly = Polygon(np.vstack((lons,lats)).T)
           flag  = refpoly.intersects(cpoly)
           if flag:
               outcid = ni + 1
               outid.append(outcid)
           else:
               c += 1
              
     #
    outid = np.array(outid)
    return outid
###############################################################################
def gamma_s1a2geopoly(s1a_iw_par,ovr=0.02):
    #
    s1a_slc_par = os.path.abspath(s1a_iw_par)
    s1a_slc_dir = os.path.dirname(s1a_slc_par)
    s1a_slc_rootname = os.path.basename(s1a_slc_par).split('.')[0]
    s1a_slc_toppar   = os.path.join(s1a_slc_dir,s1a_slc_rootname+'.tops_par')
    #
    if os.path.exists(s1a_slc_toppar):
        # numburst         = gamma_s1a2burst(s1a_slc_toppar)
        # geopoly          = gamma_s1a_geopoly(s1a_iw_par,numburst,ovr=0.02)
        geopoly,numburst = gamma_s1a_burstcors(s1a_iw_par,None)
    else:
        print(" ERROR: TOPS: %s" % s1a_slc_toppar)
        geopoly = None
        numburst = None
    #
    return geopoly,numburst
    #
###############################################################################
#    
def gamma_2pointsDIST(x1,y1,x2,y2):
    #
    return np.sqrt((x1-x2)**2 + (y1-y2)**2)
###############################################################################    
def gamma_se1_checktopsdiff(s1a_par1,s1a_par2):
    #
    geopoly1,burst1 = gamma_s1a_burstcors(s1a_par1,None)
    geopoly2,burst2 = gamma_s1a_burstcors(s1a_par2,None)
    #
    outdata = []
    if (geopoly1 is not None and geopoly2 is not None):
        #
        nburst1 = geopoly1.shape[0]
        nburst2 = geopoly2.shape[0]
        for ni in range(nburst1):
            poly1 = geopoly1[ni,:]
            poly1 = np.reshape(poly1,(5,2))
            if ni <= nburst2-1:
              poly2 = geopoly2[ni,:]
              poly2 = np.reshape(poly2,(5,2)) 
              lat1_m,lon1_m = np.mean(poly1[:,0]),np.mean(poly1[:,1])
              lat2_m,lon2_m = np.mean(poly2[:,0]),np.mean(poly2[:,1])
              dist          = gamma_2pointsDIST(lon1_m,lat1_m,lon2_m,lat2_m)
            else:
              dist = -1.
            outdata.append([ni,dist])
    #
    return outdata
###############################################################################
def gamma_slcinfo(slc_par,keyword=None):
    #
    fid = open(slc_par,'r')
    #
    width           = 0
    file_length     = 0
    range_spacing   = 0.
    azimuth_spacing = 0.
    wavelength      = 0.
    #
    slcinfo = {}
    for cline in fid:
        cline = cline.split('\n')[0]
        if "range_samples" in cline:
            tmp = cline.split()
            width=int(float(tmp[1]))
        if "azimuth_lines" in cline:
            tmp = cline.split()
            file_length=int(float(tmp[1]))
        if "range_pixel_spacing" in cline:
            tmp = cline.split()
            range_spacing=float(tmp[1])
        if "azimuth_pixel_spacing" in cline:
            tmp = cline.split()
            azimuth_spacing=float(tmp[1])
        if "radar_frequency" in cline:
            tmp = cline.split()
            wavelength=3.0e10/float(tmp[1])
        #    
        if keyword is not None:
            if keyword in cline:
                tmp              = cline.split()
                slcinfo[keyword] = float(tmp[1])
    #
    slcinfo = {"WIDTH":width,"FILE_LENGTH":file_length,"WAVELENGTH":wavelength,\
               "RANGE_SPACING":range_spacing,"AZIMUTH_SPACING":azimuth_spacing}
    return slcinfo
#
#
###############################################################################
def gamma_s1iw_cb_par_dir(slc_par1,slc_2_dir,mindist=0.02):
    #
    swath = os.path.basename(slc_par1).split('_')[1]
    slc_par2 = glob.glob(slc_2_dir+'/*%s*.slc.par' % swath)
    flag = [0,0]
    pairinfo = [2,2]
    if len(slc_par2) > 0:
       pairinfo,flag = gamma_s1iw_commonbursts(slc_par1,slc_par2[0],mindist=mindist)
    return pairinfo,flag
#
def gamma_s1iw_commonbursts(slc_par1,slc_par2,mindist=0.02):
    '''
    return the perfect common bursts between two TOPS SLC folders
    #
    '''
    # import matplotlib.pyplot as plt
    #
    geopoly1,burst1 = gamma_s1a_burstcors(slc_par1,None)
    geopoly2,burst2 = gamma_s1a_burstcors(slc_par2,None)
    #
    #
    dist_Mj = np.zeros(burst2)
    pairinfo = np.zeros([burst1,3])
    #
    for ni in range(burst1):
        cpoly_i = np.reshape(geopoly1[ni,:],[5,2])
        mean_loni = np.mean(cpoly_i[:,1])
        mean_lati = np.mean(cpoly_i[:,0])
        #
        # plt.plot(cpoly_i[:,0],cpoly_i[:,1],'-r',linewidth=5)
        #
        #
        for nj in range(burst2):
            cpoly_j = np.reshape(geopoly2[nj,:],[5,2])
            mean_lonj = np.mean(cpoly_j[:,1])
            mean_latj = np.mean(cpoly_j[:,0])
            #
            # plt.plot(cpoly_j[:,0],cpoly_j[:,1],'-b',linewidth=1)
            # plt.text(cpoly_j[:,0].mean(),cpoly_j[:,1].mean(),'%d' % nj)
            #
            # Conduct the distances of two bursts
            #
            dist_Mj[nj] = np.sqrt((mean_loni-mean_lonj)**2 + \
                                  (mean_lati-mean_latj)**2)
        #
        # plt.show()
        # make sure there is no Nan in dist_Mj
        # Why Nan could be seen here?
        #
        dist_Mj[np.isnan(dist_Mj)] = 10000.0
        index = np.where(dist_Mj == np.min(dist_Mj))[0]
        #
        # print(dist_Mj)
        # print(dist_Mj[index])
        #
        if len(index) > 1:
            print(" ERROR: please double check input slc_par carefully...")
            return None,[0,0]
        #
        pairinfo[ni,:] = [ni,index,dist_Mj[index]]
    #
    pairinfo = pairinfo[pairinfo[:,2]<mindist,:]
    if pairinfo.shape[0] > 0:
       pairinfo[:,0] = pairinfo[:,0] + 1
       pairinfo[:,1] = pairinfo[:,1] + 1
    #
    flag = [0,0]
    #
    if pairinfo.shape[0] < burst1:
        flag[0] = 1
    if pairinfo.shape[0] < burst2:
        flag[1] = 1
    if pairinfo.shape[0] == 0:
        flag = [2,2]
    #
    return pairinfo,flag
            
def gamma_se1_commonburst(slc_par1,slc_par2):
   #
   val,flag = gamma_s1a_commonbursts(slc_par1,slc_par2)
   #
   if val is None:
      val = [0,0,0,0,1]
   #
   if (val[4] > 0.02 or val is None):
      outdata = [-1,-1,-1,-1]
   else:
      #
      mstart1=val[0]
      sstart2=val[1]
      mpburts=val[2]-val[0]+1
      spburts=val[3]-val[1]+1
      #
      if mpburts > spburts:
         numburst = spburts
      else:
         numburst = mpburts
      #
      outdata = [mstart1,mstart1+numburst-1,sstart2,sstart2+numburst-1]
   return outdata
#
###############################################################################  
#
def gamma_s1burst_slave(s1a_par1,s1a_par2,master_burst,ovr=0.02,minimumdist=0.05):
    #
    val, flag = gamma_s1a_commonbursts(s1a_par1,s1a_par2,ovr=ovr,\
                                       minimumdist=minimumdist)
    return val
#
###############################################################################
#
def gamma_slcpar2date(par):
    #
    # return a date in YYYYMMDD format
    #
    info = gamma_slcpar(par)
    yr = info['date'][0]
    mn = info['date'][1]
    dd = info['date'][2]
    if len(mn) < 2:
        mn = '0%s' % mn
    if len(dd) < 2:
        dd = '0%s' % dd
    return yr+mn+dd
#
###############################################################################
#
def gamma_slcpar2CT(slcpar):
    #
    # return the central time
    # in a format like 20180215T15:01:02.5
    # fmt is %Y%m%dT%H:%M:%S.%f
    #
    fmt = '%Y%m%dT%H:%M:%S.%f'
    st,et = gamma_slcpar2timeexact(slcpar)
    st_jd = pSAR.ts.timestr2jd(st,fmt=fmt)
    et_jd = pSAR.ts.timestr2jd(et,fmt=fmt)
    ct_jd = (st_jd+et_jd) / 2.
    #
    return pSAR.ts.jd2timestr(ct_jd,fmt=fmt)
#
###############################################################################
#
def gamma_slcpar2timeexact(par):
    #
    # cdate = gamma_slcpar2date(par)
    info = gamma_slcpar(par)
    start_time = float(info['start_time'][0])
    end_time = float(info['end_time'][0])
    #
    #
    start_exac = datetime.datetime(int(info['date'][0]),\
                                   int(info['date'][1]),\
                                   int(info['date'][2])) + \
                 datetime.timedelta(seconds=start_time)
    end_exac = datetime.datetime(int(info['date'][0]),\
                                   int(info['date'][1]),\
                                   int(info['date'][2])) + \
                 datetime.timedelta(seconds=end_time)
    #
    return pSAR.ts.dt2timestr(start_exac),pSAR.ts.dt2timestr(end_exac)
#   
def gamma_s1a_commonbursts(s1a_par1,s1a_par2,ovr=0.02,minimumdist=0.05):
    #
    #
    geopoly1,burst1 = gamma_s1a_burstcors(s1a_par1,None)
    geopoly2,burst2 = gamma_s1a_burstcors(s1a_par2,None)
    # print(burst1,burst2)
    #
    flag = 0
    outval = None
    if (geopoly1 is not None and geopoly2 is not None):
        #
        nburst1 = geopoly1.shape[0]
        nburst2 = geopoly2.shape[0]
        for ni in range(nburst1):
            #
            outdata = []
            for nj in range(nburst2):
                poly1 = geopoly1[ni,:]
                poly2 = geopoly2[nj,:]
                poly1 = np.reshape(poly1,(5,2))
                poly2 = np.reshape(poly2,(5,2)) 
                lat1_m,lon1_m = np.mean(poly1[:,0]),np.mean(poly1[:,1])
                lat2_m,lon2_m = np.mean(poly2[:,0]),np.mean(poly2[:,1])
                dist          = gamma_2pointsDIST(lon1_m,lat1_m,lon2_m,lat2_m)
                outdata.append(dist)
            #
            outdata = np.array(outdata)
            ind     = outdata.argmin()
            # print("M: %f S: %f -> %f" % (ni+1,ind+1,outdata[ind]))
            if (flag == 0 and outdata[ind] < minimumdist):
                flag = 1
                outval = [ni+1,ind+1,burst1,burst2,outdata[ind]]
            #
                
    return outval, flag
###############################################################################
# convert ALOS2 to GAMMA format
#
def gamma_readalos2(leaderfile,imgfile,outslc,outslc_par):
    #
    gamma_command = (" par_EORC_PALSAR %s %s %s %s" % (leaderfile,outslc_par,\
                                                       imgfile,outslc))
    print(gamma_command)
    subinfo = subprocess.Popen(gamma_command, stdout=subprocess.PIPE,\
                                   stderr=subprocess.PIPE,shell=True)
    subinfo.wait()
    #
    # output = subinfo.communicate()[0]
    # outstr = output.split('\n')
    #
    if os.path.exists(outslc):
       return True
    else:
       return False
###############################################################################
#
def gamma_slc_cat(slc1,slc2,slcpar1,slcpar2,offpar,outslc,outslcpar,dopflag=1,\
                  iflg=0,phflg=0):
    #
    gamma_command = ("SLC_cat %s %s %s %s %s %s %s %d %d %d" % \
                     (slc1,slc2,slcpar1,slcpar2,offpar,outslc,outslcpar,\
                     dopflag,iflg,phflg))
    #
    flag,val,errors = gamma_run(gamma_command)
    #
    return flag,val,errors
#                            
###############################################################################
def gamma_s1_burstNUMs_refpolys(slcpar,refpolygon,area_thresh=0.10):
    #flag_matrix
    cpoly, bursNum = gamma_s1a_burstcors(slcpar,None)
    #
    burst_ids = np.zeros(cpoly.shape[0])
    #
    for i in range(cpoly.shape[0]):
        Poly1 = np.reshape(cpoly[i,:],[5,2])
        flag_matrix = np.zeros(refpolygon.shape[0]) 
        #
        for j in range(refpolygon.shape[0]):
            #
            Poly2 = np.reshape(refpolygon[j,:],[5,2])
            carea = pDATA.overlapOftwopolygons(Poly1,Poly2)
            #
            # debugging @SYSU, FWP, 2019/07/10
            #if '20160707' in slcpar and i == 0:
            #   print('%d - %d: %f' % (i,j,carea))
            #
            if carea > area_thresh:
                flag_matrix[j] = 1
            #
        #
        index = np.where(flag_matrix==1)[0]
        # print(index)
        if index.shape[0] > 1:
            print(' ERROR: a repeatation of bursts may be found in refpolygons')
            return None
        #
        burst_ids[i] = index[0]
        #
    return burst_ids
#
def gamma_s1_dir2geopoly(inDir,subswath=1,area_thresh=0.10):
    #
    #
    counter = 0
    slcpars = glob.glob('%s/*/*iw%d_*.slc.par' % (inDir,subswath))
    #
    if len(slcpars)<1:
        return None,None
    # all data in the same track
    info = gamma_slcpar(slcpars[0])
    #
    # in default, -12 for ascedning and -165 for descending
    # We must review this in future. by Wanpeng Feng, @SYSU,Guangzhou, 2019/07/10
    #
    heading_deg = float(info['heading'][0])
    #
    for cslcpar in slcpars:
        #
        tops_par = cslcpar.split('.slc')[0]+'.tops_par'
        #
        cpoly,numOfbursts = gamma_s1a_burstcors(cslcpar,tops_par)
        #
        nanFlag = np.isnan(cpoly)
        if nanFlag.sum() > 0:
            print( "ERROR: something wring with %s " % cslcpar)
        #
        #
        if counter == 0:
            outpoly = cpoly
        else:
            flag_matrix = np.zeros(cpoly.shape[0]) + 1
            #
            for ni in range(cpoly.shape[0]):
                Poly1 = np.reshape(cpoly[ni,:],[5,2])
                for nj in range(outpoly.shape[0]):
                    #
                    Poly2 = np.reshape(outpoly[nj,:],[5,2])
                    carea = pDATA.overlapOftwopolygons(Poly1,Poly2)
                    # print(ni,nj,carea)
                    if carea > area_thresh:
                        flag_matrix[ni] = 0
                    #
            # print(flag_matrix.shape)
            #
            if np.sum(flag_matrix) > 0:
               # print(cslcpar,flag_matrix.sum())
               outpoly = np.vstack((outpoly,cpoly[flag_matrix==1,:]))
            #
        #
        counter += 1
    #
    cenPoints = np.zeros([outpoly.shape[0],2])
    #
    for index in range(outpoly.shape[0]):
        #
        tmppoly = np.reshape(outpoly[index,:],[5,2])
        cenPoints[index,:] = [tmppoly[:,0].mean(),tmppoly[:,1].mean()]
    #
    # we must sort the polygons based on lats...
    # ascending in default
    index = np.argsort(cenPoints[:,0])
    #
    if heading_deg < -100 and heading_deg > -250:
        # descending
        index = index[::-1]
    #
    cenPoints = cenPoints[index,:]
    outpoly = outpoly[index,:]
    #
    return outpoly,cenPoints
#
def gamma_s1a_burstcors(slc_par,toppar):
    '''
    This utilizes the internal function of GAMMA to calculate coordinates 
    of each burst, SLC_burst_corners
    updated by Wanpeng Feng, @Ottawa, 2016-11-02
    '''
    #
    slc_par = os.path.abspath(slc_par)
    slc_par_dirname = os.path.dirname(slc_par)
    slc_par_rootname= os.path.basename(slc_par).split('.')[0]
    #
    if toppar is None:
        toppar = glob.glob(slc_par_dirname+'/'+slc_par_rootname+'*.tops_par')
        toppar = toppar[0]
    #
    
    gamma_command = ("SLC_burst_corners %s %s" % (slc_par,toppar)) 
    subinfo = subprocess.Popen(gamma_command, stdout=subprocess.PIPE,\
                                   stderr=subprocess.PIPE,shell=True)
    subinfo.wait()
    output = subinfo.communicate()[0]
    output = pDATA.bytestoutf8(output)
    #
    outstr = output.split('\n')
    outdata = []
    #
    for cline in outstr:
        if "Burst:" in cline:
            # print(cline)
            cline = cline.split('\n')[0]
            tmp   = cline.split()
            p1    = [np.float32(tmp[3]),np.float32(tmp[2])]
            p2    = [np.float32(tmp[5]),np.float32(tmp[4])]
            p3    = [np.float32(tmp[7]),np.float32(tmp[6])]
            p4    = [np.float32(tmp[9]),np.float32(tmp[8])]
            outdata.append([p1[1],p1[0],p2[1],p2[0],p3[1],p3[0],p4[1],\
                            p4[0],p1[1],p1[0]])
            
    outdata = np.array(outdata)
    return outdata,outdata.shape[0]
#    
###############################################################################    
#
def gamma_slcs2jointpolygon(slc_pars):
    #
    n       = 0
    polygon = []
    #
    for cpar in slc_pars:
      #
      #
      geopoly = gamma_slc_corners(cpar)
      geopoly = np.vstack((geopoly,geopoly[0,:]))
      #
      n += 1
      polygon.append(Polygon(np.copy(geopoly)))
      #
    if n > 1:
      newgeopoly = cascaded_union(polygon)
      x, y = newgeopoly.exterior.coords.xy
    else:
      newgeopoly = polygon[0]
      # print(newgeopoly)
      x,y = newgeopoly.exterior.coords.xy
    #
    
    return np.vstack((y,x)).T
#
###############################################################################
#
def gamma_slc_corners(slc_par):
    #
    gamma_command = ("SLC_corners %s" % slc_par)
    subinfo = subprocess.Popen(gamma_command, stdout=subprocess.PIPE,\
                                   stderr=subprocess.PIPE,shell=True)
    subinfo.wait()
    output = subinfo.communicate()[0]
    #
    # Updated by adoption in python 3.x 
    # output will be in bytes-like format in python 3.x 
    # bytestoutf8, a function in pDATA can help decode the variable 
    # from bytes to utf8
    #
    output = pDATA.bytestoutf8(output)
    outstr = output.split('\n')
    #
    outdata = []
    #
    for cline in outstr:
        #
        #
        if "latitude" in cline and len(cline.split())==6:
          tmp = cline.split()
          tmp = [x for x in tmp if x]
          lonlat = pSAR.roipac.listslice(tmp,[2,5])
          # print(lonlat)
          outdata.append([np.float32(lonlat[0]),np.float32(lonlat[1])]) 
    #
    outdata = np.array(outdata)
    outdata = outdata[[0,1,3,2],:]
    #
    return outdata
#
###############################################################################
#
def gamma_slcpar2time(slc_par):
    '''
    To return time information 
    
    '''
    info = gamma_slcpar(slc_par)
    dateinfo = info['date']
    #
    if len(dateinfo)>5:
       outime = '%s%s%sT%s:%s:%s' % (dateinfo[0],dateinfo[1],\
                                     dateinfo[2],dateinfo[3],\
                                     dateinfo[4],dateinfo[5])
    else:
       outime = '%s%s%s' % (dateinfo[0],dateinfo[1],dateinfo[2])
    return outime
#
def gamma_dempar2rscinfo(dem_par):
    #
    cinfo = gamma_dempar(dem_par)
    rscinfo,keys = pSAR.roipac.roipac_info()
    rscinfo['DATA_FORMAT'] = cinfo['data_format'][0]
    rscinfo['X_FIRST'] = cinfo['corner_lon'][0]
    rscinfo['X_STEP'] = cinfo['post_lon'][0]
    rscinfo['WIDTH'] = cinfo['width'][0]
    rscinfo['FILE_LENGTH'] = cinfo['nlines'][0]
    rscinfo['Y_STEP'] = cinfo['post_lat'][0]
    rscinfo['Y_FIRST'] = cinfo['corner_lat'][0]
    rscinfo['X_MAX'] = int(rscinfo['WIDTH']) - 1
    rscinfo['Y_MAX'] = int(rscinfo['FILE_LENGTH']) - 1
    #
    return rscinfo
#
def gamma_dempar(dem_par):
    return gamma_read_diffpar(dem_par)
#
def gamma_mlipar(mli_par):
    return gamma_read_diffpar(mli_par)
#
def gamma_slcpar(slc_par):
    return gamma_read_diffpar(slc_par)
#
def gamma_slcpar2wl(slc_par):
    '''
    Here the speed of light is fixed at
    3.0e+8 m/s
    #
    by Wanpeng Feng,@CCRS/NRCan, 2017-08-25
    
    '''
    freqs = gamma_slcpar2freqs(slc_par)
    return 3.0e+8 / freqs[0]
#
def gamma_slcpar2freqs(slc_par):
    #
    info = gamma_slcpar(slc_par)
    return [float(info['radar_frequency'][0]),\
            float(info['adc_sampling_rate'][0]),\
            float(info['chirp_bandwidth'][0])]

def gamma_slcpar2dim(slc_par):
    info = gamma_slcpar(slc_par)
    return [int(info['range_samples'][0]),int(info['azimuth_lines'][0])]
def gamma_read_offpar(offpar):
    outpar = gamma_read_diffpar(offpar)
    return outpar
#
###############################################################################  
#
def gamma_offset_fitm(offs,cpps,diffpar,coffsets,thres=0.15,npoly=3):
    #
    # updated by Wanpeng Feng, @NRCan, 2017-06-26
    #
    offset_fitm = 'offset_fitm %s %s %s - %s %f %d' % \
                (offs,cpps,diffpar,coffsets,thres,npoly)
    print(offset_fitm)
    os.system(offset_fitm)
    return True
#
###############################################################################
#
def gamma_offset_pwrm(sim_rdc,in_mli,diffpar,offs,cpps,patchsize,\
                      offsets,n_ovr,coffsets,npoly=3,thres=0.15):
    #
    flag = False
    #
    while not flag and patchsize > 4:
      off_pwrm = 'offset_pwrm %s %s %s %s %s %d %d %s %d' % \
               (sim_rdc,in_mli,diffpar,offs,cpps,\
                        patchsize,patchsize,offsets,n_ovr)
      print(off_pwrm)
      os.system(off_pwrm)
      flag = True
      try:
         offs_data = np.loadtxt(offsets)
         if offs_data.shape[0] < 25:
             print(" ERROR: offset_pwrm returned only %d samples @ PSIZE %d" \
                   % (offs_data.shape[0],patchsize))
             flag = False
         else:
             #
             gamma_offset_fitm(offs,cpps,diffpar,coffsets,thres=thres,\
                               npoly=npoly)
      except:
         print(" ERROR: offset_pwrm did not return a valid %s" % (offsets))
         flag = False
      #
      if not flag:
         patchsize = patchsize / 2
    #
    if os.path.exists(offsets):
       return True
    else:
       return False
#
###############################################################################
#
def gamma_read_diffpar(indiffpar):
    #
    outpar = {}
    fid    = open(indiffpar,'r')
    for cline in fid:
        tmp = cline.split()
        #
        if len(tmp) > 1:
           value = tmp[1:]
           outpar[tmp[0][0:-1]] = value 
    #
    fid.close()
    #
    return outpar 
###############################################################################
def gamma_adf(intfile,filt_int,filt_cc,width,win=[32],alp=[0.6]):
    #
    win = np.array(win)
    alp = np.array(alp)
    nwin = win.shape[0]
    for ind in range(nwin):
        cwin = win[ind]
        calp = alp[ind]
        if ind == 0:
            inint = intfile
        else:
            cind = ind - 1
            inint = ("%s.%d" % (filt_int,cind))
        #
        outint= ('%s.%d' % (filt_int,ind))
        # outcc = ('%s.%d' % (filt_cc,ind))
        if ind == nwin - 1:
            outint = filt_int
        #
        #
        gamma_command = ('adf %s %s %s %d %f %d' % (inint,outint,filt_cc,\
                                                    width,calp,cwin))
        #
        flag,info,error = gamma_run(gamma_command)
        #
        if flag >= 1:
            # print(" pGAMMA: ERROR!")
            print(" pGAMMA: ERROR!!! Check pgamma_error.inf for more details")
            sys.exit(-1)
            
    #
    gamma_command = ('rasmph %s %d' % (filt_int,width))
    flag,info,error = gamma_run(gamma_command)
    return flag
#
###############################################################################   
#
def gamma_mli2pic(inmli,outpic):
    #
    info = gamma_slcpar(inmli+'.par')
    raspwr = 'raspwr %s %s 1 0 1 1 1. 0.35 1 %s' % \
               (inmli,info['range_samples'][0],outpic)
    flag,prnts,errs = gamma_run(raspwr)
    return flag,prnts,errs
#
###############################################################################
#     
def gamma_cpx2bands(incpx,outname,width,outtype):
    #
    # outtype: 
    #       0: real part
    #       1: imaginary part
    #       2: intensity (re*re + im*im)
    #       3: magnitude (sqrt(re*re + im*im))
    #       4: phase (atan2(im, re))

    gamma_command=('cpx_to_real %s %s %s %s' % \
                   (incpx, outname, width, outtype))
    #print(' pGMT: %s' % gamma_command)
    flag,prints,errs = gamma_run(gamma_command)
    return flag

#
###############################################################################
#
def gdal_interpphs(xyz,newphs,otype='ENVI',method='nearest',\
                   xbd=[1,100],ybd=[1,100],outsize=[1,1]):
    #
    gdal_command = ('gdal_grid -of %s -txe %f %f -tye %f %f '+\
                    '-outsize %f %f -l %s -a %s %s %s' %
                    (otype,xbd[0],xbd[1],ybd[0],ybd[1],outsize[0],outsize[1],\
                     xyz,method,xyz,newphs))
    flag,prints,erro = gamma_run(gdal_command)
    return flag
#
###############################################################################
#       
def gamma_maskphs(inphs,cc,thresh,width,length):
    #
    print(' pGAMMA: reading %s in %d %d' % (inphs,width,length))
    phs = gamma_read_data(inphs,width,length)
    ccv = gamma_read_data(cc,   width,length)
    phs[np.isnan(phs)] = 0.
    phs[ccv<thresh] = 0.
    #
    xind,yind = np.where(phs!=0)
    data      = phs[phs!=0]
    #
    x  = np.reshape(xind,(xind.shape[0],1))
    y  = np.reshape(yind,(yind.shape[0],1))
    D  = np.reshape(data,(data.shape[0],1))
    xyz   = np.hstack((y,x,D))
    outxyz = inphs+'.xyz'
    #
    # np.savetxt(outxyz,xyz,fmt='%f %f %f',newline='\n')
    pSAR.roipac.roi_write(xyz,outxyz)
    #
    return outxyz
#
###############################################################################
#
# lon and lat to samples and lines
def gamma_coord_to_sarpix(slc,lon,lat,evl,slcpar=None,diffpar=''):
    #
    if slcpar is None:
       slcpar=('%s.par' % slc)
    #
    if os.path.exists(slcpar) is False:
        print('%s canot be found. Please check it first' % slcpar)
        sys.exit(-1)
    #
    gamma_command = ('coord_to_sarpix %s - - %f %f %f %s' % \
                     (slcpar,lat,lon,evl,diffpar))
    #
    # A pipe was newly added to make sure the prints on the screen can be read
    # into prints
    #
    flag,prints,errors = gamma_run(gamma_command,log='PIPE')
    #
    if flag > 0:
       #
       x=None
       y=None
    else:
       #
       for line in prints.split(os.linesep):
           if "(int):" in line:
              outline = line.split(':')[1]
              # print(outline)
              tmp = outline.split(' ')
              tmp = [x for x in tmp if x]
              x = tmp[0]
              y = tmp[1]
              # print('X: %s Y: %s' % (x,y)
    return x,y            
###############################################################################
#
def read_slc_tab(intab):
    '''
    To import SLC_tab into python
    
    '''
    slcs, slcpars = [],[]
    with open(intab,'r') as fid:
       for cline in fid:
           cline = cline.split('\n')[0]
           slcs.append(cline.split()[0])
           slcpars.append(cline.split()[1])
    #
    return slcs,slcpars
#
def gamma_read_data(indata,nx,ny,mode='unw',dtype=">f4"):
    #
    if mode.upper() == "UNW":
       #
       fmt = pSAR.roipac.img_to_fmt(indata,nx,ny)
       #
       if fmt == 'float64':
          print(" Warning: this is float64 file. Convert to float32 now...")
          dtype = fmt
         
       data = np.fromfile(indata,count=-1,dtype=dtype,sep="")
       data = data.reshape([ny,nx]).astype('float')
       # 
    #
    return data 
###############################################################################
def gamma_s1a2burst(toppar):
    #
    print(toppar)
    fid = open(toppar,'r')
    for cline in fid:
        cline = cline.split('\n')[0]
        #
        if "number_of_bursts" in cline:
            # print(cline)
            tmp = cline.split()
            #
    #
    return int(tmp[1])
#
###############################################################################
#
def gamma_extract_value(inoff,keyword):
    #
    # It is better to keep output in string.
    # Format conversion can be considered when we need.
    #
    output = 0 
    fid = open(inoff,'r')
    for line in fid:
        if re.search(keyword, line):
           output = line.split()
           output = output[1:]
    #
    return output
###############################################################################
def s1tab_topath(intab):
    #
    topdir = os.path.dirname(os.path.abspath(intab))
    fid    = open(intab,'r')
    outfiles = []
    for cline in fid:
        cline = cline.split('\n')[0].split()
        if not os.path.exists(cline[0]):
           in_slc        = os.path.join(topdir,os.path.basename(cline[0]))
           in_slc_par    = os.path.join(topdir,os.path.basename(cline[1]))
           in_slc_toppar = os.path.join(topdir,os.path.basename(cline[2]))
        else:
           in_slc = os.path.abspath(cline[0])
           in_slc_par = os.path.abspath(cline[1])
           in_slc_toppar = os.path.abspath(cline[2])
        #
        swath = os.path.basename(in_slc).split('_')[1]
        outfiles.append([in_slc,in_slc_par,in_slc_toppar,swath])
    #
    outfiles = np.array(outfiles)
    return outfiles
###############################################################################    
def tab_to_dates(intab):
    """
       Read sar information from a tab file...
    """
    dList   = []
    ID      = 0
    #
    fid = open(intab,'r')
    for line in fid:
        #
        output = line.split()
        dList.append(strStruc(ID,output[0],output[1]))
        ID += 1
    #
    return dList
###############################################################################
def gamma_SLC_interp(sslc,mslc_par,sslc_par,off_par,resslc,resslc_par):
    #
    command_sTr = 'SLC_interp ' + sslc + ' ' + mslc_par + ' ' + \
                   sslc_par + ' ' + off_par + ' ' + resslc + ' ' + resslc_par
    return command_sTr
###############################################################################
#
def gamma_SLC_copy(inslc,range0,nsample,azi0,nline,outslc):
    #
    inslcpar=inslc+'.par'
    outslcpar=outslc+'.par'
    #
    #
    command_sTr=("SLC_copy %s %s %s %s - - %s %s %s %s" % (inslc,inslcpar,\
                               outslc,outslcpar,str(range0),str(nsample),\
                               str(azi0),str(nline)))
    flag,prints,errors = gamma_run(command_sTr)
    if os.path.exists(outslc):
        return True
    else:
        return False
###############################################################################
def read_offset_list(inlist):
    #
    #
    offList = []
    counter = 0
    if os.path.exists(inlist):
       #
       fid = open(inlist,'r')
       
       for line in fid:
          tmp = line.split()
          counter += 1
          offList.append(strOffsets(tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5]))
       #
       fid.close()
       #
    else:
       print(" " + inlist + " is not found.")
    # 
    return offList,counter
#
###############################################################################
#
def nw_est_info(inlist):
    #
    #
    offList, counter = read_offset_list(inlist)
    #
    datainfo = np.zeros((counter,5))
    listOffs = []
    for ni in range(counter):
        datainfo[ni,:] = [int(offList[ni].ID),int(offList[ni].MID),\
                int(offList[ni].SID),int(offList[ni].MDATE),\
                int(offList[ni].SDATE)]
        listOffs.append(offList[ni].OFFSETS)
    #
    return datainfo, listOffs
#
###############################################################################
#
def nw_offsets_est(datainfo,listOffs,master,npoly=3,cor_thres=0.15):
    #
    # Calculate the polynominal parameters relative to a given master...
    # all offsets reference points will be considered...
    #
    dates  = np.hstack(datainfo[:,[3,4]])
    udates = np.unique(dates)
    ndates = udates.shape[0]
    #
    master_ind = np.where(udates == master)[0][0]
    counter = 0
    #
    for clist in listOffs:
        wA,dx,dy,cor = offsets_to_A(clist,npoly=npoly,cor_thres=cor_thres)
        numref       = wA.shape[0]
        numpol       = wA.shape[1]
        cA           = np.zeros((numref,numpol*(ndates)))
        # 
        mID          = datainfo[counter,1]
        sID          = datainfo[counter,2]

        #
        if udates[mID] == master:
           #
           cA[:,mID*numpol:(mID+1)*numpol] = -1 * wA
           cA[:,sID*numpol:(sID+1)*numpol] = wA
        else:
           cA[:,mID*numpol:(mID+1)*numpol] = wA
           cA[:,sID*numpol:(sID+1)*numpol] = -1 * wA
        #
        if udates[mID] != master:
           dx = dx*-1
           dy = dy*-1
        #
        if counter == 0:
           A = cA
           DDx = dx * cor
           DDy = dy * cor
        else:
           A = np.concatenate((A,cA))
           DDx = np.concatenate((DDx, dx * cor))
           DDy = np.concatenate((DDy, dy * cor))
                  
        counter += 1
    #
    A = np.delete(A,np.s_[master_ind*numpol:(master_ind+1)*numpol],1)     
    #
    cofx, resx, nrnkx, sx = np.linalg.lstsq(A,DDx)
    cofy, resy, nrnky, sy = np.linalg.lstsq(A,DDy)
    #
    resx = np.std(np.dot(A,cofx) - DDx)
    resy = np.std(np.dot(A,cofy) - DDy)
    #
    # print(nrnkx,nrnky)
    #
    out_coef = {}
    #
    for ni in range(ndates):
        if ni < master_ind or master_ind == 0:
           out_coef[str(int(udates[ni]))] = np.array(\
                [cofx[ni*numpol:(ni+1)*numpol],cofy[ni*numpol:(ni+1)*numpol]])
        if ni > master_ind and master_ind != 0:
           out_coef[str(int(udates[ni]))] = np.array(\
            [cofx[(ni-1)*numpol:(ni)*numpol],cofy[(ni-1)*numpol:(ni)*numpol]])
    #
    udates = np.delete(udates, np.where(udates==master),0).astype('int')
    #
    return out_coef,udates
#
###############################################################################
def gamma_logrun(sTr,logfile):
    #
    with open(logfile,'w') as logid:
        gamma_run(sTr,log=logid)
    return None
#
def gamma_run(sTr,log="PIPE"):
    #
    print(' pGAMMA: %s' % sTr)
    #
    if isinstance(log,str) and log.upper() == "PIPE":
        log = subprocess.PIPE
    #
    output        = subprocess.Popen(sTr,stdout=log,stderr=log,shell=True)
    flag          = output.wait()
    prints,errors = output.communicate()
    if prints is not None:
       prints = pDATA.bytestoutf8(prints)
    #
    return flag,prints,errors
#
###############################################################################
#
def gamma_create_offset(m_par,s_par,OFF_par,algorithm=1,rlks=1,azlks=1,iflg=0):
    #
    # Create a processing control file with GAMMA
    command_sTr = 'create_offset ' + m_par + ' ' + s_par + ' ' + OFF_par + \
                  ' ' + str(algorithm) + ' ' + str(rlks) + ' ' + \
                  str(azlks) + ' ' + str(iflg)

    return command_sTr 
#
###############################################################################
#
def gamma_init_offset_orbit(m_par,s_par,OFF_par):
    command_sTr = 'init_offset_orbit ' +  m_par + ' ' + s_par + ' ' + \
                  OFF_par + ' - - 1'
    return command_sTr
#
###############################################################################
#
def gamma_init_offset(m_slc,s_slc,m_par,s_par,OFF_par,cor_thres=0.2,\
                      win_size=64):
    #
    command_sTr = 'init_offset ' + m_slc + ' ' + s_slc + ' ' + \
                  m_par + ' '+ s_par + ' ' + OFF_par + ' 1 1 - - - - ' + \
                  str(cor_thres) + ' ' + str(win_size) + ' ' + \
                  str(win_size) + ' 1 '
    #              
    return command_sTr
#
###############################################################################
#
def gamma_offset_pwr(m_slc,s_slc,m_par,s_par,OFF_par_basename,\
                     rwin=32,azwin=32,nr=128,naz=128,ovr=4,thres=0.2):
    #
    OFF_par = OFF_par_basename + '.off'
    offs    = OFF_par_basename + '.offs'
    ccp     = OFF_par_basename + '.ccp'
    offsets = OFF_par_basename + '.offsets'
    #
    #
    command_sTr = 'offset_pwr ' + m_slc + ' ' + s_slc + ' ' + \
                  m_par + ' ' + s_par + ' ' + OFF_par + ' ' + offs + ' ' + \
                  ccp   + (" %d %d " % (rwin,azwin)) + offsets + \
                  (" %d " % (ovr)) + (" %d %d" % (nr,naz)) + (" %f " % thres)
    #              
    return command_sTr,offsets  
#########################################################
def gamma_offset_fit(offs,ccp,OFF_par,npoly=3,cor_thres=0.2):
    #
    OFF_par_basename = os.path.basename(OFF_par)
    outputs = OFF_par_basename.split('.')
    OFF_par_basename = outputs[0]
    coffs            = os.path.dirname(os.path.abspath(OFF_par)) + '/' + \
                       OFF_par_basename + '.coffs'
    coffsets         = os.path.dirname(os.path.abspath(OFF_par)) + '/' + \
                       OFF_par_basename + '.coffsets'
    #
    command_sTr = 'offset_fit ' + offs + ' ' + ccp + ' ' + \
                  OFF_par + ' ' + coffs + ' ' + coffsets + ' ' + \
                  str(cor_thres) + ' ' + str(npoly)
    #
    return command_sTr,coffsets
#########################################################
# 
def offsets_to_A(offsets,npoly=3,cor_thres=0.2):
    #
    refpoints = np.loadtxt(offsets)
    #
    x   = refpoints[:,0]
    y   = refpoints[:,1]
    dx  = refpoints[:,2]
    dy  = refpoints[:,3]
    cor = refpoints[:,4]
    x   = x[np.where(cor >= cor_thres)]
    y   = y[np.where(cor >= cor_thres)]
    dx  = dx[np.where(cor >= cor_thres)]
    dy  = dy[np.where(cor >= cor_thres)]
    cor = cor[np.where(cor >= cor_thres)]
    #
    dsize = x.shape[0]
    x   = np.reshape(x,(dsize,1))
    y   = np.reshape(y,(dsize,1))
    dx  = np.reshape(dx,(dsize,1))
    dy  = np.reshape(dy,(dsize,1))
    cor = np.reshape(cor,(dsize,1))
    cor = cor / np.average(cor) 
    #
    if npoly == 3:
       wA = np.concatenate(((x * 0 + 1) * cor, x * cor, y * cor), axis = 1)
    elif npoly == 4:
       wA = np.concatenate(((x * 0 + 1) * cor, x * cor, y * cor, x * y * cor),\
                           axis = 1)
    # 
    return wA, dx, dy, cor
#########################################################
def est_polycoef(offsets,npoly=3,cor_thres=0.2):
    #
    wA,dx,dy,cor = offsets_to_A(offsets,npoly=npoly,cor_thres=cor_thres);
    #
    cofx, resx, nrnkx, sx = np.linalg.lstsq(wA,dx * cor)
    cofy, resy, nrnky, sy = np.linalg.lstsq(wA,dy * cor)
    resx = np.std(np.dot(wA,cofx) - dx * cor)
    resy = np.std(np.dot(wA,cofy) - dy * cor)
    #
    return cofx,cofy,resx,resy
###############################################################################
def gamma_offset_add(offpar1,offpar2,offpar):
    #
    cmd_sTr = "offset_add %s %s %s" % (offpar1,offpar2,offpar)
    return cmd_sTr
###############################################################################
def gamma_coregionquality(m_slc,s_slc,OFF_par_basename=None,\
                          offdir='coregistration'):
    #
    m_par = m_slc+'.par'
    s_par = s_slc+'.par'
    mbase = os.path.basename(m_slc).split('.')[0]
    sbase = os.path.basename(s_slc).split('.')[0]
    #
    if not os.path.exists(offdir):
       os.makedirs(offdir)
    #
    if OFF_par_basename is None:
       OFF_par_basename = os.path.join(offdir,mbase+'-'+sbase+'-QUAD')
    #
    offsetspar = OFF_par_basename+'.offsets'
    rmu,rsig,azmu,azsig = 0.,0.,0.,0.
    # 
    flag = True
    if os.path.exists(offsetspar):
        try:
            data = np.loadtxt(offsetspar)
            if data.shape[0] < 30:
                flag = True
        except:
            flag = False
            
    if flag:
       #
       logid = open(OFF_par_basename+'.log','w')
       gamma_est_offsets(m_slc,s_slc,m_par,s_par,OFF_par_basename,\
                      init_offsets=False,npoly=3,algorithm=1,\
                      rlks=1,azlks=1,iflg=0,cor_thres=0.2,logid=logid)
       logid.close()
    #
    offsets = np.loadtxt(OFF_par_basename+'.offsets')
    #
    if offsets.shape[0] < 5:
        rmu,rsig = 99,99
        azmu,azsig = 99,99
    else:
        rmu,rsig = np.median(offsets[:,2]),np.std(offsets[:,2])
        azmu,azsig = np.median(offsets[:,3]),np.std(offsets[:,3])
   
    #
    return rmu,rsig,azmu,azsig
#
###############################################################################
def gamma_slcresample_demassist(m_slc,s_slc,m_par,s_par,OFF_par_basename,\
                      rslc,rslc_par,init_offsets=True,npoly=3,algorithm=1,\
                      rlks=1,azlks=1,iflg=0,cor_thres=0.15,ovr=4,\
                      rwin=64,azwin=64,nr=128,naz=128,demrdc=None,mmli=None,\
                      smli=None,iters=3):
    #
    mmli_par = mmli+'.par'
    smli_par = smli+'.par'
    lt0      = OFF_par_basename+'.lt0'
    lt       = OFF_par_basename+'.lt'
    diffpar  = OFF_par_basename+'.diff_par'
    logid    = OFF_par_basename+'.log'
    lt0_offs = OFF_par_basename+'.lt0.offs'
    lt0_ccp  = OFF_par_basename+'.lt0.ccp'
    #
    lt0_offsets = OFF_par_basename+'.lt0.offsets'
    mmli_par  = mmli+'.par'
    mmli_lt   = OFF_par_basename+'_lt.mli'
    #
    logid = open(logid,'a')
    #
    mliinfo  = gamma_mlipar(mmli_par)
    mwidth   = int(mliinfo["range_samples"][0])
    smliinfo = gamma_mlipar(smli_par)
    swidth   = int(smliinfo["range_samples"][0])
    #
    rdc_trans = "rdc_trans %s %s %s %s" % (mmli_par,demrdc,smli_par,lt0)
    flag,prints,errs = gamma_run(rdc_trans,log=logid)
    geocode = "geocode %s %s %d %s %d - 2 0 - - 4" % (lt0,mmli,mwidth,\
                                                      mmli_lt,swidth)
    flag,prints,errs = gamma_run(geocode,log=logid)
    create_diff_par = "create_diff_par %s - %s 1 0" % (smli_par,diffpar)
    flag,prints,errs = gamma_run(create_diff_par,log=logid)
    init_offsetm = 'init_offsetm %s %s %s' % (mmli_lt,smli,diffpar)
    flag,prints,errs = gamma_run(init_offsetm,log=logid)
    offset_pwrm = 'offset_pwrm %s %s %s %s %s - - %s' % (mmli_lt,smli,\
                                        diffpar,lt0_offs,lt0_ccp,lt0_offsets)
    flag,prints,errs = gamma_run(offset_pwrm,log=logid)
    offset_fitm = 'offset_fitm %s %s %s - - - %d' % (lt0_offs,\
                                                     lt0_ccp,diffpar,npoly)
    flag,prints,errs = gamma_run(offset_fitm,log=logid)
    gc_map_fine = 'gc_map_fine %s %d %s %s' % (lt0,mwidth,diffpar,lt)
    flag,prints,errs = gamma_run(gc_map_fine,log=logid)
    #
    SLC_interp_lt = 'SLC_interp_lt %s %s %s %s %s %s - %s %s' % \
                  (s_slc,m_par,s_par,lt,mmli_par,smli_par,rslc,rslc_par)
    flag,prints,errs = gamma_run(SLC_interp_lt,log=logid)        
    #
    for citer in range(iters):
      #
      ROFF_par_basename = OFF_par_basename+('_%d' % citer)
      gamma_est_offsets(m_slc,rslc,m_par,rslc+'.par',ROFF_par_basename,\
                      init_offsets=False,npoly=3,algorithm=1,\
                      rlks=rlks,azlks=azlks,iflg=0,\
                      cor_thres=cor_thres,ovr=ovr,\
                      rwin=rwin,azwin=azwin,nr=nr,naz=naz)
      SLC_interp_lt = 'SLC_interp_lt %s %s %s %s %s %s %s %s %s' % \
                  (s_slc,m_par,s_par,lt,mmli_par,smli_par,\
                   ROFF_par_basename+'.off',rslc,rslc_par)
      flag,prints,errs = gamma_run(SLC_interp_lt,log=logid)
    #
    # Close the log file.
    #
    logid.close()
    if (os.path.exists(rslc) and os.path.exists(rslc_par)):
        return True
    else:
        return False
#    
###############################################################################
#    
def gamma_est_offsets(m_slc,s_slc,m_par,s_par,OFF_par_basename,\
                      init_offsets=True,npoly=3,algorithm=1,update=False,\
                      rlks=1,azlks=1,iflg=0,cor_thres=0.15,ovr=4,\
                      rwin=64,azwin=64,nr=64,naz=64,logid=None):
    #
    OFF_par = OFF_par_basename + '.off'
    #
    command_create_off_par = gamma_create_offset(m_par,s_par,OFF_par,\
                        algorithm=algorithm,rlks=rlks, azlks=azlks, iflg=iflg)
    #
    #  
    if not os.path.exists(OFF_par) or os.stat(OFF_par).st_size < 2:
       #
       flag,prints,errs = gamma_run(command_create_off_par,log=logid)
       if flag != 0:
          print(prints)
    #
    if init_offsets:
       #
       # If checking the quality of resampled slave is required, 
       # inital_offsets are not needed to estimate....
       #
       command_init_offset_orbit = gamma_init_offset_orbit(m_par,s_par,OFF_par)
       flag,prints,errs = gamma_run(command_init_offset_orbit,log=logid)
       #
       command_init_offsets = gamma_init_offset(m_slc,s_slc,m_par,s_par,\
                                                 OFF_par,cor_thres=cor_thres)
       flag,prints,errs = gamma_run(command_init_offsets,log=logid) 
    #
    #
    command_offset_pwr,offsets = gamma_offset_pwr(m_slc,s_slc,m_par,s_par,\
                            OFF_par_basename,ovr=ovr,nr=nr,naz=naz,rwin=rwin,\
                            azwin=azwin) 
    #
    offs = OFF_par_basename + '.offs'
    ccp  = OFF_par_basename + '.ccp'
    if not os.path.exists(offs) or os.stat(offs).st_size < 1:
       flag,prints,errs = gamma_run(command_offset_pwr,log=logid)
    #
    command_offsets_fit,coffsets = gamma_offset_fit(offs,ccp,OFF_par,\
                                            npoly=npoly,cor_thres=cor_thres)
    flag,prints,errs = gamma_run(command_offsets_fit,log=logid)
    #
    return flag,prints,errs
