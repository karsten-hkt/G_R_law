#! /usr/bin/env python
#
# A Python Module to handle different data in different formats
# The data formats supported in this module include ndk(GCMT), kml, geojson,
# specified coseismic GPS data. This can be growing quickly and bigger 
# with time, which hope to be very supportive in the future academic 
# activities...
# 
# ALOS2 CEOS format can also be supported since May 2017. 
# led2info() can return information 
# Created and Updated by Wanpeng Feng, @CCRS/NRCan, since 2016-11-14
#
# Since 2017-05-25, all python scripts, modules and any other python related 
# codes will be implemented with python 3.x...
# Some bugs due to the huge changes cannot be avoided. 
# Btw, a module, math was also obsoluted since this version. All funtions in
# math will be replaced with numpy. Noted by Wanpeng Feng, @CCRS/NRCan, 
# 2017-05-25
#
###############################################################################
import numpy as np
import pSAR
import os
import sys
import subprocess
import datetime
try:
  import geojson
except:
  1
#
import csv
from lxml import etree
from datetime import datetime as dt
import xml.dom.minidom
import glob
import xml.etree.ElementTree as ET
import zipfile
import shutil
from datetime import timedelta
#
###############################################################################
def alos2zip2mode(inzip):
    #
    # zipid = open(inzip,'r')
    ledfilename = None
    modetype = None
    if zipfile.is_zipfile(inzip):
       zf = zipfile.ZipFile(inzip, 'r')
       #
       for cfile in zf.namelist():
           if 'LED' in cfile:
               ledfilename = cfile
           #
       #
    #
    if ledfilename is not None:
       tmp = ledfilename.split('_')[0].split('-')[-1]
       modetype = tmp[0:3]
    #
    return modetype,ledfilename
def bytestoutf8(t):
    '''
    To convert a bytes-like variable to utf8
    '''
    if isinstance(t,bytes):
        return t.decode("utf-8")
    return t
#
###############################################################################
#
def led2info(ledfile):
    '''
    '''
    with open(ledfile,'r') as fid:
      fid.seek(780,os.SEEK_SET)
      time = fid.read(30).strip()
      fid.seek(156,0)
      info = fid.read(40)
      #
      #  mission number and data level
      #
      fid.seek(48,0)
      info = fid.read(16)
      #
      # print(info)
      missionnum   = info[2]
      processlevel = info[7]
    #
    # Updated by Wanpeng Feng, @NRCan, 2017-10-16
    # Two delimieters, - and _ can both be considered now for names used 
    # in ALOS2 and ALOS.
    #
    tmps = os.path.basename(ledfile).split('.')[0].replace('-','_').split('_')
    #
    if len(tmps)< 2:
        corzip = ledfile.replace('.led','.zip')
        if os.path.exists(corzip):
            corled = alos2ledfile(corzip,outled=None,orgled=True)
            tmps = corled.split('-')
    #
    if missionnum == "2":
       for cfile in tmps:
         if "ALOS" in cfile:
            info = cfile 
       orbnum = info[5:10]
       # print(orbnum)
       track  = ((int(float(orbnum)) - 989)*14) % 207 + 1
       frame  = info[10:]
    else:
       #
       for cfile in tmps:
         if "ALP" in cfile:
            info = cfile
       orbnum = info[6:11]
       track  = (int(float(orbnum))*46+84) % 671 + 1
       frame  = info[11:]
    return {"time":time,"track":str(track),"frame":frame,\
            "mission":missionnum,"Processlevel":processlevel,\
            'yyyymmdd':time[:8]}
#
###############################################################################
def alos2zip2file(inzip,givenstr):
    #
    # extract files from an ALOS-2 zip file
    #
    try:
      zf            = zipfile.ZipFile(inzip)
      zip_files     = zf.namelist()
      #
      for cname in zip_files:
         #
         if givenstr in cname:
            #
            zf.extract(cname)
      return True  
    except:
      pass
      return False 
#
def alos2zip2imgs(inzip):
    #
    HHs = []
    VVs = []
    HVs = []
    VHs = []
    try:
      zf = zipfile.ZipFile(inzip)
      zip_files     = zf.namelist()
      for cname in zip_files:
          if 'IMG' in cname and 'HH' in cname:
              HHs.append(cname)
          if 'IMG' in cname and 'VV' in cname:
              VVs.append(cname)
          if 'IMG' in cname and 'HV' in cname:
              HVs.append(cname)
          if 'IMG' in cname and 'VH' in cname:
              VHs.append(cname)   
      #
      
    except:
      #
      pass
      #
    return HHs,VVs,HVs,VHs
    #  
def alos2ledfile(inzip,outled=None,orgled=False):
    '''
    To extract LED files from zip files direclty...
    '''
    outled = []
    try:
      zipbasename   = os.path.basename(inzip).split('.')[0]
      zf            = zipfile.ZipFile(inzip)
      zip_files     = zf.namelist()
      czip_dir_name = os.path.dirname(os.path.abspath(inzip))
      #
      for cname in zip_files:
         #
         # root_name    = os.path.basename(cname)
         dir_name     = os.path.dirname(cname)
         cur_dir_name = os.getcwd()
         #
         if "LED" in cname:
            #
            if (outled is None or len(outled) == 0):
               outled = zipbasename + ".led"
            #
            output_manifest = os.path.join(czip_dir_name,outled)
            zf.extract(cname)
            # 
            if orgled:
                outled = cname
            else:
                if (outled is None or len(outled) == 0):
                   outled = zipbasename + ".led"
                   #
                output_manifest = os.path.join(czip_dir_name,outled)
                zf.extract(cname)
                #
                shutil.move(cname,output_manifest)
                print(" %s : %s" % (inzip,os.path.basename(output_manifest)))
                if cur_dir_name != dir_name:
                   shutil.rmtree(dir_name,ignore_errors=False)
              
    except:
      #print(' Error during opening %s' % inzip)
      pass
    return outled
#    
###############################################################################
def auig2csv2data(incsv):
    #
    # Read ALOS2 data into separate KML with beam and track information
    # by Wanpeng Feng, @CCRS/NRCan, 2017-05-02
    #
    outdata = []
    counter = 0
    atts    = []
    with open(incsv,'r') as csvid:
        csvreader = csv.reader(csvid,delimiter=',')
        for row in csvreader:
            #
            counter += 1
            if counter == 1:
                atts = np.array(row).T
            if (len(row) > 0 and counter > 1):
               nrow = np.array(row).T
               outdata.append(nrow)
        #
    #
    return np.array(outdata),atts
#
#   
def alos2ledinfo(ledfile):
    '''
    To read ALOS2 metadata from a LED-* file
    '''
    #
    fid = open(ledfile,'r')
    fid.seek(780,os.SEEK_SET)
    time = fid.read(30).strip()
    fid.seek(64,0)
    info = fid.read(4)
    #
    #  mission number and data level
    #
    fid.seek(48,0)
    info = fid.read(16)
    missionnum   = info[2]
    processlevel = info[7]
    #
    fid.close()
    #
    #print(os.path.basename(ledfile))
    info = os.path.basename(ledfile).split('.')[0].split('-')
    tinfo = None
    for cinfo in info:
        if "ALPS" in cinfo:
            tinfo = cinfo
    if tinfo is not None:
        info = tinfo
    else:
        info = info[1]
    # print(missionnum)
    if missionnum == "2":
       # print(info)
       #
       orbnum = info[5:10]
       # print(orbnum)
       track  = ((int(float(orbnum)) - 989)*14) % 207 + 1
       frame  = info[10:]
    else:
       # print(info,ledfile)
       orbnum = info[6:11]
       track  = (int(float(orbnum))*46+84) % 671 + 1
       frame  = info[11:]
       
    return {"time":time,"track":str(track),"frame":frame,"mission":missionnum,\
            "Processlevel":processlevel,'yyyymmdd':time[:8]}
#
############################################################################### 
#
def asfapifromled(cled,roi=None,outcsv=None):
    #
   info = led2info(cled)
   strtime = info['time']
   track = info['track']
   frame = info['frame']
   if frame[0] == '0':
      frame = frame[1::]
   #
   if roi is not None:
      polygons = pSAR.roipac.ext2polygon(roi)
   else:
      polygons = None
   #
   cdt = datetime.datetime(int(strtime[0:4]),int(strtime[4:6]),\
                          int(strtime[6:8]),hour=int(strtime[8:10]),\
                          minute=int(strtime[10:12]),\
                          second=int(strtime[12:14]),\
                          microsecond=int(strtime[14::]))
   # 
   # print(dt.strftime("%Y-%m-%dT%H:%M:%S.%f UTC"),strtime)
   #
   dtstart = cdt - datetime.timedelta(minutes=100)
   dtend   = cdt + datetime.timedelta(minutes=100)
   startstr = dtstart.strftime('%Y-%m-%dT%H:%M:%SUTC')
   endstr = dtend.strftime('%Y-%m-%dT%H:%M:%SUTC')
   #
   if outcsv is None:
      outcsv = os.path.dirname(os.path.abspath(cled))+'/API-' + \
            os.path.basename(cled).split('.led')[0]+'.csv'
   #
   outcsv = asfapi(polygons,platform='ALOS',start=startstr,end=endstr,\
             relativeOrbit=track,frame=None,outcsv=outcsv)
   return outcsv
#
###############################################################################
#
def asfapi(searchpoly,platform='ALOS',start=None,end=None,\
           relativeOrbit=None,frame=None,processingLevel='L1.0',\
           outcsv='asf-test.csv',verbose=True):
    # 
    # searchpoly is an numpy array of N * 2 
    #
    searchstr='curl https://api.daac.asf.alaska.edu/services/search/param?'
    if searchpoly is not None:
       polygonstr = 'polygon='+','.join(searchpoly.ravel().astype('str'))
    else:
       polygonstr = ''  
    if platform is not None:
       if len(polygonstr) > 1:
          platformstr = '\&platform=' + platform
       else:
          platformstr = 'platform=' + platform   
    else:
       platformstr = ''
    if start is not None:
        startstr = '\&start='+start
    else:
        startstr = ''
    if end is not None:
        endstr = '\&end='+end
    else:
        endstr = ''    
    if relativeOrbit is not None:
        trackstr = '\&relativeOrbit='+relativeOrbit
    else:
        trackstr = ''
    if frame is not None:
        framestr = '\&frame=' + frame
    else:
        framestr = ''
    # This should be platform-related.
    # 
    if processingLevel is not None:
        processingLevelstr = '\&processingLevel=%s' % processingLevel
    else:
        processingLevelstr = ''
    goCMD = searchstr + \
            polygonstr + \
            platformstr + \
            startstr + \
            trackstr + \
            endstr + \
            framestr + \
            processingLevelstr+\
            '\&output=csv > %s' % outcsv
    #
    if verbose:
        print(goCMD)
    os.system(goCMD)
    if os.path.exists(outcsv):
        return outcsv
    else:
        return None
    
def asfcsv(incsv):
    #
    atts = []
    data = []
    with open(incsv,'r') as fid:
        atts = fid.readline().split('\n')[0].split(',')
        for cline in fid:
            tmp = cline.split('\n')[0].split(',')
            tmp = [ctmp.split('"')[1] for ctmp in tmp]
            data.append(tmp)
    #
    return data,atts
