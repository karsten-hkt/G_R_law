#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 29 09:21:51 2018

@author: wafeng
"""
import os
import glob
import pSAR
#
def asa_findorb(starttime,endtime,orbdir=None):
    #
    # return a fullpath of orbit file for a single envisat data
    #
    if orbdir is None:
        orbdir = os.environ['VOR_DIR']
    #
    st_time_data = pSAR.ts.timestr2jd(starttime,fmt='%Y%m%dT%H:%M:%S.%f')
    sp_time_data = pSAR.ts.timestr2jd(endtime,fmt='%Y%m%dT%H:%M:%S.%f')
    #
    outorb = []
    files = glob.glob(orbdir+'/DOR_VOR*')
    for cfile in files:
        tmp = os.path.basename(cfile).split('_')
        st_time = pSAR.ts.timestr2jd(tmp[4]+'T'+tmp[5],fmt='%Y%m%dT%H%M%S')
        sp_time = pSAR.ts.timestr2jd(tmp[6]+'T'+tmp[7],fmt='%Y%m%dT%H%M%S')
        if st_time <=  st_time_data and sp_time >= sp_time_data:
           outorb = cfile
           break
    #
    return outorb
