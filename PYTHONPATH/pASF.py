#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 23:16:21 2019

@author: wafeng
"""
import os
import sys
#
def asfapi(searchpoly,platform='ALOS',start=None,end=None,\
           relativeOrbit=None,frame=None,processingLevel='L1.0',\
           output='asf-test',verbose=True,of='CSV'):
    # 
    # searchpoly is an numpy array of N * 2 
    #
    output = output 
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
    #
    goCMD = searchstr + \
            polygonstr + \
            platformstr + \
            startstr + \
            trackstr + \
            endstr + \
            framestr + \
            processingLevelstr+\
            '\&output=%s > %s' % (of,output)
    #
    if verbose:
        print(goCMD)
    #
    os.system(goCMD)
    if os.path.exists(output):
        return output
    else:
        return None