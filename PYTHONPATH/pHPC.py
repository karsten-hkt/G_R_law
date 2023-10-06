#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 28 09:29:42 2017

@author: wafeng
"""
import os
import sys
#
#
def dorvalhpc2servers():
    '''
    This only works for the HPC system in Dorval, Quebec. 
    
    '''
    servers = []
    if "GECOSHEP_TMP" in os.environ:
      tmp_path = os.environ["GECOSHEP_TMP"]
      server_list = tmp_path+'/pnet.hosts'
      #
      # servers = []
      #
      with open(server_list,'r') as fid:
           servers = fid.readlines()
      #
      servers = [cserver.split('\n')[0] for cserver in servers]
      servers.sort()
      #
    return servers

