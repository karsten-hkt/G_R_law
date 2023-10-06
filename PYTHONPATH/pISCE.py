#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 08:35:13 2017

@author: wafeng
"""

import pSAR
import xml.etree.ElementTree as ET
#
###############################################################################
def iscexml(inxml):
    #
    info, keys = pSAR.roipac.roipac_info()
    # DOMTree    = xml.dom.minidom.parse(inxml)
    Tree       = ET.parse(inxml)
    root       = Tree.getroot()
    #
    for child in root.iter('component'):
        # print(child.tag,child.get('name'))
        if child.get('name').upper() == 'COORDINATE1':
           for cid in child.iter('property'):
              #
              if cid.get('name').upper() == 'DELTA':
                 subid = cid.find('value')
                 info['X_STEP'] = subid.text
              if cid.get('name').upper() == 'SIZE':
                 subid = cid.find('value')
                 info['WIDTH'] = subid.text 
              if cid.get('name').lower() == 'startingvalue':
                 subid = cid.find('value')
                 info['X_FIRST'] = subid.text     
        if child.get('name').upper() == 'COORDINATE2':
           for cid in child.iter('property'):
              #
              if cid.get('name').upper() == 'DELTA':
                 subid = cid.find('value')
                 info['Y_STEP'] = subid.text
              if cid.get('name').upper() == 'SIZE':
                 subid = cid.find('value')
                 info['FILE_LENGTH'] = subid.text 
              if cid.get('name').lower() == 'startingvalue':
                 subid = cid.find('value')
                 info['Y_FIRST'] = subid.text               
    return info