#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 24 08:08:34 2017

@author: wafeng
"""
import netCDF4
import numpy as np
#
def import_grd(in_grd):
    #
    grdid = netCDF4.Dataset(in_grd,'r')
    db    = grdid.variables
    z = db['z']
    # a necessary operation to convert netCDF4 array to np.ndarray
    #
    z = np.copy(z)
    z = np.flipud(z)
    return z
