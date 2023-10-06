# -*- coding: utf-8 -*-
"""
Created on Mon Nov  1 06:21:14 2021

@author: wafeng
"""

from . import tde
#
def tde2G(tri,sx,sy,sz,pr=0.25,r1=0,r2=90):
    #
    # r1, rake1, in degree
    # r2, rake2, in degree
    # normally, we need two orthogonal slip vectors in practice.
    #