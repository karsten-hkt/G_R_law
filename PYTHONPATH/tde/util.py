# -*- coding: utf-8 -*-
"""
Created on Mon Nov  1 06:21:14 2021

@author: wafeng
"""
import os
from . import tde
import datetime as dt
import numpy as np
import multiprocessing as mpg
from functools import partial
# from numba import jit
#
# global tris,sx,sy,sz,rad_fac
#
# @jit(nopython=True)
def tri2out(tri,outname,un='NULL',ul='NULL'):
    #
    #
    fmt = ['%f' for i in range(12)]
    fmt = fmt + ['\n']
    fmt = ' '.join(fmt)
    with open(outname,'w') as fid:
        #
        fid.write('# Updated on %s\n' % dt.datetime.now().strftime('%Y-%m-%d'))
        fid.write('# UTM ZONE: %s %s\n' % (un,ul))
        fid.write('# Number of faults: %d  ModelType: TDE\n' % tri.shape[0])
        fid.write('# x1(km) y1(km) z1(km) x2(km) y2(km) z2(km) x3(km) y3(km) z3(km) s_slip(m) d_slip(m) o_slip(m)\n')
        for i in range(tri.shape[0]):
            #
            fid.write(fmt % (tri[i,0],tri[i,1],tri[i,2],\
                             tri[i,3],tri[i,4],tri[i,5],\
                             tri[i,6],tri[i,7],tri[i,8],\
                             tri[i,9],tri[i,10],tri[i,11]))
            #
        #
    #
    if os.path.exists(outname):
        return True
    else:
        return False
def import_tri(inname):
    #
    un = 'NULL'
    ul = 'NULL'
    tri = []
    append = tri.append
    #
    with open(inname,'r') as fid:
        for cline in fid:
            #
            cline = cline.split('\n')[0]
            if 'UTM' in cline:
                tmp = cline.split('ZONE:')[1].split()
                un = tmp[0]
                ul = tmp[1]
            if 'Number of faults:' in cline:
                tmp = cline.split('faults:')[1].split()
                nfault = int(tmp[0])
                #
                fid.readline()
                for i in range(nfault):
                    #
                    tmp = fid.readline()
                    append([float(ctmp) for ctmp in tmp.split('\n')[0].split()])
                
    return np.array(tri),un,ul
#
# @jit(nopython=True)
def tde2G_individual(i,tris,sx,sy,sz,rake=0,pr=0.25,ts=0,rad_fac=np.pi/180):
    #
    #
    # print(" ... calculating %d TDE..." % i)
    xyz = np.reshape(tris[i,0:9],[3,3])
    u = tde.calc_tri_displacements(sx=sx, sy=sy, sz=sz, \
                               x=xyz[:,0], y=xyz[:,1], z=xyz[:,2], 
                               pr=pr, ts=ts,\
                               ss=np.cos(rake*rad_fac), ds=np.sin(rake*rad_fac))
    return u
#
# @jit(nopython=True)
def tde2G(tris,sx,sy,sz,pr=0.25,r1=0,r2=90,njob=4):
    #
    # r1, rake1, in degree
    # r2, rake2, in degree
    # normally, we need two orthogonal slip vectors in practice.
    #
    if True:
      rad_fac = np.pi/180
      G_E_r1 = np.zeros([sx.shape[0],tris.shape[0]])
      G_N_r1 = np.zeros([sx.shape[0],tris.shape[0]])
      G_U_r1 = np.zeros([sx.shape[0],tris.shape[0]])
      if r2 != r1:
        G_E_r2 = np.copy(G_E_r1)
        G_N_r2 = np.copy(G_N_r1)
        G_U_r2 = np.copy(G_U_r1)
      else:
        G_E_r2 = None
        G_N_r2 = None
        G_U_r2 = None
      #
      for i in range(tris.shape[0]):
          #
          print(" ... now calculating displacements for a unit slip at %d fault" % i)
          u = tde2G_individual(i,tris,sx,sy,sz,rake=r1,rad_fac=rad_fac)
          G_E_r1[:,i] = u['x'] * -1
          G_N_r1[:,i] = u['y'] * -1
          G_U_r1[:,i] = u['z']
          if r2 != r1:
             u = tde2G_individual(i,tris,sx,sy,sz,rake=r2,rad_fac=rad_fac)
             G_E_r2[:,i] = u['x'] * -1
             G_N_r2[:,i] = u['y'] * -1
             G_U_r2[:,i] = u['z']
      #
    #
    return G_E_r1,G_N_r1,G_U_r1,G_E_r2,G_N_r2,G_U_r2
