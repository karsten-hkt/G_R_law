#!/usr/bin/env python
#
import os
import sys
import numpy as np
#
def losvec(inc,azi):
    #
    rad_factor = np.pi / 180.
    #
    n_rng =      np.sin(azi * rad_factor) * np.sin(inc * rad_factor)
    e_rng = -1 * np.cos(azi * rad_factor) * np.sin(inc * rad_factor)
    u_rng =      np.cos(inc * rad_factor)
    return n_rng, e_rng, u_rng
#
helpstr = \
      '''
         %s <inc> <azi>
         +++++++++++++++++++++++++++++++
         To calculate projection vectors along SAR viewing geometry 
         based on given incidence angle (inc) and azimuth (azi)
         Note that inc and azi are all in degree.
         
         Written by Wanpeng Feng, @ SYSU, 2018-12-20

      '''
if __name__ == "__main__":
   if len(sys.argv)<3:
     print(helpstr % sys.argv[0])
     sys.exit(-1)
   #
   inc = float(sys.argv[1])
   azi = float(sys.argv[2])
   #
   #
   print("INC = %f and AZI =  %f" % (inc,azi))
   #
   e_rng, n_rng, u_rng = losvec(inc,azi)
   #
   print("Vector(e/n/u)_for_RNG_changes: %5.5f %5.5f %5.5f" % (e_rng,n_rng,u_rng))
