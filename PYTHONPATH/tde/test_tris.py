import sys
import numpy as np
import matplotlib.pyplot as plt
import tde

from importlib import reload
reload(tde)

# tri dislocation parameters
pr = 0.25
ss = -1.
ts = 0.
ds = 0.
N = 30

sx, sy, sz = np.meshgrid( np.linspace(0, 100, N), np.linspace(0, 100, N), 0)

sxr = sx.ravel(order='F')
syr = sy.ravel(order='F')
szr = sz.ravel(order='F')

X = np.array([40., 60., 40.])
Y = np.array([50., 50., 30.])
Z = np.array([0., 0., 20.])

S = tde.calc_tri_strains(sxr, syr, szr, X, Y, Z, pr, ss, ts, ds)

U = tde.calc_tri_displacements(sxr, syr, szr, X, Y, Z, pr, ss, ts, ds)


plt.figure()

plt.gca().add_patch(plt.Polygon( np.array([X,Y]).T))

plt.quiver(sxr, syr, U['x']*0.5, U['y']*0.5)

plt.show()
