#!/usr/bin/python3
import matplotlib.pyplot as plt
from matplotlib import cm, _cm, colors, gridspec, rcParams
from matplotlib.ticker import NullLocator, NullFormatter
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from numpy import pi, sin, cos
from invert_angles import G_angles_q12
from my_cms import husl_hot

# Parameters
N = 2**8
epsilon = 1

Theta, Phi = np.mgrid[0:pi/2:complex(0, N), 0:2*pi:complex(0, 2*N + 1)]
# Theta should be in the open interval (0, pi/2)
Theta = Theta[1:-1,:]
Phi = Phi[1:-1,:]
# Radial coordinates for a Lambert azimuthal equal-area projection
R = 2*cos(pi/2 - Theta/2)

X = sin(Theta) * cos(Phi)
Y = sin(Theta) * sin(Phi)
Z = cos(Theta)

angles = np.array([cos(Theta), Phi])

colorfunction=G_angles_q12(angles, epsilon)

norm = colors.Normalize()

#cmap = colors.LinearSegmentedColormap(name='cubehelixhot',
#        segmentdata=_cm.cubehelix(gamma=1, s=.9, r=.2, h=1.0))

cmap = husl_hot

fontsize = 72
fig = plt.figure(figsize=(9,12))
# fig = plt.figure()
gs = gridspec.GridSpec(2, 1, height_ratios=[3, 2])#, width_ratios=[5, 2])
fig.suptitle(r'$\epsilon=' + str(epsilon) + '$', fontsize=fontsize)

ax1 = plt.subplot(gs[0], projection='polar')
# ax1 = plt.subplot(projection='polar')
# Plot using an equal-area azimuthal projection
ax1.pcolormesh(Phi, R, colorfunction, cmap=cmap, norm=norm,
               shading='gouraud', rasterized=True)
ax1.set_ylim(0, np.max(R))
ax1.xaxis.set_major_formatter(NullFormatter())
ax1.yaxis.set_major_formatter(NullFormatter())
# Having problems with my custom colormap and facecolors.
ax2 = plt.subplot(gs[1], projection='3d')
surf = ax2.plot_surface(X, Y, Z,  rstride=1, cstride=1,
                facecolors=cmap(norm(colorfunction)), shade=False)
surf.set_rasterized(True)
ax2.set_zlim3d(0, 1)
ax2.w_zaxis.set_major_formatter(NullFormatter())
ax2.w_zaxis.set_major_locator(NullLocator())
ax2.set_xlim3d(-1, 1)
ax2.w_xaxis.set_major_formatter(NullFormatter())
ax2.w_xaxis.set_major_locator(NullLocator())
ax2.set_ylim3d(-1, 1)
ax2.w_yaxis.set_major_formatter(NullFormatter())
ax2.w_yaxis.set_major_locator(NullLocator())
'''
ax2.set_xlabel(r'$z$', fontsize=fontsize)
ax2.set_ylabel(r'$x$', fontsize=fontsize)
ax2.set_zlabel(r'$y$', fontsize=fontsize)
'''
m = cm.ScalarMappable(cmap=cmap, norm=norm)
m.set_array(colorfunction)
cbar = plt.colorbar(m, ticks=NullLocator(), format=NullFormatter())
cbar.solids.set_edgecolor("face")
#fig.colorbar(surf, shrink=0.5, aspect=10)
rcParams['savefig.dpi'] = 300
plt.savefig('plots/squint_talk_q12_bloch_e' + str(epsilon) + '.svg', transparent=True)
