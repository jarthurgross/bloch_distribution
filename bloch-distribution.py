#!/usr/bin/python3
import matplotlib.pyplot as plt
from matplotlib import cm, colors, gridspec
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from invert_angles import *

cmap = cm.hot

fontsize = 20
fig = plt.figure(figsize=(24,8))
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 2])
fig.suptitle(r'$\epsilon=' + str(epsilon) + '$', fontsize=fontsize)

ax1 = plt.subplot(gs[0], projection='polar')
ax1.pcolormesh(Phi, Theta, colorfunction, cmap=cmap, norm=norm,
               shading='gouraud')

ax2 = plt.subplot(gs[1], projection='3d')
surf = ax2.plot_surface(X, Y, Z,  rstride=1, cstride=1,
                facecolors=cmap(norm(colorfunction)), shade=False)
ax2.set_zlim3d(0, 1)
ax2.set_xlim3d(-1, 1)
ax2.set_ylim3d(-1, 1)
ax2.set_xlabel(r'$z$', fontsize=fontsize)
ax2.set_ylabel(r'$x$', fontsize=fontsize)
ax2.set_zlabel(r'$y$', fontsize=fontsize)

m = cm.ScalarMappable(cmap=cmap, norm=norm)
m.set_array(colorfunction)
plt.colorbar(m)
#fig.colorbar(surf, shrink=0.5, aspect=10)
plt.savefig('plots/e' + str(epsilon) + '.png')
