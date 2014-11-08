#!/usr/bin/python3
import matplotlib.pyplot as plt
from matplotlib import cm, colors, gridspec
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from numpy import pi, sin, cos
from invert_angles import G_angles_qpm

# Parameters
N = 16
epsilon = 1

Theta, Phi = np.mgrid[0:pi/2:complex(0, N), 0:2*pi:complex(0, 2*N + 1)]
# Theta should be in the open interval (0, pi/2)
Theta = Theta[1:-1,:]
Phi = Phi[1:-1,:]

X = sin(Theta) * cos(Phi)
Y = sin(Theta) * sin(Phi)
Z = cos(Theta)

angles = np.array([cos(Theta), Phi])

colorfunction=G_angles_qpm(angles, epsilon)

norm = colors.Normalize()

cmap = cm.hot

fontsize = 20
fig = plt.figure(figsize=(24,8))
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 2])
fig.suptitle(r'$\epsilon=' + str(epsilon) + '$', fontsize=fontsize)

ax1 = plt.subplot(gs[0], projection='polar')
ax1.pcolormesh(Phi, Theta, colorfunction, cmap=cmap, norm=norm,
               shading='gouraud')
ax1.set_ylim(0, Theta[-1,-1])

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
