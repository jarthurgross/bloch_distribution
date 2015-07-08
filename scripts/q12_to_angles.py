#!/usr/bin/python3
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FixedFormatter
from matplotlib import cm, colors, gridspec
import numpy as np
from numpy import pi
from bloch_distribution.invert_angles import map_q12_to_sphere
from my_cms import huslp

# Parameters
epsilon = .05
q1_min = -4
q1_max = 4
q1_samples = 256
q2_min = -4
q2_max = 4
q2_samples = 256

Q1 = np.linspace(q1_min, q1_max, q1_samples)
Q2 = np.linspace(q2_min, q2_max, q2_samples)
Q1, Q2 = np.meshgrid(Q1, Q2)

CosTheta, Phi = map_q12_to_sphere(np.array([Q1, Q2]), epsilon)

norm = colors.Normalize()

cmap = huslp

fontsize = 20
fig = plt.figure(figsize=(18,6))
fig.suptitle(r'$\epsilon=' + str(epsilon) + '$', fontsize=fontsize)
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 2])

ax0 = plt.subplot(gs[0])
ax0.pcolormesh(Q1, Q2, Phi, cmap=cmap, norm=norm, shading='gouraud')
ax0.contour(Q1, Q2, Phi, colors='w', linewidths=2,
            levels=np.linspace(-pi, pi, 9))

ax0.set_xlabel(r'$q_1$', fontsize=fontsize)
ax0.set_ylabel(r'$q_2$', fontsize=fontsize)

ax1 = plt.subplot(gs[1], projection='3d')
ax1.plot_surface(Q1, Q2, CosTheta, rstride=1, cstride=1,
                 facecolors=cmap(norm(Phi)))
ax1.set_zlim3d(0, 1)
ax1.set_xlabel(r'$q_1$', fontsize=fontsize)
ax1.set_ylabel(r'$q_2$', fontsize=fontsize)
ax1.set_zlabel(r'$\cos\theta$', fontsize=fontsize)

m = cm.ScalarMappable(cmap=cmap, norm=norm)
m.set_array(Phi)
plt.colorbar(m, ticks=MultipleLocator(pi/4),
             format=FixedFormatter([r'$-\pi$', r'$-3\pi/4$', r'$-\pi/2$',
                                    r'$-\pi/4$', r'$0$', r'$\pi/4$',
                                    r'$\pi/2$', r'$3\pi/4$', r'$\pi$']))
plt.savefig('plots/q12_to_angles_e' + str(epsilon) + '.png')
