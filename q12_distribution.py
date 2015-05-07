#!/usr/bin/python3
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from matplotlib.ticker import NullLocator, NullFormatter
import numpy as np
from invert_angles import G_q12
from my_cms import husl_hot

# Parameters
epsilon = 0.05

q1_min = -4
q1_max = 4
q1_samples = 512
q2_min = -4
q2_max = 4
q2_samples = 512

Q1 = np.linspace(q1_min, q1_max, q1_samples)
Q2 = np.linspace(q2_min, q2_max, q2_samples)
Q1, Q2 = np.meshgrid(Q1, Q2)

# Set the bounds so that epsilon*q does not exceed a certain value
eq_min = -3
eq_max = 3
eq_samples = 512

EQ1, EQ2 = np.mgrid[eq_min:eq_max:eq_samples*1j, eq_min:eq_max:eq_samples*1j]

G = G_q12(EQ1/epsilon, EQ2/epsilon, epsilon)/(epsilon**2)

norm = colors.Normalize()

cmap = husl_hot

fontsize = 36
fig = plt.figure()
fig.suptitle(r'$\epsilon=' + str(epsilon) + '$', fontsize=fontsize)

ax = fig.add_subplot(1, 1, 1, aspect='equal')
ax.pcolormesh(EQ1, EQ2, G, cmap=cmap, norm=norm, shading='gouraud',
              rasterized=True)

ax.set_xlabel(r'$\epsilon q_1$', fontsize=fontsize)
ax.set_ylabel(r'$\epsilon q_2$', fontsize=fontsize)
ax.xaxis.set_major_formatter(NullFormatter())
ax.yaxis.set_major_formatter(NullFormatter())
ax.xaxis.set_major_locator(NullLocator())
ax.yaxis.set_major_locator(NullLocator())

'''
m = cm.ScalarMappable(cmap=cmap, norm=norm)
m.set_array(G)
plt.colorbar(m, ticks=NullLocator(), format=NullFormatter())
'''
plt.savefig('plots/eq_bound_3_husl_hot_q12_e' + str(epsilon) + '.svg')
