#!/usr/bin/python3
import matplotlib.pyplot as plt
from matplotlib import cm, colors
import numpy as np
from invert_angles import parallelogram_area_q12
from my_cms import husl_hot

# Parameters
epsilon = 0.575
q1_min = -4
q1_max = 4
q1_samples = 512
q2_min = -4
q2_max = 4
q2_samples = 512

Q1 = np.linspace(q1_min, q1_max, q1_samples)
Q2 = np.linspace(q2_min, q2_max, q2_samples)
Q1, Q2 = np.meshgrid(Q1, Q2)

Area = parallelogram_area_q12(Q1, Q2, epsilon)
Area_symmetry = Area - parallelogram_area_q12(Q1, -Q2, epsilon)

norm0 = colors.Normalize()
norm1 = colors.Normalize()

cmap0 = husl_hot
cmap1 = cm.coolwarm

fontsize = 20
fig = plt.figure(figsize=(16, 8))
fig.suptitle(r'$\epsilon=' + str(epsilon) + '$', fontsize=fontsize)

ax0 = plt.subplot(1, 2, 1)
ax0.pcolormesh(Q1, Q2, Area, cmap=cmap0, norm=norm0, shading='gouraud')
ax0.plot([-4, 4], [0, 0], color='w')

ax0.set_xlabel(r'$q_1$', fontsize=fontsize)
ax0.set_ylabel(r'$q_2$', fontsize=fontsize)
m0 = cm.ScalarMappable(cmap=cmap0, norm=norm0)
m0.set_array(Area)
plt.colorbar(m0)

ax1 = plt.subplot(1, 2, 2)
ax1.pcolormesh(Q1, Q2, Area_symmetry, cmap=cmap1, norm=norm1,
               shading='gouraud')

ax1.set_xlabel(r'$q_1$', fontsize=fontsize)
ax1.set_ylabel(r'$q_2$', fontsize=fontsize)

m1 = cm.ScalarMappable(cmap=cmap1, norm=norm1)
m1.set_array(Area_symmetry)
plt.colorbar(m1)
plt.savefig('plots/parallelogram_area_q12_e' + str(epsilon) + '.png')
