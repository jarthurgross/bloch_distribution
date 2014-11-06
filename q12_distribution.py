#!/usr/bin/python3
import matplotlib.pyplot as plt
from matplotlib import cm, colors
import numpy as np
from invert_angles import G_q12
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

G = G_q12(Q1, Q2, epsilon)

norm = colors.Normalize()

cmap = husl_hot

fontsize = 20
fig = plt.figure()
fig.suptitle(r'$\epsilon=' + str(epsilon) + '$', fontsize=fontsize)

ax = plt.subplot(1, 1, 1)
ax.pcolormesh(Q1, Q2, G, cmap=cmap, norm=norm, shading='gouraud')

ax.set_xlabel(r'$q_1$', fontsize=fontsize)
ax.set_ylabel(r'$q_2$', fontsize=fontsize)

m = cm.ScalarMappable(cmap=cmap, norm=norm)
m.set_array(G)
plt.colorbar(m)
plt.savefig('plots/husl_hot_q12_e' + str(epsilon) + '.png')
