#!/usr/bin/python3
import matplotlib.pyplot as plt
from matplotlib import cm, colors
import numpy as np
from invert_angles import G_qpm

# Parameters
epsilon = 0.575
qp_min = -4
qp_max = 4
qp_samples = 512
qm_min = -4
qm_max = 4
qn_samples = 512

Qp = np.linspace(qp_min, qp_max, qp_samples)
Qm = np.linspace(qp_min, qp_max, qp_samples)
Qp, Qm = np.meshgrid(Qp, Qm)

G = G_qpm(Qp, Qm, epsilon)

norm = colors.Normalize()

cmap = cm.hot

fontsize = 20
fig = plt.figure()
fig.suptitle(r'$\epsilon=' + str(epsilon) + '$', fontsize=fontsize)

ax = plt.subplot(1, 1, 1)
ax.pcolormesh(Qp, Qm, G, cmap=cmap, norm=norm, shading='gouraud')

ax.set_xlabel(r'$q_+$', fontsize=fontsize)
ax.set_ylabel(r'$q_-$', fontsize=fontsize)

m = cm.ScalarMappable(cmap=cmap, norm=norm)
m.set_array(G)
plt.colorbar(m)
#fig.colorbar(surf, shrink=0.5, aspect=10)
plt.savefig('plots/qpm_e' + str(epsilon) + '.png')
