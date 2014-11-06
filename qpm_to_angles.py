#!/usr/bin/python3
import matplotlib.pyplot as plt
from matplotlib import cm, colors
import numpy as np
from numpy import pi
from invert_angles import map_qpm_to_sphere
from my_cms import huslp

# Parameters
epsilon = 0.05
qp_min = -4
qp_max = 4
qp_samples = 512
qm_min = -4
qm_max = 4
qm_samples = 512

Qp = np.linspace(qp_min, qp_max, qp_samples)
Qm = np.linspace(qm_min, qm_max, qm_samples)
Qp, Qm = np.meshgrid(Qp, Qm)

CosTheta, Phi = map_qpm_to_sphere(np.array([Qp, Qm]), epsilon)

norm = colors.Normalize()

cmap = huslp

fontsize = 20
fig = plt.figure()
fig.suptitle(r'$\epsilon=' + str(epsilon) + '$', fontsize=fontsize)

ax = plt.subplot(1, 1, 1)
ax.pcolormesh(Qp, Qm, Phi, cmap=cmap, norm=norm, shading='gouraud')
ax.contour(Qm, Qp, Phi, colors='w', linewidths=2,
           levels=np.linspace(-pi, pi, 9))

ax.set_xlabel(r'$q_+$', fontsize=fontsize)
ax.set_ylabel(r'$q_-$', fontsize=fontsize)

m = cm.ScalarMappable(cmap=cmap, norm=norm)
m.set_array(Phi)
plt.colorbar(m)
plt.savefig('plots/qpm_to_angles_e' + str(epsilon) + '.png')
