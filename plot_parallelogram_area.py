#!/usr/bin/python3
import matplotlib.pyplot as plt
from matplotlib import cm, colors
import numpy as np
from invert_angles import parallelogram_area_qpm
from my_cms import husl_hot

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

Area = parallelogram_area_qpm(Qp, Qm, epsilon)
Area_symmetry = Area - parallelogram_area_qpm(Qm, Qp, epsilon)

norm0 = colors.Normalize()
norm1 = colors.Normalize()

cmap0 = husl_hot
cmap1 = cm.coolwarm

fontsize = 20
fig = plt.figure(figsize=(16, 8))
fig.suptitle(r'$\epsilon=' + str(epsilon) + '$', fontsize=fontsize)

ax0 = plt.subplot(1, 2, 1)
ax0.pcolormesh(Qp, Qm, Area, cmap=cmap0, norm=norm0, shading='gouraud')
ax0.plot([-4, 4], [-4, 4], color='w')

ax0.set_xlabel(r'$q_p$', fontsize=fontsize)
ax0.set_ylabel(r'$q_m$', fontsize=fontsize)
m0 = cm.ScalarMappable(cmap=cmap0, norm=norm0)
m0.set_array(Area)
plt.colorbar(m0)

ax1 = plt.subplot(1, 2, 2)
ax1.pcolormesh(Qp, Qm, Area_symmetry, cmap=cmap1, norm=norm1,
               shading='gouraud')

ax1.set_xlabel(r'$q_p$', fontsize=fontsize)
ax1.set_ylabel(r'$q_m$', fontsize=fontsize)

m1 = cm.ScalarMappable(cmap=cmap1, norm=norm1)
m1.set_array(Area_symmetry)
plt.colorbar(m1)
plt.savefig('plots/parallelogram_area_qpm_e' + str(epsilon) + '.png')
