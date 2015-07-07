#!/use/bin/python
from __future__ import division
from sampling import build_tree, get_samples
import matplotlib.pyplot as plt
from matplotlib import cm, colors
import numpy as np
from collections import Counter
from my_cms import husl_hot

# Parameters
min = -8
max = 8
linsamps = 2**8 + 1
epsilon = 0.575
distsamps = 2**24
cmap = husl_hot
fontsize = 20

tree = build_tree(min, max, linsamps, epsilon)
samples = get_samples(distsamps, tree)
frequencies = Counter(samples)
step = (max - min)/(linsamps - 1)
area = step**2
Q1 = np.linspace(min + step/2, max - step/2, linsamps - 1)
Q2 = np.linspace(min + step/2, max - step/2, linsamps - 1)
QQ1, QQ2 = np.meshgrid(Q1, Q2)
densities = np.array([[frequencies[(q1, q2)]/(area*distsamps) for q1 in Q1] for
                      q2 in Q2])

norm = colors.Normalize()

fig = plt.figure()
fig.suptitle(r'$\epsilon=' + str(epsilon) + '$', fontsize=fontsize)
ax = fig.add_subplot(1, 1, 1)
ax.pcolor(QQ1, QQ2, densities, cmap=cmap, norm=norm)
ax.set_xlim(-4, 4)
ax.set_ylim(-4, 4)
m = cm.ScalarMappable(cmap=cmap, norm=norm)
m.set_array(densities)
plt.colorbar(m)

plt.savefig('plots/q12_dist_sampled_e' + str(epsilon) + '.png')
