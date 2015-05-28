import matplotlib.pyplot as plt
from matplotlib import gridspec, rcParams
from matplotlib.ticker import NullLocator, NullFormatter
from mpl_toolkits.mplot3d import Axes3D, proj3d
from matplotlib.patches import FancyArrowPatch

class Arrow3D(FancyArrowPatch):
    """From http://stackoverflow.com/a/11156353/1236650

    """

    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

fontsize = 72
fig = plt.figure(figsize=(9,12))
gs = gridspec.GridSpec(2, 1, height_ratios=[3, 2])#, width_ratios=[5, 2])

ax1 = plt.subplot(gs[0], aspect='equal', alpha=0, frame_on=False)
ax1.set_xlim([-1.25, 1.25])
ax1.set_ylim([-1.25, 1.25])
ax1.xaxis.set_major_formatter(NullFormatter())
ax1.xaxis.set_major_locator(NullLocator())
ax1.yaxis.set_major_formatter(NullFormatter())
ax1.yaxis.set_major_locator(NullLocator())
ax1.arrow(-1, 0, 2, 0, width=.005, color='k')
ax1.arrow(0, -1, 0, 2, width=.005, color='k')
ax1.text(0.55, 1, r'$x$', fontsize=fontsize, horizontalalignment='left',
         verticalalignment='top', transform=ax1.transAxes)
ax1.text(1, 0.55, r'$z$', fontsize=fontsize, horizontalalignment='right',
         verticalalignment='bottom', transform=ax1.transAxes)

ax2 = plt.subplot(gs[1], projection='3d')
ax2.set_zlim3d(0, 1)
ax2.w_zaxis.set_major_formatter(NullFormatter())
ax2.w_zaxis.set_major_locator(NullLocator())
ax2.set_xlim3d(-1, 1)
ax2.w_xaxis.set_major_formatter(NullFormatter())
ax2.w_xaxis.set_major_locator(NullLocator())
ax2.set_ylim3d(-1, 1)
ax2.w_yaxis.set_major_formatter(NullFormatter())
ax2.w_yaxis.set_major_locator(NullLocator())
x_arrow = Arrow3D([0,0],[-1,1],[0,0], mutation_scale=20, lw=2, arrowstyle='-|>',
                  color='k')
y_arrow = Arrow3D([0,0],[0,0],[0,1], mutation_scale=20, lw=2, arrowstyle='-|>',
                  color='k')
z_arrow = Arrow3D([-1,1],[0,0],[0,0], mutation_scale=20, lw=2, arrowstyle='-|>',
                  color='k')
ax2.add_artist(x_arrow)
ax2.add_artist(y_arrow)
ax2.add_artist(z_arrow)

ax2.text(0, 1, 0, r'$x$', fontsize=fontsize)
ax2.text(0.05, 0, 1, r'$y$', fontsize=fontsize)
ax2.text(1, 0, 0.05, r'$z$', fontsize=fontsize)

plt.savefig('plots/xz_and_xyz_axes.svg', transparent=True)
