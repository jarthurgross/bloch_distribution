#!/usr/bin/python3
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import NullLocator, NullFormatter
from mpl_toolkits.mplot3d import Axes3D, proj3d
from matplotlib.patches import FancyArrowPatch
import numpy as np
import h5py
import argparse

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

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Plot probability densities ' +
                                     'for projective measurements on the ' +
                                     'Bloch sphere in accordance with ' +
                                     'Das-Arvind.')
    parser.add_argument('data_filenames', metavar='files', nargs='+',
                        help='List of hdf5 files containing densities')
    parser.add_argument('--colormap', '-c', default='CMRmap',
                        help='String indicating colormap to use')
    parser.add_argument('--filetype', '-f', default='pdf',
                        help='Extension for output file')
    args = parser.parse_args()

    data_files = [h5py.File(data_filename, 'r') for data_filename in
                  args.data_filenames]
    epsilons = [data_file['densities'].attrs['epsilon'] for data_file in
                data_files]
    # The vstacks fill in the indeterminate polar density with the average of
    # the neighboring densities
    Density_meshes = [data_file['densities'][:] for data_file in data_files]
    Density_meshes = [np.vstack((np.mean(Density_mesh[0,:]) *
                                 np.ones(Density_mesh.shape[1]), Density_mesh))
                      for Density_mesh in Density_meshes]

    R_meshes = [data_file['R'][:] for data_file in data_files]
    R_meshes = [np.vstack((np.zeros(R_mesh.shape[1]), R_mesh)) for R_mesh in
                R_meshes]
    Phi_meshes = [data_file['Phi'][:] for data_file in data_files]
    Phi_meshes = [np.vstack((Phi_mesh[0,:], Phi_mesh)) for Phi_mesh in
                  Phi_meshes]
    X_meshes = [data_file['X'][:] for data_file in data_files]
    X_meshes = [np.vstack((np.zeros(X_mesh.shape[1]), X_mesh)) for X_mesh in
                X_meshes]
    Y_meshes = [data_file['Y'][:] for data_file in data_files]
    Y_meshes = [np.vstack((np.zeros(Y_mesh.shape[1]), Y_mesh)) for Y_mesh in
                Y_meshes]
    Z_meshes = [data_file['Z'][:] for data_file in data_files]
    Z_meshes = [np.vstack((np.ones(Z_mesh.shape[1]), Z_mesh)) for Z_mesh in
                Z_meshes]

    fontsize = 24
    mpl.rcParams['axes.labelsize'] = fontsize
    mpl.rcParams['axes.titlesize'] = fontsize

    norm = mpl.colors.Normalize()

    cmap = plt.get_cmap(args.colormap)

    fig = plt.figure(figsize=(2.5*(len(data_files) + 1),4))
    gs = mpl.gridspec.GridSpec(2, len(data_files) + 1, height_ratios=[3, 2])

    ax = plt.subplot(gs[0], aspect='equal', alpha=0, frame_on=False)
    ax.set_xlim([-1.25, 1.25])
    ax.set_ylim([-1.25, 1.25])
    ax.xaxis.set_major_formatter(NullFormatter())
    ax.xaxis.set_major_locator(NullLocator())
    ax.yaxis.set_major_formatter(NullFormatter())
    ax.yaxis.set_major_locator(NullLocator())
    ax.arrow(-1, 0, 2, 0, width=.005, color='k')
    ax.arrow(0, -1, 0, 2, width=.005, color='k')
    ax.text(0.55, 1, r'$x$', fontsize=fontsize, horizontalalignment='left',
            verticalalignment='top', transform=ax.transAxes)
    ax.text(1, 0.55, r'$z$', fontsize=fontsize, horizontalalignment='right',
            verticalalignment='bottom', transform=ax.transAxes)

    ax = plt.subplot(gs[len(data_files) + 1], projection='3d')
    ax.set_zlim3d(0, 1)
    ax.w_zaxis.set_major_formatter(NullFormatter())
    ax.w_zaxis.set_major_locator(NullLocator())
    ax.set_xlim3d(-1, 1)
    ax.w_xaxis.set_major_formatter(NullFormatter())
    ax.w_xaxis.set_major_locator(NullLocator())
    ax.set_ylim3d(-1, 1)
    ax.w_yaxis.set_major_formatter(NullFormatter())
    ax.w_yaxis.set_major_locator(NullLocator())
    x_arrow = Arrow3D([0,0],[-1,1],[0,0], mutation_scale=20, lw=2,
                      arrowstyle='-|>', color='k')
    y_arrow = Arrow3D([0,0],[0,0],[0,1], mutation_scale=20, lw=2,
                      arrowstyle='-|>', color='k')
    z_arrow = Arrow3D([-1,1],[0,0],[0,0], mutation_scale=20, lw=2,
                      arrowstyle='-|>', color='k')
    ax.add_artist(x_arrow)
    ax.add_artist(y_arrow)
    ax.add_artist(z_arrow)

    ax.text(0, 1, 0, r'$x$', fontsize=fontsize)
    ax.text(0.05, 0, 1, r'$y$', fontsize=fontsize)
    ax.text(1, 0, 0.05, r'$z$', fontsize=fontsize)

    for n in range(len(data_files)):
        ax = plt.subplot(gs[n + 1], projection='polar')
        ax.pcolormesh(Phi_meshes[n], R_meshes[n], Density_meshes[n], cmap=cmap,
                      norm=norm, shading='gouraud', rasterized=True)
        ax.set_ylim(0, np.max(R_meshes[n]))
        ax.xaxis.set_major_formatter(NullFormatter())
        ax.yaxis.set_major_formatter(NullFormatter())

        ax.set_title(r'$\epsilon=' + str(epsilons[n]) + '$')

        # Having problems with my custom colormap and facecolors.
        ax = plt.subplot(gs[len(data_files) + 1 + n + 1], projection='3d')
        surf = ax.plot_surface(X_meshes[n], Y_meshes[n], Z_meshes[n],
                               rstride=1, cstride=1,
                               facecolors=cmap(norm(Density_meshes[n])),
                               shade=False)
        surf.set_rasterized(True)
        ax.set_zlim3d(0, 1)
        ax.w_zaxis.set_major_formatter(NullFormatter())
        ax.w_zaxis.set_major_locator(NullLocator())
        ax.set_xlim3d(-1, 1)
        ax.w_xaxis.set_major_formatter(NullFormatter())
        ax.w_xaxis.set_major_locator(NullLocator())
        ax.set_ylim3d(-1, 1)
        ax.w_yaxis.set_major_formatter(NullFormatter())
        ax.w_yaxis.set_major_locator(NullLocator())

    mpl.rcParams['savefig.dpi'] = 300
    plt.savefig('../../data/plots/bloch_densities_e' + str(epsilons[0]) + '_' +
                args.colormap + '.' + args.filetype, transparent=True)
