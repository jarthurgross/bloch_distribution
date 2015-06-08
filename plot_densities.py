#!/usr/bin/python3
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import NullLocator, NullFormatter
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import h5py
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Plot probability densities ' +
                                     'for projective measurements on the ' +
                                     'Bloch sphere in accordance with ' +
                                     'Das-Arvind.')
    parser.add_argument('data_filenames', metavar='files', nargs='+',
                        help='List of hdf5 files containing densities')
    parser.add_argument('--colormap', '-c', default='CMRmap',
                        help='List of hdf5 files containing densities')
    args = parser.parse_args()

    data_file = h5py.File(args.data_filenames[0], 'r')
    epsilon = data_file['densities'].attrs['epsilon']
    # The vstacks fill in the guessed limit of epsilon^3*exp(-epsilon/2)/2*pi
    # for theta = 0
    Densities = data_file['densities'][:]
    Densities = np.vstack((epsilon**3 * np.exp(-epsilon/2) *
                           np.ones(Densities.shape[1])/(2 * np.pi), Densities))
    R = data_file['R'][:]
    R = np.vstack((np.zeros(R.shape[1]), R))
    Phi = data_file['Phi'][:]
    Phi = np.vstack((Phi[0,:], Phi))
    X = data_file['X'][:]
    X = np.vstack((np.zeros(X.shape[1]), X))
    Y = data_file['Y'][:]
    Y = np.vstack((np.zeros(Y.shape[1]), Y))
    Z = data_file['Z'][:]
    Z = np.vstack((np.ones(Z.shape[1]), Z))

    mpl.rcParams['axes.labelsize'] = 72
    mpl.rcParams['axes.titlesize'] = 72

    norm = mpl.colors.Normalize()

    cmap = plt.get_cmap(args.colormap)

    fig = plt.figure(figsize=(9,12))
    gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[3, 2])
    fig.suptitle(r'$\epsilon=' + str(epsilon) + '$')

    ax1 = plt.subplot(gs[0], projection='polar')
    ax1.pcolormesh(Phi, R, Densities, cmap=cmap, norm=norm,
                   shading='gouraud', rasterized=True)
    ax1.set_ylim(0, np.max(R))
    ax1.xaxis.set_major_formatter(NullFormatter())
    ax1.yaxis.set_major_formatter(NullFormatter())

    # Having problems with my custom colormap and facecolors.
    ax2 = plt.subplot(gs[1], projection='3d')
    surf = ax2.plot_surface(X, Y, Z,  rstride=1, cstride=1,
                            facecolors=cmap(norm(Densities)), shade=False)
    surf.set_rasterized(True)
    ax2.set_zlim3d(0, 1)
    ax2.w_zaxis.set_major_formatter(NullFormatter())
    ax2.w_zaxis.set_major_locator(NullLocator())
    ax2.set_xlim3d(-1, 1)
    ax2.w_xaxis.set_major_formatter(NullFormatter())
    ax2.w_xaxis.set_major_locator(NullLocator())
    ax2.set_ylim3d(-1, 1)
    ax2.w_yaxis.set_major_formatter(NullFormatter())
    ax2.w_yaxis.set_major_locator(NullLocator())

    mpl.rcParams['savefig.dpi'] = 300
    plt.savefig('../../data/plots/bloch_densities_e' + str(epsilon) + '_' +
                args.colormap + '.pdf', transparent=True)
