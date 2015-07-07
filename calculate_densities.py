from __future__ import division, print_function
import numpy as np
from numpy import pi, sin, cos
import h5py
from datetime import datetime as dt
import argparse
from invert_angles import G_angles_q12, construct_grid

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Calculate probability ' +
                                     'densities on the Bloch sphere for ' +
                                     'given eigenbasis measurements ' +
                                     'in equivalence with Das-Arvind.')
    parser.add_argument('-e', '--epsilon', type=float, default=0.575,
                        help='Strength of the weak measurements')
    parser.add_argument('-N', type=int, default=256,
                        help='Resolution at which to calculate densities')
    parser.add_argument('-f', '--folder', type=str, default='../../data/',
                        help='Folder to write results to')
    args = parser.parse_args()

    data_file = h5py.File(args.folder + 'bloch_densities_e' +
                          str(args.epsilon) + '_N' + str(args.N) + '_' +
                          dt.utcnow().strftime('%Y%m%dT%H%M%SZ') + '.hdf5',
                          'w-')

    Theta, Phi = construct_grid(args.N)

    # Radial coordinates for a Lambert azimuthal equal-area projection
    R = 2*np.cos(pi/2 - Theta/2)

    X = np.sin(Theta) * np.cos(Phi)
    Y = np.sin(Theta) * np.sin(Phi)
    Z = np.cos(Theta)

    Angles = np.array([cos(Theta), Phi])

    Densities = G_angles_q12(Angles, args.epsilon)

    dataset = data_file.create_dataset('densities', data=Densities)
    dataset.attrs['epsilon'] = args.epsilon
    dataset.attrs['N'] = args.N
    data_file.create_dataset('R', data=R)
    data_file.create_dataset('Theta', data=Theta)
    data_file.create_dataset('Phi', data=Phi)
    data_file.create_dataset('X', data=X)
    data_file.create_dataset('Y', data=Y)
    data_file.create_dataset('Z', data=Z)

    data_file.close()

