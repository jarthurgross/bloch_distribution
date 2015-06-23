import argparse
import numpy as np
import h5py

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Integrate probability ' +
                                     'densities to verify that they are ' +
                                     'normalized')
    parser.add_argument('data_filenames', metavar='files', nargs='+',
                        help='List of hdf5 files containing densities')
    args = parser.parse_args()

    data_files = [h5py.File(data_filename, 'r') for data_filename in
                  args.data_filenames]
    epsilons = [data_file['densities'].attrs['epsilon'] for data_file in
                data_files]
    Density_meshes = [data_file['densities'][:] for data_file in data_files]
    Phi_meshes = [data_file['Phi'][:] for data_file in data_files]
    Theta_meshes = [-2*np.arccos(data_file['R'][:]/2) + np.pi for data_file in
                    data_files]

    Total_probs = []

    for Density_mesh, Phi_mesh, Theta_mesh in zip(Density_meshes, Phi_meshes,
                                                  Theta_meshes):
        # Scale Density_mesh so that the integration can be thought of as on a
        # rectangle rather than a hemisphere
        Scaled_density_mesh = Density_mesh*np.sin(Theta_mesh)
        Total_probs.append(np.trapz(np.trapz(Scaled_density_mesh, Phi_mesh),
                                    Theta_mesh[:,0]))

    for epsilon, prob in zip(epsilons, Total_probs):
        print(epsilon, prob)
