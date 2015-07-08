from setuptools import setup, find_packages

requires = [
        'numpy',
        'scipy',
        'matplotlib',
        'BTrees',
        'h5py',
        ]

setup(name='bloch_distribution',
      version='0.0',
      py_modules=['invert_angles', 'sampling'],
      install_requires=requires,
      packages=['bloch_distribution'],
      package_dir={'bloch_distribution': 'src/bloch_distribution'},
     )
