from setuptools import setup

requires = [
        'numpy',
        'scipy',
        'matplotlib',
        'BTrees',
        'h5py',
        'Cython',
        ]

setup(name='bloch_distribution',
      version='0.0',
      py_modules=['invert_angles', 'sampling'],
      install_requires=requires,
      # Workaround from
      # https://github.com/numpy/numpy/issues/2434#issuecomment-65252402
      # and
      # https://github.com/h5py/h5py/issues/535#issuecomment-79158166
      setup_requires=['numpy', 'Cython'],
      packages=['bloch_distribution'],
      package_dir={'bloch_distribution': 'src/bloch_distribution'},
     )
