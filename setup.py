# In case the user doesn't have setuptools installed
import ez_setup
ez_setup.use_setuptools()

from setuptools import setup

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
      install_requires=requires
     )
