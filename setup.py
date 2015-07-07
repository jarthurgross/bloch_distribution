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
      # py_modules=['bloch_distribution'],
      install_requires=requires,
      packages=find_packages(),
     )
