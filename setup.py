from setuptools import setup

requires = [
        'numpy',
        'scipy',
        'matplotlib',
        'BTrees',
        ]

setup(name='bloch_distribution',
      version='0.0',
      py_modules=['invert_angles', 'sampling'],
      install_requires=requires
     )
