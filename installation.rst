.. Installation instructions

Installation
============

Linux
-----

Installation is easiest if you already have a working `Anaconda`_ installation.
In this case, simply navigate to the root directory of the project where the
``setup.py`` file is and execute ``python setup.py install``.

If you don't have Anaconda, then you may need to also install the HDF5 shared
library with development headers (1.8.4 or newer, packaged as ``libhdf5-dev`` or
similar). Once you have that installed you should be able to successfully run
``python setup.py install``.

Issues
------

Building seems to work fine with Python 3.4 and 2.7 when starting with a base
Anaconda installation. When building with a minimal Python installation (i.e. in
an environment created by ``conda create --name new_env python=2.7`` or
``conda create --name new_env python=3``, there seem to be problems building
matplotlib. In Python 3 I get ``TypeError: unorderable types: str() < int()`` at
the line reading ``if self.version < other.version:``, and in Python 2 I get::

  * The following required packages can not be built:
  * freetype

.. _Anaconda: http://continuum.io/downloads
