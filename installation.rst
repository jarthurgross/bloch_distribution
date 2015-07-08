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

.. _Anaconda: http://continuum.io/downloads
