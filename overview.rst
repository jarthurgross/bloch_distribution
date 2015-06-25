.. Overview of the module

Overview
========

:code:`invert_angles.py`
------------------------

:code:`invert_angles.py` is the meat of this module. Everything else is a
collection of scripts making use of this module in some way. This module can map
points from the hemisphere to the plane and vice-versa, calculate the Jacobian
of this mapping, and calculate probability densities on both manifolds.

:code:`calculate_densities.py`
------------------------------

This script calculates probability density distributions on the Bloch sphere for
given measurement strengths and saves the results to an :code:`hdf5` file.

:code:`plot_densities.py`
-------------------------

This script plots the probability density distribution stored in an :code:`hdf5`
file. Allows you to fiddle with the appearance of the graphic without having to
recalculate densities each time.

:code:`samling.py`
------------------

This module allows the user to sample from the probability distributions on the
plane. This is useful for performing sequential Monte Carlo simulations.
Currently the sampler has better performance if you get all the samples you need
with one function call rather than splitting the sampling across function calls,
since some data structures have to be set up with each function call.

:code:`q12` versus :code:`qpm`
------------------------------

There are several files and functions in this module that have identical names
save for the difference between :code:`q12` and :code:`qpm`. The :code:`qpm`
files and functions were written using a parametrization of the plane
:math:`q_\pm=(q_1\pm q_2)/2`. These functions gave incorrect results that are
explored in the scripts that plot parallelogram areas (that is, the determinates
of the Jacobian at different points). Since in my analytic work I was moving
away from that parameterization back to simply using :math:`q_1` and
:math:`q_2`, and I could not easily determine where I had made my mistake in the
functions, I determined to rewrite everything using the :math:`(q_1,q_2)`
parameterization. As far as I can tell, this eliminated whatever bug I had
introduced.

Ideally, one should be able to go back and fix the :code:`qpm` code to agree
with the :code:`q12` code, but I haven't bothered to do that.

**tl;dr**: Don't use :code:`qpm` scripts or functions. Use :code:`q12` scripts
and functions instead.

:code:`plot_parallelogram_area_q*.py`
-------------------------------------

Plots the value of the Jacobian determinant at different points on the plane.
Useful for debugging purposes.

:code:`q*_to_angles.py`
-----------------------

Visualization code designed to help me understand how to invert the map from the
plane to the hemisphere.

:code:`plot_q12_dist_samples.py`
--------------------------------

Plots a histogram of samples from a distribution on the plane for a quick visual
check to make sure the sampler is drawing from the distribution it's supposed
to.

:code:`tests.py`
----------------

A few test cases to increase confidence that the code is working as it should.
