.. Explanation of the PDF inversion process

PDF change of variables
=======================

The goal of this module is to numerically evaluate the probability density at a
point on the Bloch sphere. The probability density is known analytically on
:math:`\mathbb{R}^2` for pairs of real variables :math:`(q_1,q_2)` as:

.. math::

   \begin{align}
   G(q_1,q_2)&=\frac{\epsilon}{4\pi}\exp\left\{-\frac{\epsilon(q_1^2+q_2^2+1)}
   {2}\right\}(\cosh(\epsilon(q_1+q_2))+\cosh(\epsilon(q_1-q_2)))
   \end{align}

We also have an analytic expression for mapping points in :math:`\mathbb{R}^2`
to points on the Bloch sphere:

.. math::

   \begin{align}
   \cos\theta&=\frac{2}{\cosh(\epsilon(q_1+q_2))+\cosh(\epsilon(q_1-q_2))} \\
   \phi&=\operatorname{atan2}\Big(2\sinh(\epsilon q_2),\,
   \sinh(\epsilon(q_1+q_2))+\sinh(\epsilon(q_1-q_2))\Big)
   \end{align}

where :math:`(\theta,\phi)` are spherical coördinates about the y-axis:

.. math::

   \begin{align}
   z&=\sin\theta\cos\phi \\
   x&=\sin\theta\sin\phi \\
   y&=\cos\theta
   \end{align}

Let :math:`\vec{\omega}:=(\cos\theta,\phi)` and :math:`\vec{q}:=(q_1,q_2)`, and
define :math:`\Omega:\vec{q}\mapsto\vec{\omega}`. The differentials in these two
sets of variables transform as:

.. math::

   \begin{align}
   d(\cos\theta)d\phi&=\vert\operatorname{det}(\mathrm{D}\Omega)\vert dq_1dq_2
   \end{align}

where :math:`\mathrm{D}\Omega` is the Jacobian:

.. math::

   \begin{align}
   \begin{bmatrix}
   \frac{d(\cos\theta)}{dq_1} & \frac{d(\cos\theta)}{dq_2} \\
   \frac{d\phi}{dq_1}         & \frac{d\phi}{dq_2}
   \end{bmatrix}
   \end{align}

The probability density function we want is :math:`\tilde{G}` such that:

.. math::

   \begin{align}
   \tilde{G}(\cos\theta,\phi)d(\cos\theta)d\phi&=
   G(\Omega^{-1}(\cos\theta,\phi))dq_1dq_2
   \end{align}

Using the change of variables formula we can write:

.. math::

   \begin{align}
   \tilde{G}(\cos\theta,\phi)&=\frac{G(\Omega^{-1}(\cos\theta,\phi))}
   {\vert\operatorname{det}(\mathrm{D}\Omega)\vert}
   \end{align}

Jacobian
--------

We calculate the matrix elements of the Jacobian to be:

.. math::

   \begin{align}
   \frac{d(\cos\theta)}{dq_1}&=-2\epsilon\frac{\sinh(\epsilon q_+)+
   \sinh(\epsilon q_-)}{\big(\cosh(\epsilon q_+)+\cosh(\epsilon q_-)\big)^2} \\
   \frac{d(\cos\theta)}{dq_2}&=-2\epsilon\frac{\sinh(\epsilon q_+)-
   \sinh(\epsilon q_-)}{\big(\cosh(\epsilon q_+)+\cosh(\epsilon q_-)\big)^2} \\
   \frac{d\phi}{dq_1}&=-2\epsilon\frac{\sinh(\epsilon q_2)\big(
   \cosh(\epsilon q_+)+\cosh(\epsilon q_-)\big)}{\big(\sinh(\epsilon q_+)+
   \sinh(\epsilon q_-)\big)^2+4\sinh^2(\epsilon q_2)} \\
   \frac{d\phi}{dq_2}&=-2\epsilon\frac{\sinh(\epsilon q_2)\big(
   \cosh(\epsilon q_+)-\cosh(\epsilon q_-)\big)-\cosh(\epsilon q_2)\big(
   \sinh(\epsilon q_+)+\sinh(\epsilon q_-)\big)}{\big(\sinh(\epsilon q_+)+
   \sinh(\epsilon q_-)\big)^2+4\sinh^2(\epsilon q_2)}
   \end{align}

using the shorthand :math:`q_\pm:=q_1\pm q_2`.
