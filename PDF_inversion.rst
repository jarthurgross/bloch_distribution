.. Explanation of the PDF inversion process

PDF inversion
=============

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
   \phi&=\operatorname{atan2}(2\sinh(\epsilon q_2),\sinh(\epsilon(q_1+q_2))+
   \sinh(\epsilon(q_1-q_2)))
   \end{align}

where :math:`(\theta,\phi)` are spherical coördinates about the y-axis:

.. math::

   \begin{align}
   z&=\sin\theta\cos\phi
   x&=\sin\theta\sin\phi
   y&=\cos\theta
   \end{align}
