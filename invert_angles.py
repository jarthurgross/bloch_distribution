#!/usr/bin/python3
from __future__ import division
from scipy.optimize import root
import numpy as np
from numpy import cosh, arccosh, sinh, arctan2, cos, sin, exp, pi
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D

atan2 = arctan2


def map_qpm_to_sphere_y(qs, epsilon):
    """Takes values for qp and qm and maps them to points on the surface of the
    Bloch sphere given as pairs (cos theta, phi) for spherical coordinates
    around the y-axis with phi starting on the z-axis.

    """
    cos_theta = 2/(cosh(2*epsilon*qs[0]) + cosh(2*epsilon*qs[1]))
    phi = atan2(sinh(epsilon*(qs[0] - qs[1])),
                (sinh(2*epsilon*qs[0]) + sinh(2*epsilon*qs[1]))/2)
    return np.array([cos_theta, phi])


def map_q12_to_sphere_y(qs, epsilon):
    """Takes values for q1 and q2 and maps them to points on the surface of the
    Bloch sphere given as pairs (cos theta, phi) for spherical coordinates
    around the y-axis with phi starting on the z-axis.

    """
    cos_theta = 2/(cosh(epsilon*(qs[0] + qs[1])) + 
                   cosh(epsilon*(qs[0] - qs[1])))
    phi = atan2(sinh(epsilon*qs[1]),
                (sinh(epsilon*(qs[0] + qs[1])) +
                 sinh(epsilon*(qs[0] - qs[1])))/2)
    return np.array([cos_theta, phi])


def map_q12_to_sphere_z(qs, epsilon):
    """Takes values for q1 and q2 and maps them to points on the surface of the
    Bloch sphere given as pairs (cos theta, phi) for spherical coordinates
    around the z-axis with phi starting on the x-axis (i.e. the standard
    arrangement).

    """
    cos_theta = ((sinh(epsilon*(qs[0] + qs[1])) +
                  sinh(epsilon*(qs[0] - qs[1])))/
                 (cosh(epsilon*(qs[0] + qs[1])) + 
                  cosh(epsilon*(qs[0] - qs[1]))))
    phi = atan2(2, sinh(epsilon*qs[1]))
    return np.array([cos_theta, phi])


# angles = np.array([cos(Theta), Phi])
def guess_qpm_vals(angles, epsilon):
    """Generates an initial guess of qp and qm for the rootfinder.

    """
    return np.where(angles[0] == 1, np.zeros(angles.shape),
                    np.array([
                    cos(pi/4 - angles[1])*arccosh(2/angles[0] - 1)/(2*epsilon),
                    sin(pi/4 - angles[1])*arccosh(2/angles[0] - 1)/(2*epsilon)
                    ]))


# angles = np.array([cos(Theta), Phi])
def array_get_qpm_vals(angles, epsilon):
    """Takes points on the upper hemisphere of the bloch sphere and maps them to
    points on the qp qm plane.

    """
    costheta = np.array(angles[0])
    phi = np.array(angles[1])
    # Put phi in the range [-pi/2, 3pi/2).
    phi = (phi + pi/2)%(2*pi) - pi/2
    # If phi in range [pi/2, 3pi/2), solve for qs using phi - pi and negate the
    # q values obtained.
    sign = np.where(phi > pi/2, -1, 1)
    phi = np.where(phi > pi/2, phi - pi, phi)
    # Currently solving for values individually, but when passed an array of
    # of angle values it won't converge (I think the rootfinder thinks it has
    # to solve for all the roots simultaneously as a large system of
    # equations).
    inverted = np.zeros(np.array(angles).shape)
    for m in range(angles.shape[1]):
        for n in range(angles.shape[2]):
            if costheta[m,n] == 1:
                inverted[:,m,n] = np.array([0, 0])
            else:
                try:
                    sol = root(lambda qs: map_qpm_to_sphere_y(qs, epsilon) -
                                np.array([costheta[m,n], phi[m,n]]),
                                guess_qpm_vals(np.array([costheta[m,n],
                                phi[m,n]]), epsilon))
                    inverted[:,m,n] = sign[m,n]*sol.x
                except Exception as e:
                    print('cos(theta) = ' + str(costheta[m,n]) +
                          ', phi = ' + str(phi[m,n]) + ', sign = ' +
                          str(sign[m,n]) + '\n' + str(e))

    return inverted
                

def G_qpm(qp, qm, epsilon):
    """The probability density function on the qp qm plane. (might be improperly
    normalized)

    """
    return (epsilon/(4*pi)*exp(-epsilon*(qp**2 + qm**2 + 1))*
            (cosh(2*epsilon*qp) + cosh(2*epsilon*qm)))


def G_q12(q1, q2, epsilon):
    """The probability density function on the q1 q2 plane. (might be improperly
    normalized)

    """
    return (epsilon/(4*pi)*exp(-epsilon*(q1**2 + q2**2 + 2)/2)*(cosh(epsilon*(q1
        + q2)) + cosh(epsilon*(q1 - q2))))


def parallelogram_area_qpm(qp, qm, epsilon):
    """Calculate the area of dOmega in units of dqp*dqm.

    """
    s_pp = sinh(2*epsilon*qp)
    s_pm = sinh(epsilon*(qp + qm))
    s_mp = sinh(epsilon*(qp - qm))
    s_mm = sinh(2*epsilon*qm)
    c_pp = cosh(2*epsilon*qp)
    c_mp = cosh(epsilon*(qp - qm))
    c_mm = cosh(2*epsilon*qm)
    return (np.abs((s_pp + s_mm)**2*c_mp + 2*(s_pp - s_mm)*s_mp*c_pp) /
            ((c_pp + c_mm)**2*(4*s_pm**2 + s_mp**2 + s_pp**2 + s_mm**2)))


def parallelogram_area_q12(q1, q2, epsilon):
    """Calculate the area of dOmega in units of dq1*dq2

    """
    cp = cosh(epsilon*(q1 + q2))
    cm = cosh(epsilon*(q1 - q2))
    c2 = cosh(epsilon*q2)
    sp = sinh(epsilon*(q1 + q2))
    sm = sinh(epsilon*(q1 - q2))
    s2 = sinh(epsilon*q2)

    dcostheta = (2*epsilon/(cp + cm)**2)*np.array([sp + sm, sp - sm])
    dphi = ((2*epsilon/((sp + sm)**2 + 4*s2**2)) *
            np.array([-s2*(cp + cm), c2*(sp + sm) - s2*(cp - cm)]))
    return np.abs(dcostheta[0]*dphi[1] - dcostheta[1]*dphi[0])


def G_angles(angles, epsilon):
    """The probability density function on the upper hemisphere of the bloch
    sphere.

    """
    qp, qm = array_get_qpm_vals(angles, epsilon)
    return G_qpm(qp, qm, epsilon)/parallelogram_area(qp, qm, epsilon)
