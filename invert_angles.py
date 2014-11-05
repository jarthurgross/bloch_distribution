#!/usr/bin/python3
from scipy.optimize import root
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D

cosh = np.cosh
arccosh = np.arccosh
sinh = np.sinh
atan2 = np.arctan2
cos = np.cos
sin = np.sin
exp = np.exp
pi = np.pi

def map_qpm_to_sphere(qs, epsilon):
    """Takes values for qp and qm and maps them to the surface of the bloch
    sphere.

    """
    cos_theta = 2/(cosh(2*epsilon*qs[0]) + cosh(2*epsilon*qs[1]))
    phi = atan2(sinh(epsilon*(qs[0] - qs[1])),
                (sinh(2*epsilon*qs[0]) + sinh(2*epsilon*qs[1]))/2)
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
                    sol = root(lambda qs: map_qpm_to_sphere(qs, epsilon) -
                                np.array([costheta[m,n], phi[m,n]]),
                                guess_qpm_vals(np.array([costheta[m,n],
                                phi[m,n]]), epsilon))
                    inverted[:,m,n] = sign[m,n]*sol.x
                except Exception as e:
                    print('cos(theta) = ' + str(costheta[m,n]) +
                          ', phi = ' + str(phi[m,n]) + ', sign = ' +
                          str(sign[m,n]) + '\n' + str(e))

    return inverted
                
#    return np.where(angles[0] == 1, np.zeros(np.array(angles).shape),
#                   sign*broyden1(lambda qs: prob_dens_q(qs, epsilon) - angles,
#                                  guess_q_vals(angles, epsilon),
#                                  f_tol=1e-14))

def G_qpm(qp, qm, epsilon):
    """The probability density function on the qp qm plane.

    """
    return (epsilon/(4*pi)*exp(-epsilon*(qp**2 + qm**2 + 1))*
            (cosh(2*epsilon*qp) + cosh(2*epsilon*qm)))

def G_angles(angles, epsilon):
    """The probability density function on the upper hemisphere of the bloch
    sphere.

    """
    qp, qm = array_get_qpm_vals(angles, epsilon)
    s_pp = sinh(2*epsilon*qp)
    s_pm = sinh(epsilon*(qp + qm))
    s_mp = sinh(epsilon*(qp - qm))
    s_mm = sinh(2*epsilon*qm)
    c_pp = cosh(2*epsilon*qp)
    c_mp = cosh(epsilon*(qp - qm))
    c_mm = cosh(2*epsilon*qm)
    parallelogram_area = (np.abs((s_pp + s_mm)**2*c_mp +
                          2*(s_pp - s_mm)*s_mp*c_pp) /
                          ((c_pp + c_mm)**2*(4*s_pm**2 + s_mp**2 +
                                             s_pp**2 + s_mm**2)))
    return G_qpm(qp, qm, epsilon)/parallelogram_area

N = 16
Theta, Phi = np.mgrid[0:pi/2:complex(0, N), 0:2*pi:complex(0, 2*N + 1)]
# Theta should be in the open interval (0, pi/2)
Theta = Theta[1:-1,:]
Phi = Phi[1:-1,:]

X = sin(Theta) * cos(Phi)
Y = sin(Theta) * sin(Phi)
Z = cos(Theta)

epsilon = 1

pairs = np.dstack((cos(Theta), Phi))
angles = np.array([cos(Theta), Phi])

colorfunction=G_angles(angles, epsilon)

norm = colors.Normalize()
