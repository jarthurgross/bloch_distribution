from nose.tools import assert_almost_equal, assert_equal, assert_greater
from nose.tools import assert_less_equal
import numpy as np
import bloch_distribution.invert_angles as inv
import bloch_distribution.calculate_densities as den

def check_inversions(Q1, Q2, epsilon):
    Qs = np.array([Q1, Q2])
    angles = inv.map_q12_to_sphere_y(Qs, epsilon)
    inverted = inv.array_get_q12_vals(angles, epsilon)
    diffs = Qs - inverted
    assert_almost_equal(np.max(np.abs(diffs)), 0, 7)

def check_bloch_normalization(N, epsilon):
    Theta, Phi = den.construct_grid(N)
    Angles = np.array([np.cos(Theta), Phi])
    Densities = inv.G_angles_q12(Angles, epsilon)
    Scaled_densities = Densities*np.sin(Theta)
    Total_prob = np.trapz(np.trapz(Scaled_densities, Phi), Theta[:,0])
    assert_greater(Total_prob, 0.95)
    assert_less_equal(Total_prob, 1)

def check_plane_normalization(Q1, Q2, epsilon):
    Densities = inv.G_q12(Q1, Q2, epsilon)
    Total_prob = np.trapz(np.trapz(Densities, Q2), Q1[:,0])
    assert_greater(Total_prob, 0.95)
    assert_less_equal(Total_prob, 1)

# This currently fails if N = 2**3, and when N = 2**2 it fails for epsilon = 1

def test_inversion():
    N = 2**2    # Make mgrids without q1 = q2 = 0, since that case can be
                # problematic
    Q1, Q2 = np.mgrid[-8:8:complex(0, N), -8:8:complex(0, N)]
    epsilons = np.linspace(.1, .7, N)
    for epsilon in epsilons:
        check_inversions(Q1, Q2, epsilon)

def test_bloch_normalization():
    """Test to make sure distributions on the Bloch sphere integrate to 1.

    """

    N = 2**6
    epsilons = np.linspace(.1, 1, 4)
    for epsilon in epsilons:
        check_bloch_normalization(N, epsilon)

def test_plane_normalization():
    """Test to make sure distributions on the plane integrate to 1.

    """

    N = 2**8 + 1
    # Need different integration bounds for different epsilons, since the
    # distribution spreads out more for smaller epsilons
    bounds = [16, 8, 8, 4]
    epsilons = np.linspace(.1, .7, 4)
    for bound, epsilon in zip(bounds, epsilons):
        Q1, Q2 = np.mgrid[-bound:bound:complex(0, N),
                          -bound:bound:complex(0, N)]
        check_plane_normalization(Q1, Q2, epsilon)
