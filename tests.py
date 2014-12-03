from nose.tools import assert_almost_equal, assert_equal
import numpy as np
import bloch_distribution.invert_angles as inv

def check_inversions(Q1, Q2, epsilon):
    Qs = np.array([Q1, Q2])
    angles = inv.map_q12_to_sphere_y(Qs, epsilon)
    inverted = inv.array_get_q12_vals(angles, epsilon)
    diffs = Qs - inverted
    assert_almost_equal(np.max(np.abs(diffs)), 0, 7)

# This currently fails if N = 2**3, and when N = 2**2 it fails for epsilon = 1

def test_inversion():
    N = 2**2    # Make mgrids without q1 = q2 = 0, since that case can be
                # problematic
    Q1, Q2 = np.mgrid[-8:8:complex(0, N), -8:8:complex(0, N)]
    epsilons = np.linspace(.1, .7, N)
    for epsilon in epsilons:
        check_inversions(Q1, Q2, epsilon)
