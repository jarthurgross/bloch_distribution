#!/usr/bin/python
from numpy import sqrt, cosh, exp, pi, arccos, cos, sin
import numpy as np
from scipy.special import erf
from itertools import product
from BTrees.OOBTree import OOBTree
from .invert_angles import map_q12_to_sphere_z

def rp(q_range, epsilon):
    return (erf(sqrt(epsilon/2)*(q_range[1] - 1)) -
            erf(sqrt(epsilon/2)*(q_range[0] - 1)) +
            erf(sqrt(epsilon/2)*(q_range[1] + 1)) -
            erf(sqrt(epsilon/2)*(q_range[0] + 1)))/2


def bin_probability(q_ranges, epsilon):
    """Calculate the probability of making measurement corresponding to q1 and
    q2 being in the specified range

    """
    return rp(q_ranges[0], epsilon)*rp(q_ranges[1], epsilon)/4


def prob_dens(qs, epsilon):
    """The probability density function on the q1 q2 plane. (might be improperly
    normalized)

    """
    q1 = qs[:,:,0]
    q2 = qs[:,:,1]
    return (epsilon/(4*pi)*exp(-epsilon*(q1**2 + q2**2 + 2)/2)*(cosh(epsilon*(q1
        + q2)) + cosh(epsilon*(q1 - q2))))


def avg_bin_val(q_ranges, epsilon):
    #TODO: Implement taking the expected value in the bin
    '''
    # Take pairs of intervals and expand into the points at corners
    corners = np.array([list(product(*bins[:,:,n])) for n in
                        range(bins.shape[2])])
    densities = G(corners, epsilon)
    norms = np.sum(densities)
    return norms
    '''

def center_bin_point(q_ranges, epsilon):
    return zip((q_ranges[0,0] + q_ranges[0,1])/2,
               (q_ranges[1,0] + q_ranges[1,1])/2)


def get_bins(min_q, max_q, lin_samps):
    points = np.linspace(min_q, max_q, lin_samps)
    starts = points[:-1]
    stops = points[1:]
    return np.array(list(zip(zip(*product(starts, starts)),
                             zip(*product(stops, stops)))))


def calc_cumulative_dist(min_q, max_q, lin_samps, epsilon):
    bins = get_bins(min_q, max_q, lin_samps)
    probs = bin_probability(bins, epsilon)
    return dict(zip(np.cumsum(probs), center_bin_point(bins, epsilon)))

def build_tree(min_q, max_q, lin_samps, epsilon):
    return OOBTree(calc_cumulative_dist(min_q, max_q, lin_samps, epsilon))

def get_samples(n, tree):
    """Return n points in q1 q2 space sampled from the specified
    distribution over that space.

    """

    rand_vals = np.random.random(n)*tree.maxKey()
    return [tree.get(tree.minKey(rand_val)) for rand_val in rand_vals]


def get_state_samples(min_q, max_q, lin_samps, epsilon, dist_samps):
    """High level function that constructs all the necessary objects for
    sampling from the appropriate distribution and then does the sampling, with
    samples returned as normalized vectors in C^2 (represented in the
    computational z basis).

    :param min_q:       The minimum value of q1 and q2 to sample from
    :param max_q:       The maximum value of q1 and q2 to sample from
    :param lin_samps:   The number of bin edges in each direction (one more
                        than the number of bins in each direction)
    :param epsilon:     The strength of the weak measurement
    :param dist_samps:  The number of times to sample from the distribution
    :returns:           numpy.array with shape (2,dist_samps) with the first row
                        containing the |0> coefficients and the second row
                        containingthe |1> coefficients

    """

    tree = build_tree(min_q, max_q, lin_samps, epsilon)
    samples = get_samples(dist_samps, tree)
    angles = map_q12_to_sphere_z(np.array(list(zip(*samples))), epsilon)
    Theta = arccos(angles[0])
    Phase = exp(1.j*angles[1])
    return np.array([cos(Theta/2), sin(Theta/2)*Phase])


# TODO: Figure out why calling linspace makes python think sampler needs an
# attribute '__float__'
'''
class sampler:
    def __init__(min, max, lin_samps, epsilon, self):
        self.tree = OOBTree(calc_cumulative_dist(min, max, lin_samps, epsilon))

    def get_samples(n, self):
        """Return n points in q1 q2 space sampled from the specified
        distribution over that space.

        """

        rand_vals = np.random.random(n)*self.tree.maxKey()
        return [self.tree.get(self.tree.minKey(rand_val)) for rand_val in
                rand_vals]
'''
