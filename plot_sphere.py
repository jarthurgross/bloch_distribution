from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.tri as mtri
from invert_angles import G_tilde

def points_on_sphere(N):
    """ Generate N evenly distributed points on the unit sphere centered at
    the origin. Uses the 'Golden Spiral'.
    Code by Chris Colbert from the numpy-discussion list.
    """
    phi = (1 + np.sqrt(5)) / 2 # the golden ratio
    long_incr = 2*np.pi / phi # how much to increment the longitude

    dz = 2.0 / float(N) # a unit sphere has diameter 2
    bands = np.arange(N) # each band will have one point placed on it
    z = bands * dz - 1 + (dz/2) # the height z of each band/point
    r = np.sqrt(1 - z*z) # project onto xy-plane
    az = bands * long_incr # azimuthal angle of point modulo 2 pi
    x = r * np.cos(az)
    y = r * np.sin(az)
    return x, y, z

def average_g(triples):
    return np.mean([triple[2] for triple in triples])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X, Y, Z = points_on_sphere(2**12)

upper_x = [x for x, z in zip(X, Z) if z > 0]
upper_y = [y for y, z in zip(Y, Z) if z >0 ]
upper_z = [z for z in Z if z > 0]
upper_triples = np.array(list(zip(upper_x, upper_y, upper_z)))

Triples = np.array(list(zip(X, Y, Z)))

from scipy.spatial import ConvexHull

hull = ConvexHull(Triples)
triangles = hull.simplices
upper_triangles = np.array([triangle for triangle in triangles if
                            np.all(Z[triangle] > 0)])

# triangles = mtri.Triangulation(upper_x, upper_y).triangles
colors = np.array([average_g([Triples[idx] for idx in triangle]) for
                   triangle in upper_triangles])

collec = ax.plot_trisurf(mtri.Triangulation(X, Y, upper_triangles),
        Z, shade=False, cmap=plt.get_cmap('Blues'), array=colors,
        edgecolors='none')
#collec.set_array(colors)
collec.autoscale()

plt.show()
