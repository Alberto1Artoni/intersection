import os
import subprocess
import numpy as np
from scipy import linalg as la
from scipy.spatial import Delaunay
from matplotlib.tri import Triangulation
# import Ofpp
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
os.chdir(os.path.dirname(os.path.abspath(__file__)))

def writeFacet(fl, v1, v2, v3, normal):
    s1, s2, s3, s4 = ('', '', '', '') 
    for i in range(3):
        s1 += str(normal[i]) + ' '
        s2 += str(v1[i]) + ' '
        s3 += str(v2[i]) + ' '
        s4 += str(v3[i]) + ' '

    fl.write(' facet normal ' + s1 +'\n')
    fl.write('  outer loop\n')
    fl.write('   vertex ' + s2 +'\n')
    fl.write('   vertex ' + s3 +'\n')
    fl.write('   vertex ' + s4 +'\n')
    fl.write('  endloop\n')
    fl.write(' endfacet\n')


def writeSTL(name, points, triList):
    file = open(name + '.stl', 'w')
    file.write('solid\n')
    for i in range(np.shape(triList)[0]):
        v1, v2, v3 = (points[triList[i, 0], :],
                      points[triList[i, 1], :],
                      points[triList[i, 2], :])
        normal = np.cross(v2-v1, v3-v1)
        writeFacet(file, v1, v3, v2, normal)
    file.write('endsolid\n')
    file.close()


def parametrizedSurf(t, z):
    r = 1
    x = r * np.cos(t)
    y = r * np.sin(t)
    return np.vstack([np.ravel(x), np.ravel(y), np.ravel(z)]).T


nT, nZ = (1000, 1)
tVec  = np.linspace(0, 2*np.pi, nT)
hT = 2*np.pi / nT
zVec  = np.array([-hT/2, hT/2])


THETA, ZETA = np.meshgrid(tVec, zVec)
THETA[1, :] += hT / 2
THETA, ZETA = np.ravel(THETA), np.ravel(ZETA)

tri = Triangulation(THETA, ZETA)
connList = tri.triangles
points = parametrizedSurf(THETA, ZETA / hT)

writeSTL('cylinder', points, connList)

