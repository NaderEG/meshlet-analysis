#
# GeomProc: geometry processing library in python + numpy
#
# Copyright (c) 2008-2024 Oliver van Kaick <ovankaic@gmail.com>
# under the MIT License.
#
# See file LICENSE.txt for details on copyright licenses.
#
"""This module contains the meshlet mesh class of the GeomProc geometry
processing library.
"""

import numpy as np
import copy
import math
import random
import os
import unittest

from .creation import *
from .mesh import *
from .meshlet import *

VERTEX_BUFFER_CAP = 64
PRIMITIVE_BUFFER_CAP = 126

class meshlet_mesh:
    """A class that represents a mesh subdivided into meshlets.
    Notes
    -----
    The class stores information which can be used in conjunction with a Mesh object
    in order to treat clusters of triangles as meshlets.

     Attributes
    ----------
    meshlet : numpy.array_like
    Stores objects of a meshlet class which can be used to reference 
    meshlets from the vertex buffer while being minimally compact.

    vertex_buffer : numpy.array_like
    Stores the vertices of a mesh in a specific order such that
    they can be accessed by the meshlet descriptors.

    primitive_buffer : numpy.array_like
    Stores the primitives (triangles) of a mesh in a specific order
    such that they can be access by the meshlet descriptors.
    """

    def __init__(self, mesh, algorithm='nvidia'):
        self.mesh = mesh
        global_buffer = mesh.vertex

        if algorithm == 'nvidia':
            meshlet_triangles = self.tipsify(126)

            seen = set()
            meshlet_vertices = [x for x in meshlet_triangles if not (x in seen or seen.add(x))]






    def skipDeadEnd(self, L, D, i):
        while D:
            d = D.pop()
            if L[d] > 0:
                return d
        while i < len(self.mesh.vertex) - 1:
            i += 1
            if L[i] > 0:
                return i
        return -1

    def getNextVertex(self, i, k, N, C, s, L, D):
        p = -1
        n = -1
        m = -1
        for v in N:
            if L[v] > 0:
                p = 0
                if s - C[v] + 2 * L[v] <= k:
                    p = s - C[v]
                if p > m:
                    m = p
                    n = v
        if n == -1:
            n = self.skipDeadEnd(L, D, i)
        return n

    def tipsify(self, k):
        self.mesh.compute_vif()
        L = [len(adjacent_triangles) for adjacent_triangles in self.mesh.vif]
        C = [0 for _ in self.mesh.vertex]
        D = []
        E = [False for _ in self.mesh.face]
        O = []
        f = 0
        s = k + 1
        i = 1

        while f >= 0:
            N = set()
            for t in self.mesh.vif[f]:
                if not E[t]:
                    for v in self.mesh.face[t]:
                        O.append(v)
                        D.append(v)
                        N.add(v)
                        L[v] -= 1
                        if s - C[v] > k:
                            C[v] = s
                            s += 1
                    E[t] = True
            f = self.getNextVertex(i, k, N, C, s, L, D)
        return O


    

    
class test_meshlet_mesh(unittest.TestCase):
    def test_tipsify_vertices(self):
        flag = True

        tm = create_torus(1.0, 0.33, 90, 30)
        tm.compute_connectivity()
        tm_meshlet = meshlet_mesh(tm)
        tipsified_vertex_list = [v for v in tm_meshlet.tipsify(256)]
        for i in range(len(tm.vertex)):
            if i not in tipsified_vertex_list:
                flag = False
        
        self.assertEqual(True, flag)

        
        



if __name__ == "__main__":
    unittest.main()




