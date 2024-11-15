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
        if algorithm == 'nvidia':
            return False
        
    
    def tipsify(self, I, k):
        self.mesh.compute_vif()
        L = [len(adjacent_triangles) for adjacent_triangles in self.mesh.vif]
        C = [0 for vertex in self.mesh.vertex]
        D = []
        E = [False for faces in self.mesh.face]
        O = []
        f = 0
        s = k+1
        i = 1

        while f >= 0:
            N = {}
            for t in self.mesh.vif[f]:
                if not E[t]:
                    for v in self.mesh.face[t]:
                        O.append(v)
                        D.append(v)
                        N.add(v)
                        L[v] = L[v]-1
                        if s - C[v] > k:
                            C[v] = s
                            s +=1
            f = getNextVertex(i, k, N, C, s, L, D)
        return O

    def getNextVertex(i, k, N, C, s, L, D):
        n = -1
        p = -1
        for v in N:
            if L[v] > 0:
                p = 0
                if s-C[v]+2*L[v] <= k:
                    p = s-C[v]
                if p > m:
                    m = p
                    n = v
        if n == -1:
            n = skipDeadEnd(L, D, i)
        return n

    def skipDeadEnd(L, D, i):
        while D:
            d = D.pop()
            if L[d] > 0:
                return d
        while i < len(self.mesh.vertex):
            i = i+1
            if L[i] > 0:
                return i
        return -1

