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

    def __init__(self, mesh, method, max_vertices, max_triangles):
        self.mesh = mesh

        if method == 'nvidia':
            meshlet_vertices, meshlet_triangles = self.tipsify(15)
            self.meshlet = []

            vertex_buffer = set()
            triangle_buffer = []
            current_triangle = 0

            while current_triangle < len(meshlet_triangles):
                # Get vertices of the current triangle
                triangle_vertices = set(mesh.face[meshlet_triangles[current_triangle]])
                
                # Check if adding this triangle violates meshlet limits
                if len(vertex_buffer | triangle_vertices) <= max_vertices and len(triangle_buffer) < max_triangles:
                    vertex_buffer |= triangle_vertices  # Add triangle vertices
                    triangle_buffer.append(meshlet_triangles[current_triangle])  # Add triangle
                    current_triangle += 1
                else:
                    # Store the current meshlet
                    self.meshlet.append(meshlet(list(vertex_buffer), triangle_buffer))
                    # Reset buffers for the next meshlet
                    vertex_buffer = set()
                    triangle_buffer = []

            # Handle the last meshlet
            if triangle_buffer:
                self.meshlet.append(meshlet(list(vertex_buffer), triangle_buffer))

        if method == 'greedy':
            vertex_index = self.sort_by_bounding_box_axis()
            self.meshlet = []
            visited_triangles = set()
            visited_vertices = set()

            for idx in vertex_index:
                if idx in visited_vertices:
                    continue

                vertex_buffer = set()
                triangle_buffer = set()
                border_triangles = set(self.mesh.vif[idx])

                vertex_buffer.add(idx)
                visited_vertices.add(idx)

                while border_triangles:
                    tri_idx = border_triangles.pop()
                    if tri_idx in visited_triangles:
                        continue

                    triangle = self.mesh.face[tri_idx]
                    new_vertices = [v for v in triangle if v not in vertex_buffer]

                    if (len(vertex_buffer) + len(new_vertices) <= max_vertices and
                        len(triangle_buffer) + 1 <= max_triangles):
                        triangle_buffer.add(tri_idx)
                        visited_triangles.add(tri_idx)

                        for v in triangle:
                            vertex_buffer.add(v)
                            visited_vertices.add(v)

                        for v in triangle:
                            for neighbour_tri in self.mesh.vif[v]:
                                if neighbour_tri not in visited_triangles:
                                    border_triangles.add(neighbour_tri)

                self.meshlet.append(meshlet(list(vertex_buffer), list(triangle_buffer)))

            # Handle unvisited triangles
            for tri_idx in range(len(self.mesh.face)):
                if tri_idx not in visited_triangles:
                    triangle = self.mesh.face[tri_idx]
                    vertex_buffer = set(triangle)
                    triangle_buffer = {tri_idx}
                    visited_triangles.add(tri_idx)

                    for v in triangle:
                        for neighbour_tri in self.mesh.vif[v]:
                            if neighbour_tri not in visited_triangles:
                                triangle_buffer.add(neighbour_tri)
                                visited_triangles.add(neighbour_tri)
                                vertex_buffer.update(self.mesh.face[neighbour_tri])

                                if len(vertex_buffer) > max_vertices or len(triangle_buffer) > max_triangles:
                                    break

                    self.meshlet.append(meshlet(list(vertex_buffer), list(triangle_buffer)))






            






            

            

                            



                
                
                






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
        T = []
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
                    T.append(t)
            f = self.getNextVertex(i, k, N, C, s, L, D)
        return O, T



    def sort_by_bounding_box_axis(self):
        min_bounds = self.mesh.vertex.min(axis=0)
        max_bounds = self.mesh.vertex.max(axis=0)
    
        # Find the axis with the largest length
        axis_lengths = max_bounds - min_bounds
        longest_axis = np.argmax(axis_lengths)

        return np.argsort(self.mesh.vertex[:, longest_axis])
