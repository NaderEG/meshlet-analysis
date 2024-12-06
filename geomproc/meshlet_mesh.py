# All code was initially written by Nader El-Ghotmi
# But was then sometimes reviewed, refined, or editied by ChatGPT.
# This file also contains the initial implementation of a bounding sphere
# Algorithm but I was not able to complete it in time.



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
import sys


from collections import deque
from collections import defaultdict
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
        self.meshlet = []

        if method == 'nvidia':
            meshlet_vertices, meshlet_triangles = self.tipsify(25)
            

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
            vertex_list = self.sort_by_bounding_box_axis()

            used_vertices = set()
            used_triangles = set()

            vertex_buffer = set()
            triangle_buffer = set()
            local_border = set()

            for idx in vertex_list:
                if idx in used_vertices:
                    continue

                # Start a new meshlet
                local_border.add(idx)

                while local_border:
                    # Select a vertex from the border (randomized choice for variety)
                    curr = random.choice(list(local_border))
                    local_border.remove(curr)

                    if curr in used_vertices:
                        continue

                    # Mark the vertex as used
                    used_vertices.add(curr)

                    # Process triangles connected to the current vertex
                    for triangle in mesh.vif[curr]:
                        if triangle in used_triangles:
                            continue

                        triangle_vertices = mesh.face[triangle]
                        new_vertices = [v for v in triangle_vertices if v not in vertex_buffer]

                        # Check if adding this triangle exceeds buffer limits
                        if len(vertex_buffer) + len(new_vertices) > max_vertices or len(triangle_buffer) + 1 > max_triangles:
                            # Attempt to fill remaining space with border triangles
                            if len(triangle_buffer) < max_triangles:
                                for v in vertex_buffer:
                                    for border_triangle in mesh.vif[v]:
                                        if border_triangle in used_triangles:
                                            continue
                                        if set(mesh.face[border_triangle]).issubset(vertex_buffer):
                                            triangle_buffer.add(border_triangle)
                                            used_triangles.add(border_triangle)

                            # Finalize the current meshlet
                            self.meshlet.append(meshlet(vertex_buffer.copy(), triangle_buffer.copy()))
                            vertex_buffer.clear()
                            triangle_buffer.clear()
                            local_border.clear()

                            # Start a new meshlet with the current vertex
                            local_border.add(curr)
                            break

                        # Add the triangle and its vertices
                        triangle_buffer.add(triangle)
                        used_triangles.add(triangle)
                        vertex_buffer.update(triangle_vertices)

                        # Add new vertices to the border
                        for vertex in triangle_vertices:
                            if vertex not in used_vertices:
                                local_border.add(vertex)

                # Finalize any remaining data in the buffers after the loop
                if vertex_buffer or triangle_buffer:
                    self.meshlet.append(meshlet(vertex_buffer.copy(), triangle_buffer.copy()))
                    vertex_buffer.clear()
                    triangle_buffer.clear()

        if method == 'min_curve':
            meshlet_vertices, meshlet_triangles = self.tipsify(25)
            visited_triangles = set()
            border_triangles = []

            for triangle in meshlet_triangles:
                # If we have border triangles, prioritize them
                if border_triangles:
                    triangle = border_triangles.pop(0)
                    if triangle in visited_triangles:
                        continue
                elif triangle in visited_triangles:
                    # Fallback to the next unvisited triangle if no valid border exists
                    continue

                # Initialize buffers for the new meshlet
                vertex_buffer = set(self.mesh.face[triangle])
                triangle_buffer = {triangle}
                visited_triangles.add(triangle)

                # Calculate the average normal for the initial triangle
                avg_normal = self.mesh.fnormal[triangle]

                # Maintain candidate triangles and normals
                candidate_triangles = []
                candidate_normals = []

                while len(vertex_buffer) < max_vertices and len(triangle_buffer) < max_triangles:
                    # Collect new candidates from the current meshlet boundary
                    for tri in triangle_buffer:
                        for adj_tri in self.mesh.fif[tri]:
                            if adj_tri not in visited_triangles and adj_tri not in candidate_triangles:
                                candidate_triangles.append(adj_tri)
                                candidate_normals.append(self.mesh.fnormal[adj_tri])

                    # Break if no valid candidates are found
                    if not candidate_triangles:
                        break

                    # Find the closest normal to the current average
                    closest_normal, closest_idx = self.find_closest_normal(avg_normal, candidate_normals)
                    selected_triangle = candidate_triangles.pop(closest_idx)
                    candidate_normals.pop(closest_idx)

                    # Update the average normal
                    avg_normal = self.weighted_mean_vector(avg_normal, closest_normal, len(triangle_buffer))

                    # Add the selected triangle and its vertices to the buffers
                    triangle_buffer.add(selected_triangle)
                    visited_triangles.add(selected_triangle)
                    for v in self.mesh.face[selected_triangle]:
                        vertex_buffer.add(v)

                    # Update the border: Add neighbors of the new triangle if not visited
                    for adj_tri in self.mesh.fif[selected_triangle]:
                        if adj_tri not in visited_triangles and adj_tri not in candidate_triangles:
                            border_triangles.append(adj_tri)

                # Finalize the current meshlet
                self.meshlet.append(meshlet(list(vertex_buffer), list(triangle_buffer)))

                # Clean up the border: Remove triangles that were added to the meshlet
                border_triangles = [tri for tri in border_triangles if tri not in triangle_buffer]
            
        if method == 'bounding_sphere':
            epsilon = 0.0001
            vertex_list = self.sort_by_bounding_box_axis()
            radius = 0
            center = np.array([0, 0, 0])
            used_vertices = set()
            used_triangles = set()

            for vertex in vertex_list:
                if vertex in used_vertices:
                    continue

                vertex_buffer = set()
                triangle_buffer = set()
                border = set(self.mesh.vif[vertex])
                best_triangle = None
                new_vertex = -1
                new_radius = sys.float_info.max
                best_new_radius = new_radius - 1
                vertex_score = 0
                best_vertex_score = 0

                for triangle in border:
                    if triangle in used_triangles:
                        continue
                    for vert in self.mesh.face[triangle]:
                        if vert in vertex_buffer:
                            vertex_score+=1
                        else:
                            new_vertex = vertex
                    
                    if vertex_score == 3:
                        new_radius = radius
                    elif vertex_score == 1:
                        continue
                    else:
                        new_radius = 0.5*np.linalg.norm(radius + radius - self.mesh.vertex[vertex])
                    triangle_count = 0
                    for tri in self.mesh.fif[triangle]:
                        if tri in triangle_buffer:
                            triangle_count+=1
                    if len(self.mesh.fif[triangle]) == triangle_count:
                        vertex_score+=1
                    if vertex_score >= best_vertex_score or new_radius <= best_new_radius:
                        best_vertex_score = vertex_score
                        best_new_radius = new_radius
                        best_triangle = triangle
                if best_triangle == None:
                    for triangle in self.mesh.vif[vertex]:
                        if triangle in used_triangles:
                            continue
                        best_triangle = triangle
                        center = 0
                        for v in self.mesh.face[best_triangle]:
                            center+=self.mesh.vertex[v]
                        center = center/3
                        best_new_radius = max((max(center - self.mesh.vertex[self.mesh.face[best_triangle][0]]), 
                        max(center - self.mesh.vertex[self.mesh.face[best_triangle][1]]), 
                        max(center - self.mesh.vertex[self.mesh.face[best_triangle][2]])))

                if best_triangle ==  None:
                    continue
                radius = best_new_radius
                center = vertex_list[new_vertex] + (radius / (epsilon + np.linalg.norm(center - self.mesh.vertex[new_vertex]))) * (self.mesh.vertex[new_vertex] - center)
                if len(vertex_buffer) >= max_vertices:
                    if len(triangle_buffer) < max_triangles:
                        for triangle in border:
                            if set(self.mesh.face[triangle]).issubset(vertex_buffer):
                                triangle_buffer.add(triangle)
                                border.remove(triangle)
                                used_triangles.add(triangle)
                                border.update(self.mesh.fif[triangle])
                                vertex_buffer.update(self.mesh.face[triangle])
                    self.meshlet.append(meshlet(list(vertex_buffer), list(triangle_buffer)))
                    continue
                triangle_buffer.add(triangle)
                vertex_buffer.update(self.mesh.face[triangle])
                used_triangles.add(triangle)
                used_vertices.add()
                border.remove(triangle)
                border.update(self.mesh.fif[triangle])










                


            







            


    def to_obj(self, output_dir):
        os.makedirs(output_dir, exist_ok=True)

        meshlets_dir = os.path.join(output_dir, "meshlets")
        os.makedirs(meshlets_dir, exist_ok=True)

        for file in os.listdir(meshlets_dir):
            file_path = os.path.join(meshlets_dir, file)
            if os.path.isfile(file_path):
                os.remove(file_path)

        for i, meshlet in enumerate(self.meshlet):
            global_vertex_indices = meshlet.vertex_buffer
            global_face_indices = meshlet.prim_buffer

            local_vertex_mapping = {v: idx + 1 for idx, v in enumerate(global_vertex_indices)}

            local_vertices = [self.mesh.vertex[v] for v in global_vertex_indices]

            local_triangles = []
            for face_idx in global_face_indices:
                global_face = self.mesh.face[face_idx]
                local_triangles.append(tuple(local_vertex_mapping[v] for v in global_face))

            filename = os.path.join(meshlets_dir, f"meshlet_{i+1}.obj")

            with open(filename, 'w') as obj_file:
                for vertex in local_vertices:
                    obj_file.write(f"v {vertex[0]} {vertex[1]} {vertex[2]}\n")

                for triangle in local_triangles:
                    obj_file.write(f"f {triangle[0]} {triangle[1]} {triangle[2]}\n")


    def to_ply(self, output_dir):
        os.makedirs(output_dir, exist_ok=True)
        for meshlet in self.meshlet:
            color = self.generate_random_rgb()
            for face in meshlet.prim_buffer:
                self.mesh.fcolor[face] = color
        wo = write_options()
        wo.write_face_colors = True
        self.mesh.save(output_dir+'/output.ply', wo)


        
    
    def generate_random_rgb(self):
        """
        Generate a random RGB value and return it as a NumPy array of size 3.
        """
        return np.array([random.randint(0, 255) for _ in range(3)])
    
    def average_triangles(self):
        total_tris = 0
        for meshlet in self.meshlet:
            total_tris += len(meshlet.prim_buffer)

        return total_tris/len(self.meshlet)
    
    def average_vertices(self):
        total_verts = 0
        for meshlet in self.meshlet:
            total_verts += len(meshlet.vertex_buffer)

        return total_verts/len(self.meshlet)
    
    def shared_vertices(self):
        vertex_occurrences = defaultdict(int)

        # Count how many times each vertex appears in the meshlets
        for meshlet in self.meshlet:
            for vertex in meshlet.vertex_buffer:
                vertex_occurrences[vertex] += 1

        # Count vertices that appear in more than one meshlet
        shared_vertices = sum(1 for count in vertex_occurrences.values() if count > 1)

        return shared_vertices





            
                            







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

    def find_closest_normal(self, avg_normal, candidate_normals):
        # Normalize candidate normals (in case they are not unit vectors)
        normalized_candidates = candidate_normals / np.linalg.norm(candidate_normals, axis=1)[:, np.newaxis]

        # Compute dot products
        dot_products = np.dot(normalized_candidates, avg_normal)

        # Find the index of the maximum dot product
        closest_index = np.argmax(dot_products)

        # Return the closest normal index
        return candidate_normals[closest_index], closest_index

    def weighted_mean_vector(self, v1, v2, n):
        weighted_sum = n * np.array(v1) + np.array(v2)
        mean_vector = weighted_sum / np.linalg.norm(weighted_sum)
        return mean_vector
