#
# GeomProc: geometry processing library in python + numpy
#
# Copyright (c) 2008-2024 Oliver van Kaick <ovankaic@gmail.com>
# under the MIT License.
#
# See file LICENSE.txt for details on copyright licenses.
#
"""This module contains the meshlet class of the GeomProc geometry
processing library.
"""

import numpy as np
import copy
import math
import random
import os
import unittest
import json


class meshlet:
    def __init__(self, vertex_buffer=[], prim_buffer=[]):
        self.vertex_buffer = vertex_buffer      # contains the indices of vertices from mesh.vertex
        self.prim_buffer = prim_buffer      # contains the indices of faces from mesh.face
        
    def is_full(self, max_vertices, max_triangles):
        return len(self.vertex_buffer) == max_vertices, len(self.prim_buffer) == max_triangles


    def print(self):
        print(json.dumps(self.__dict__, indent=4))
        

