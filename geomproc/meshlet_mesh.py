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
