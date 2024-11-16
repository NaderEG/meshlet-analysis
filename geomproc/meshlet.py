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


class meshlet:
    def __init__(self, vertex_count, prim_count, vertex_begin, prim_begin):
        self.vertex_count = vertex_count
        self.prim_count = prim_count
        self.vertex_begin = vertex_begin
        self.prim_begin = prim_begin

        

