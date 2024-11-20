import geomproc as gp
import numpy as np
import random
import math

def generate_random_rgb():
    """
    Generate a random RGB value and return it as a NumPy array of size 3.
    """
    return np.array([random.randint(0, 255) for _ in range(3)])

tm = gp.load('firelink.obj')
tm.fcolor = np.zeros((len(tm.face), 3), dtype=np.single)
tm.normalize()
tm.compute_connectivity()
tm.compute_vertex_and_face_normals()

tm_meshlet = tm.generate_meshlets('nvidia')

for meshlet in tm_meshlet.meshlet:
    color = generate_random_rgb()
    for face in meshlet.prim_buffer:
        tm.fcolor[face] = color

wo = gp.write_options()
wo.write_face_colors = True

tm.save('output.ply', wo)


