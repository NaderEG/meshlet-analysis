import geomproc as gp
import numpy as np
import random
import math
import time
import sys
import os

if len(sys.argv) > 1:
    mesh_name = sys.argv[1]
    if len(sys.argv) > 2:
        method = sys.argv[2]
        if method not in ["greedy", "nvidia", "min_curve"]:
            print("The clustering method provided is not valid. Using default: nvidia.")
            method = "nvidia"
    else:
        print("A clustering method was not provided. Using default: nvidia.")

    file_name = os.path.join("meshes", f"{mesh_name}.obj")
    if os.path.isfile(file_name):
        tm = gp.load(file_name)
        tm.fcolor = np.zeros((len(tm.face), 3), dtype=np.single)
        tm.normalize()
        tm.compute_connectivity()
        tm.compute_vertex_and_face_normals()

        start = time.time()
        tm_meshlet = tm.generate_meshlets(method)
        end = time.time()

        print("Processing time: ", end - start)
        print("Avg triangles: ", tm_meshlet.average_triangles())
        print("Avg vertices: ", tm_meshlet.average_vertices())
        print("Shared vertices: ", tm_meshlet.shared_vertices())
        print("Num Meshlets: ", len(tm_meshlet.meshlet))

        tm_meshlet.to_ply("output")
        tm_meshlet.to_obj("output")


    else:
        print("Mesh "+mesh_name+" does not exist in meshes folder. Terminating...")


else:
    print("No mesh specified. Terminating...")