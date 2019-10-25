import os
import shutil
import argparse
import numpy as np
import precice_future as precice
import subprocess
import os.path

import re
import numpy as np
import matplotlib.pyplot as plt

import fileinput

configuration_file_name = "precice-config.xml"
participant_name = "dummy_Solver_Two"
mesh_name = "dummy_Mesh_Two"

n = 3
nodeVertexIDs = n*[0.0]

solver_process_index = 0
solver_process_size = 1

dt = 1

interface = precice.Interface(participant_name, solver_process_index, solver_process_size)
interface.configure(configuration_file_name)
dimensions = interface.get_dimensions()
    
mesh_id = interface.get_mesh_id(mesh_name)

vertices = np.zeros((n, dimensions))

i = 0
for x in range(0,n):
    for y in range(0,3):
        vertices[x,y] = i + 1
    i = i + 1

print("vertices are ",vertices)

solver_One_Data = np.zeros(n)
solver_Two_Data = np.zeros(n)

for x in range(0, n):
    solver_Two_Data[x]=1

data_indices = interface.set_mesh_vertices(mesh_id, vertices)
print("data_indices are ",data_indices)

solver_One_Data_ID = interface.get_data_id("solver_One_Data",mesh_id)
solver_Two_Data_ID = interface.get_data_id("solver_Two_Data",mesh_id)

dt = interface.initialize()

while interface.is_coupling_ongoing():
   
    print("DUMMY: Advancing in time")
    for i in range(0,n): 
        solver_One_Data[i] = interface.read_scalar_data(solver_One_Data_ID, data_indices[i])            
        print("solver_One_Data = ", solver_One_Data[i])
        print("data_indices[i] = ", data_indices[i])
        solver_Two_Data[i] = solver_One_Data[i] + 1

    print("DUMMY: Writing iteration checkpoint")
    interface.fulfilled_action(precice.action_write_iteration_checkpoint())
    for i in range(0,n):    
        print("solver_Two_Data = ", solver_Two_Data[i])
        print("data_indices[i] = ", data_indices[i])
        interface.write_scalar_data(solver_Two_Data_ID, data_indices[i], solver_Two_Data[i])

    input("Press Enter to continue...")
    dt = interface.advance(dt)

interface.finalize()
print("DUMMY: Closing python solver dummy...")
