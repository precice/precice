from __future__ import division

import argparse
import numpy as np
import precice

parser = argparse.ArgumentParser()
parser.add_argument("configurationFileName", help="Name of the xml config file.", type=str)
parser.add_argument("participantName", help="Name of the solver.", type=str)
parser.add_argument("meshName", help="Name of the mesh.", type=str)

try:
    args = parser.parse_args()
except SystemExit:
    print("")
    print("Usage: python ./solverdummy precice-config participant-name mesh-name")    
    quit()

configuration_file_name = args.configurationFileName
participant_name = args.participantName
mesh_name = args.meshName

n = 3

solver_process_index = 0
solver_process_size = 1

interface = precice.Interface(participant_name, solver_process_index, solver_process_size)
interface.configure(configuration_file_name)
    
mesh_id = interface.get_mesh_id(mesh_name)

dimensions = interface.get_dimensions()
vertex = np.random.rand((n + 1) * dimensions)
data_indices = np.zeros(n)

# add n vertices
interface.set_mesh_vertices(mesh_id, n, vertex[:dimensions*n], data_indices)
for i, idx in enumerate(data_indices):
    assert(idx == i)  # data_indices are initialized in ascending order: [0, 1, ..., n-1]

# add one more vertex
idx = interface.set_mesh_vertex(mesh_id, vertex[dimensions*n:dimensions*(n+1)])
assert(idx == n)  # if we add one more vertex, we get id = n

# make sure that vertex has been appended
n_vertices = interface.get_mesh_vertex_size(mesh_id)
assert(n_vertices == n+1)

# get all vertex positions
position = np.zeros(dimensions * (n+1))
interface.get_mesh_vertices(mesh_id, n+1, np.array(range(n+1)), position)
assert(np.array_equal(position, vertex))   

# get individual positions
for idx in range(n+1):
    position = np.zeros(dimensions)
    interface.get_mesh_vertices(mesh_id, 1, [idx], position)  # TODO: Here we have a failing assertion. This is behavior is incorrect!
    assert(np.array_equal(position, vertex[idx*dimensions:(idx+1)*dimensions]))   

dt = interface.initialize()
    
while interface.is_coupling_ongoing():
   
    if interface.is_action_required(precice.action_write_iteration_checkpoint()):
        print("DUMMY: Writing iteration checkpoint")
        interface.fulfilled_action(precice.action_write_iteration_checkpoint())
    
    dt = interface.advance(dt)
    
    if interface.is_action_required(precice.action_read_iteration_checkpoint()):
        print("DUMMY: Reading iteration checkpoint")
        interface.fulfilled_action(precice.action_read_iteration_checkpoint())
    else:
        print("DUMMY: Advancing in time")
    
interface.finalize()
print("DUMMY: Closing python solver dummy...")

