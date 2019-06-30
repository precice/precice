from __future__ import division

import argparse
import numpy as np
import precice_future as precice

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
vertex_x = np.linspace(0, 1, n + 1)
vertex_y = np.zeros(n+1)
vertex = np.vstack((vertex_x,vertex_y)).ravel(order='F')
data_indices = np.zeros(n)

if participant_name == "SolverOne":
    write_data_name = "Forces"
    read_data_name = "Velocities"
elif participant_name == "SolverTwo":
    write_data_name = "Velocities"
    read_data_name = "Forces"
else:
    raise Exception("unknown participant!")

# add n vertices
data_indices = interface.set_mesh_vertices(mesh_id, vertex[:dimensions*n])
for i, idx in enumerate(data_indices):
    assert(idx == i)  # data_indices are initialized in ascending order: [0, 1, ..., n-1]

print(data_indices.dtype)
# add one more vertex
data_indices = np.append(data_indices, interface.set_mesh_vertex(mesh_id, vertex[dimensions*n:dimensions*(n+1)]))
assert(data_indices[-1] == n)  # if we add one more vertex, we get id = n
print(data_indices.dtype)
# make sure that vertex has been appended
n_vertices = interface.get_mesh_vertex_size(mesh_id)
assert(n_vertices == n+1)

# get all vertex positions
position = interface.get_mesh_vertices(mesh_id, data_indices[:-1])
assert(np.array_equal(position, vertex[:dimensions*n]))

# get individual positions
for i in data_indices:
    position = interface.get_mesh_vertices(mesh_id, np.array([data_indices[i]]))  # TODO: Here we have a failing assertion. This is behavior is incorrect!
    assert(np.array_equal(position, vertex[i*dimensions:(i+1)*dimensions]))   

write_data_id = interface.get_data_id(write_data_name, mesh_id)
read_data_id = interface.get_data_id(read_data_name, mesh_id)

dt = interface.initialize()

write_data = n+2
write_block_data = np.arange(n, dtype=np.double)
write_block_data[-1] = 10
    
while interface.is_coupling_ongoing():
   
    if interface.is_action_required(precice.action_write_iteration_checkpoint()):
        print("DUMMY: Writing iteration checkpoint")
        interface.fulfilled_action(precice.action_write_iteration_checkpoint())

    interface.write_block_scalar_data(write_data_id, data_indices[:-1], write_block_data)
    interface.write_scalar_data(write_data_id, data_indices[-1], write_data)
    dt = interface.advance(dt)
    read_block_data = interface.read_block_scalar_data(read_data_id, data_indices[:-1])
    read_data = interface.read_scalar_data(read_data_id, data_indices[-1])
    for i in range(n):
        assert(read_block_data[i] == write_block_data[i])
    assert(read_data == write_data)
    
    if interface.is_action_required(precice.action_read_iteration_checkpoint()):
        print("DUMMY: Reading iteration checkpoint")
        interface.fulfilled_action(precice.action_read_iteration_checkpoint())
    else:
        print("DUMMY: Advancing in time")
interface.finalize()
print("DUMMY: Closing python solver dummy...")

