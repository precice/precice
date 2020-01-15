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

n = 1

solver_process_index = 0
solver_process_size = 1

interface = precice.Interface(participant_name, configuration_file_name, solver_process_index, solver_process_size)
    
mesh_id = interface.get_mesh_id(mesh_name)

dimensions = interface.get_dimensions()
vertices = np.zeros((n, dimensions))

data_indices = interface.set_mesh_vertices(mesh_id, vertices)

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

