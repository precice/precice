#import os
#import shutil
import argparse
import numpy as np
import precice_future as precice
#import subprocess
#import os.path

#import re
#import numpy as np
#import matplotlib.pyplot as plt

#import fileinput

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

if (participant_name == 'SolverOne'):
  DataOne_Name='Forces'
  DataTwo_Name='Velocities'

if (participant_name == 'SolverTwo'):
  DataTwo_Name='Forces'
  DataOne_Name='Velocities'

n = 3				#Number of vertices
dt = 1				#Time step
nodeVertexIDs = n*[0.0]		#
solver_process_index = 0	#Rank number
solver_process_size = 1		#Total number of ranks

interface = precice.Interface(participant_name, solver_process_index, solver_process_size)
interface.configure(configuration_file_name)

dimensions = interface.get_dimensions()    
mesh_id = interface.get_mesh_id(mesh_name)
vertices = np.zeros((n, dimensions))

i = 0
for x in range(0,n):
    for y in range(0,dimensions):
        vertices[x,y] = i + 1
    i = i + 1

DataOne = np.zeros(n)
DataTwo = np.zeros(n)

for x in range(0, n):
	DataTwo[x]=1
	DataOne[x]=1

data_indices = interface.set_mesh_vertices(mesh_id, vertices)

DataOne_ID = interface.get_data_id(DataOne_Name,mesh_id)
DataTwo_ID = interface.get_data_id(DataTwo_Name,mesh_id)

dt = interface.initialize()

while interface.is_coupling_ongoing():
	   
	print("DUMMY: Reading iteration checkpoint")
	for i in range(0,n): 
		DataTwo[i] = interface.read_scalar_data(DataTwo_ID, data_indices[i])            

	DataOne = DataTwo + 1

	print("DUMMY: Writing iteration checkpoint")
	interface.fulfilled_action(precice.action_write_iteration_checkpoint())
	for i in range(0,n):    
		interface.write_scalar_data(DataOne_ID, data_indices[i], DataOne[i])

	print("DUMMY: Advancing in time")
	dt = interface.advance(dt)

interface.finalize()
print("DUMMY: Closing python solver dummy...")
