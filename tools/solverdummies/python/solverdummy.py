from __future__ import division

import os
import sys
import argparse
import numpy as np
 
from mpi4py import MPI
import precice
from precice import *

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

configFileName = args.configurationFileName
participantName = args.participantName
meshName = args.meshName

N = 1

solverProcessIndex = 0
solverProcessSize = 1

interface = Interface(participantName, solverProcessIndex, solverProcessSize)
interface.configure(configFileName)
    
meshID = interface.get_mesh_id(meshName)

dimensions = interface.get_dimensions()
vertex = np.zeros(dimensions)
dataIndices = np.zeros(N)

interface.set_mesh_vertices(meshID, N, vertex, dataIndices)

dt = interface.initialize()
    
while interface.is_coupling_ongoing():
   
    if interface.is_action_required(action_write_iteration_checkpoint()):
        print("DUMMY: Writing iteration checkpoint")
        interface.fulfilled_action(action_write_iteration_checkpoint())
    
    dt = interface.advance(dt)
    
    if interface.is_action_required(action_read_iteration_checkpoint()):
        print("DUMMY: Reading iteration checkpoint")
        interface.fulfilled_action(action_read_iteration_checkpoint())
    else:
        print("DUMMY: Advancing in time")
    
interface.finalize()
print("DUMMY: Closing python solver dummy...")

