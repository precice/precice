from __future__ import division

import os
import sys
import argparse
import numpy as np
 
from mpi4py import MPI
import PySolverInterface
from PySolverInterface import *

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

interface = PySolverInterface(participantName, solverProcessIndex, solverProcessSize)
interface.configure(configFileName)
    
meshID = interface.getMeshID(meshName)

dimensions = interface.getDimensions()
vertex = np.zeros(dimensions)
dataIndices = np.zeros(N)

interface.setMeshVertices(meshID, N, vertex, dataIndices)

dt = interface.initialize()
    
while interface.isCouplingOngoing():
   
    if interface.isActionRequired(PyActionWriteIterationCheckpoint()):
        print("DUMMY: Writing iteration checkpoint")
        interface.fulfilledAction(PyActionWriteIterationCheckpoint())
    
    dt = interface.advance(dt)
    
    if interface.isActionRequired(PyActionReadIterationCheckpoint()):
        print("DUMMY: Reading iteration checkpoint")
        interface.fulfilledAction(PyActionReadIterationCheckpoint())
    else:
        print("DUMMY: Advancing in time")
    
interface.finalize()
print("DUMMY: Closing python solver dummy...")

