#import os
from libcpp.string cimport string

cdef extern from "./src/precice/Constants.hpp"  namespace "precice::constants":

   const string& nameConfiguration()

   const string& dataDisplacements()
   const string& dataForces()
   const string& dataVelocities()

   const string& actionWriteInitialData()
   const string& actionWriteSimulationCheckpoint()
   const string& actionReadSimulationCheckpoint()
   const string& actionWriteIterationCheckpoint()
   const string& actionReadIterationCheckpoint()
   const string& actionPlotOutput()

   int exportVTK()
   int exportVRML()
   int exportAll()

cdef extern from "./src/cplscheme/Constants.hpp" namespace "precice::cplscheme::constants":
   const string& actionWriteIterationCheckpoint()
   const string& actionReadIterationCheckpoint()
   const string& actionWriteInitialData()

cdef extern from "./src/io/Constants.hpp" namespace "precice::io::constants":
   int exportVTK()
   int exportVRML()
   int exportAll()

#cdef class PyConstants:
#   cdef Constants *thisptr # hold a C++ instance being wrapped
#   cdef Spacetree *spacetreeptr

   # constructor
#   def __cinit__ (self):
#      self.thisptr = new Constants ()
#      self.spacetreepts = new Spacetree ()

#   # destructor
#   def __dealloc__ (self):
#      del self.thisptr

def PyNameConfiguration ():
   return nameConfiguration()

def PyDataDisplacements ():
   return dataDisplacements()

def PyDataForces ():
   return dataForces()

def PyDataVelocities ():
   return dataVelocities()

def PyActionWriteInitialData ():
   return actionWriteInitialData()
   
def PyActionWriteSimulationCheckpoint ():
   return actionWriteSimulationCheckpoint()

def PyActionReadSimulationCheckpoint ():
   return actionReadSimulationCheckpoint()

def PyActionWriteIterationCheckpoint ():
   return actionWriteIterationCheckpoint()

def PyActionReadIterationCheckpoint ():
   return actionReadIterationCheckpoint()

def PyActionPlotOutput ():
   return actionPlotOutput()

def PyExportVTK ():
   return exportVTK()

def PyExportVRML ():
   return exportVRML()

def PyExportAll ():
   return exportAll()
