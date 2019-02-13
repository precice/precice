#import os
from libcpp.string cimport string

cdef extern from "./src/precice/Constants.hpp"  namespace "precice::constants":

   const string& nameConfiguration()

   const string& dataDisplacements()
   const string& dataForces()
   const string& dataVelocities()

   const string& actionWriteInitialData()
   const string& actionWriteIterationCheckpoint()
   const string& actionReadIterationCheckpoint()
   const string& actionPlotOutput()

   int exportVTK()
   int exportAll()

cdef extern from "./src/cplscheme/Constants.hpp" namespace "precice::cplscheme::constants":
   const string& actionWriteIterationCheckpoint()
   const string& actionReadIterationCheckpoint()
   const string& actionWriteInitialData()

cdef extern from "./src/io/Constants.hpp" namespace "precice::io::constants":
   int exportVTK()
   int exportAll()

def nameConfiguration ():
   return nameConfiguration()

def dataDisplacements ():
   return dataDisplacements()

def dataForces ():
   return dataForces()

def dataVelocities ():
   return dataVelocities()

def actionWriteInitialData ():
   return actionWriteInitialData()
   
def actionWriteIterationCheckpoint ():
   return actionWriteIterationCheckpoint()

def actionReadIterationCheckpoint ():
   return actionReadIterationCheckpoint()

def actionPlotOutput ():
   return actionPlotOutput()

def exportVTK ():
   return exportVTK()

def exportAll ():
   return exportAll()
