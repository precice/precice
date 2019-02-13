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

def NameConfiguration ():
   return nameConfiguration()

def DataDisplacements ():
   return dataDisplacements()

def DataForces ():
   return dataForces()

def DataVelocities ():
   return dataVelocities()

def ActionWriteInitialData ():
   return actionWriteInitialData()
   
def ActionWriteIterationCheckpoint ():
   return actionWriteIterationCheckpoint()

def ActionReadIterationCheckpoint ():
   return actionReadIterationCheckpoint()

def ActionPlotOutput ():
   return actionPlotOutput()

def ExportVTK ():
   return exportVTK()

def ExportAll ():
   return exportAll()
