#import os
from libcpp.string cimport string

cdef extern from "precice/Constants.hpp"  namespace "precice::constants":

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

def name_configuration ():
   return nameConfiguration()

def data_displacements ():
   return dataDisplacements()

def data_forces ():
   return dataForces()

def data_velocities ():
   return dataVelocities()

def action_write_initial_data ():
   return actionWriteInitialData()
   
def action_write_iteration_checkpoint ():
   return actionWriteIterationCheckpoint()

def action_read_iteration_checkpoint ():
   return actionReadIterationCheckpoint()

def action_plot_output ():
   return actionPlotOutput()

def export_vtk ():
   return exportVTK()

def export_all ():
   return exportAll()
