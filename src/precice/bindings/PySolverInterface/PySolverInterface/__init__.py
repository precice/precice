# -*- coding: utf-8 -*-
"""PySolverInterface

PySolverInterface is deprecated and will be removed in preCICE version 2.0.0, please use the python module precice instead.
"""

import warnings

warnings.warn("PySolverInterface is deprecated and will be removed in preCICE version 2.0.0, please use the python module precice instead.", Warning, stacklevel=2)

from PySolverInterface.PySolverInterface import PySolverInterface as PySolverInterface  # black magic, but works and will be deprecated in 2.0.0
from precice import action_write_initial_data as PyActionWriteInitialData
from precice import action_write_iteration_checkpoint as PyActionWriteIterationCheckpoint
from precice import action_read_iteration_checkpoint as PyActionReadIterationCheckpoint
from precice import action_plot_output as PyActionPlotOutput
from precice import export_vtk as PyExportVTK
from precice import export_all as PyExportAll
