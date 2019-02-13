import warnings

warnings.warn("PySolverInterface is deprecated and will be removed in preCICE version 2.0.0, please use the python module precice instead", DeprecationWarning, stacklevel=2)

from precice import Interface as PySolverInterface
from precice import nameConfiguration as PyNameConfiguration
from precice import dataDisplacements as PyDataDisplacements
from precice import dataForces as PyDataForces
from precice import dataVelocities as PyDataVelocities
from precice import actionWriteInitialData as PyActionWriteInitialData
from precice import actionWriteIterationCheckpoint as PyActionWriteIterationCheckpoint
from precice import actionReadIterationCheckpoint as PyActionReadIterationCheckpoint
from precice import actionPlotOutput as PyActionPlotOutput
from precice import exportVTK as PyExportVTK
from precice import exportAll as PyExportAll
