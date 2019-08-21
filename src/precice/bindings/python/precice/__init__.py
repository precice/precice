# -*- coding: utf-8 -*-
"""precice

this version of the python bindings is deprecated and will be removed in preCICE version 2.0.0. You may import the updated python module using \"import precice_future as precice\" instead. Refer to the changelog for more information.
"""

import warnings

warnings.warn("This version of the precice bindings is deprecated and will be removed in preCICE version 2.0.0. You may import the updated python module using \"import precice_future as precice\" instead. Refer to the changelog for more information.", Warning, stacklevel=2)

from .precice import Interface
from precice_future import action_write_initial_data
from precice_future import action_write_iteration_checkpoint
from precice_future import action_read_iteration_checkpoint

