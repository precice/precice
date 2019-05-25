# distutils: language = c++

# comments on test layout: https://docs.pytest.org/en/latest/goodpractices.html
# run with python -m unittest tests.test_fenicsadapter

cimport precice_future
from unittest.mock import MagicMock, patch
from unittest import TestCase
import numpy as np


class TestBindings(TestCase):
    """
    Test suite to check correct behaviour of python bindings.
    """
    
    def test_constructor(self):
        solver_interface = precice_future.Interface("test", 0, 1)
        self.assertTrue(True)
        return_value = solver_interface.read_block_scalar_data(1, np.array([2, 3]))
