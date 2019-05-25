# distutils: language = c++

# comments on test layout: https://docs.pytest.org/en/latest/goodpractices.html
# run with python -m unittest tests.test_fenicsadapter

cimport precice_future
from unittest.mock import MagicMock, patch
from unittest import TestCase

class TestBindings(TestCase):
    """
    Test suite to check correct behaviour of python bindings.
    """

    def test_constructor(self):
        precice_future.Interface("test", 0, 1)
        self.assertTrue(True)
