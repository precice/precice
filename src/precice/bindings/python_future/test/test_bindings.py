# comments on test layout: 
# here https://docs.pytest.org/en/latest/goodpractices.html
# and here https://stackoverflow.com/questions/42259741/how-to-test-functions-cdefd-in-cython
# TODO run with python -m unittest tests.test_bindings

import pyximport; pyximport.install()
from test_bindings_module import *

import unittest


