"""
Functions to check the type of input argument against what is expected
and raise an Exception if mismatch found
"""
import numpy as np

def check1D(self, func_name, usr_input, true_val):
    if (usr_input != true_val):
        raise Exception("Function {} received position of type {} but expects a 1D numpy \
                         array with {} as position input".format(func_name, usr_input, true_val))

def check2D(self, func_name, usr_input, true_val):
    if (usr_input != true_val):
        raise Exception("Function {} received positions of dimension {}expects a numpy array \
                         [no. of vertices x {}] as positions input".format(func_name, usr_input, true_val))