"""
Functions to check the type of input argument against what is expected
and raise an Exception if mismatch found
"""


def check1D(func_name, usr_input, true_val):
    if usr_input != true_val:
        raise Exception("Function {} received position of type {} but expects a 1D numpy \
                         array with {} values as position input".format(func_name, usr_input, true_val))


def check2D(func_name, usr_input, true_val):
    if usr_input != true_val:
        raise Exception("Function {} received positions of dimension {}expects a numpy array of dimensions \
                         [no. of vertices x {}] as positions input".format(func_name, usr_input, true_val))
