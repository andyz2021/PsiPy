import numpy
from psipy.fortran.mapfl import mapfl


def wrap_mapfl(input_file=None):

    # Error if py_rho is not specified or file isn't found
    # Error if both py_b is None and py_pb is None
    # Get the value arrays from the Variable objects
    return mapfl.mapfl(input_file)



