import numpy
from psipy.fortran import get_va


def wrap_get_va(br, bt, bp, rho):

    # Get the value arrays from the Variable objects
    br = br.data.values
    bt = bt.data.values
    bp = bp.data.values
    rho = rho.data.values

    #Convert Real units back to Mas units
    b_unit = 2.2068914
    # b_unit = 1
    br = br / b_unit
    bt = bt / b_unit
    bp = bp / b_unit

    rho_unit = 1.6726*10**-16
    # rho_unit = 1
    rho = rho / rho_unit
    print(br.shape)
    print(br[269][99][199])
    # Set sizes
    # Dimension order is reversed in Python
    nr = br.shape[2]
    nt = bt.shape[1]
    np = bp.shape[0]

    # Dimensions for va are nr-1, nt-1, np-1
    nr = nr-1
    nt = nt-1
    np = np-1

    # Return va from the getva function
    v_unit = 481.37107
    # v_unit = 1
    return v_unit * get_va.getva(br, bt, bp, rho, nr, nt, np)

def test_func():
    print(get_va.getva.__doc__)
    br = numpy.random.rand(3, 4, 5)
    bt = numpy.random.rand(3, 4, 5)
    bp = numpy.random.rand(3, 4, 5)
    rho = numpy.random.rand(3, 4, 5)
    va = numpy.zeros((3, 4, 5))
    # print(get_va.getva(br, bt, bp, rho, 3, 4, 5))
    # print(va)
    return va



    # nr = min(br.shape[0], bt.shape[0], bp.shape[0], rho.shape[0])
    # nt = min(br.shape[1], bt.shape[1], bp.shape[1], rho.shape[1])
    # np = min(br.shape[2], bt.shape[2], bp.shape[2], rho.shape[2])

    # Reshape Arrays so that they're all the same size
    # Not sure if this is correct
    # br = numpy.resize(br, (min_shape))
    # bt = numpy.resize(bt, (min_shape))
    # bp = numpy.resize(bp, (min_shape))
    # rho = numpy.resize(rho, (min_shape))