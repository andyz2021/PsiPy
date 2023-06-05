import numpy
from psipy.fortran import get_va


def wrap_get_va(br, bt, bp, rho):

    # Get the value arrays from the Variable objects
    # Transpose them so that they work for Fortran
    br = br.data.values
    # br = numpy.squeeze(br, axis=3)
    # br = numpy.transpose(br, (2, 1, 0))
    bt = bt.data.values
    # bt = numpy.squeeze(bt, axis=3)
    # bt = numpy.transpose(bt, (2, 1, 0))
    bp = bp.data.values
    # bp = numpy.squeeze(bp, axis=3)
    # bp = numpy.transpose(bp, (2, 1, 0))
    rho = rho.data.values
    # rho = numpy.squeeze(rho, axis=3)
    # rho = numpy.transpose(rho, (2, 1, 0))

    #Convert Real units back to Mas units
    b_unit = 2.2068914
    # b_unit = 1
    br = br / b_unit
    bt = bt / b_unit
    bp = bp / b_unit

    rho_unit = 1.6726*10**-16
    # rho_unit = 1
    rho = rho / rho_unit
    # print(br.shape)
    # print(bt.shape)
    # print(bp.shape)
    # print(rho.shape)
    # Set sizes
    # Dimension order is reversed in Python
    nr = br.shape[2]
    nt = bt.shape[1]
    np = bp.shape[0]

    # Dimensions for va are nr-1, nt-1, np-1
    nr = nr-1
    nt = nt-1
    np = np-1
    #
    # print("nr: ", nr)
    # print("nt: ", nt)
    # print("np: ", np)
    # a = 29
    # b = 24
    # c = 19
    # x = 0.5 * (br[29][24][19] + br[a][b-1][c])
    # print(x)
    # y = 0.5 * (bt[a][b][c] + bt[a][b][c-1])
    # print(y)
    # z = 0.25 * (bp[a][b][c] + bp[a][b-1][c] + bp[a][b][c-1] + bp[a][b-1][c-1])
    # print(z)
    # bsq = x**2 + y**2 + z**2
    # print("bsq: ", bsq)
    # print("rho: ", rho[a][b][c])
    # print("va: ", numpy.sqrt(bsq / abs(rho[a][b][c])))


    # print(type(get_va.getva(br, bt, bp, rho, np, nt, nr)))

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