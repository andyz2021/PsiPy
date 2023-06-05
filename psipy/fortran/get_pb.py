import numpy
from psipy.fortran import getpb


def wrap_get_pb(py_rho=None, py_b = None, py_pb = None, py_help=False, py_verbose=False, py_cubic=False, py_oldmas=False, py_long=0, py_p=0, py_b0=0, py_r=0, py_nx=201, py_ny=201, py_x0=-3, py_x1=3, py_y0=-3, py_y1=3, py_wispr1=False, py_wispr2=False, py_rocc=1, py_dsmult=1, py_power=0, py_disk=-1, py_vf=None, py_vr=None, py_vt=None, py_vp=None, py_scalar=None, py_avg_scalar=None, py_avg_los_angle=None, py_avg_vlos=None, py_avg_vx=None, py_avg_vy=None, py_avg_using_b=False, py_he_frac=0, py_mu = 0.63):

    # Error if py_rho is not specified or file isn't found
    # Error if both py_b is None and py_pb is None
    # Get the value arrays from the Variable objects
    return getpb.getpb(py_rho, py_b, py_pb, py_help, py_verbose, py_cubic, py_oldmas, py_long, py_p, py_b0, py_r, py_nx, py_ny, py_x0, py_x1, py_y0, py_y1, py_wispr1, py_wispr2, py_rocc, py_dsmult, py_power, py_disk, py_vf, py_vr, py_vt, py_vp, py_scalar, py_avg_scalar, py_avg_los_angle, py_avg_vlos, py_avg_vx, py_avg_vy, py_avg_using_b, py_he_frac, py_mu)



