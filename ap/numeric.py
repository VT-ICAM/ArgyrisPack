#! /usr/bin/env python
import numpy as np
import ctypes as ct

# TODO get this path to work in a more general way.
_ap = np.ctypeslib.load_library('libargyris_pack.so', '.')

array_1d_double = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='C_CONTIGUOUS')
array_2d_double = np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='C_CONTIGUOUS')

_ap.ap_ref_values.restype  = None
_ap.ap_ref_values.argtypes = [array_1d_double, array_1d_double,
                              ct.c_int, array_2d_double]

_ap.ap_ref_gradients.restype  = None
_ap.ap_ref_gradients.argtypes = [array_1d_double, array_1d_double,
                                 ct.c_int, array_2d_double,
                                 array_2d_double]

_ap.ap_ref_hessians.restype  = None
_ap.ap_ref_hessians.argtypes = [array_1d_double, array_1d_double,
                                ct.c_int, array_2d_double,
                                array_2d_double, array_2d_double]

_ap.ap_physical_maps.restype  = None
_ap.ap_physical_maps.argtypes = [array_1d_double, array_1d_double,
                                 array_2d_double, array_2d_double,
                                 array_1d_double]

_ap.ap_physical_values.restype  = None
_ap.ap_physical_values.argtypes = [array_2d_double, array_2d_double,
                                  ct.c_int, array_2d_double]

_ap.ap_physical_gradients.restype  = None
_ap.ap_physical_gradients.argtypes = [array_2d_double, array_2d_double,
                                      array_2d_double, array_2d_double,
                                      ct.c_int,
                                      array_2d_double, array_2d_double]

_ap.ap_physical_hessians.restype  = None
_ap.ap_physical_hessians.argtypes = [array_2d_double, array_2d_double,
                                     array_2d_double, array_2d_double,
                                     array_2d_double, ct.c_int,
                                     array_2d_double, array_2d_double,
                                     array_2d_double]

_ap.ap_matrix_mass.restype  = None
_ap.ap_matrix_mass.argtypes = [array_2d_double, array_2d_double,
                               array_2d_double, array_1d_double,
                               ct.c_int, array_2d_double]

_ap.ap_matrix_stiffness.restype  = None
_ap.ap_matrix_stiffness.argtypes = [array_2d_double, array_2d_double,
                                    array_2d_double, array_2d_double,
                                    array_1d_double, ct.c_int,
                                    array_2d_double]

_ap.ap_matrix_betaplane.restype  = None
_ap.ap_matrix_betaplane.argtypes = [array_2d_double, array_2d_double,
                                    array_2d_double, array_2d_double,
                                    array_2d_double, array_1d_double,
                                    ct.c_int, array_2d_double]

_ap.ap_matrix_biharmonic.restype  = None
_ap.ap_matrix_biharmonic.argtypes = [array_2d_double, array_2d_double,
                                     array_2d_double, array_2d_double,
                                     array_2d_double, array_1d_double,
                                     ct.c_int, array_2d_double]

def ref_values(x, y):
    """
    Calculate the values of the Argyris basis functions at given reference
    points.

    Arguments:
    - `x` : 1-dimensional matrix of x-coordinates.
    - `y` : 1-dimensional matrix of y-coordinates.
    """
    check_evaluation_points(x, y)
    values = np.empty((21,x.shape[0]))
    _ap.ap_ref_values(x, y, x.shape[0], values)
    return values

def ref_gradients(x, y):
    """
    Calculate the first derivatives of the Argyris basis functions at given
    reference points.

    Arguments:
    - `x` : 1-dimensional matrix of x-coordinates.
    - `y` : 1-dimensional matrix of y-coordinates.
    """
    check_evaluation_points(x, y)
    ref_dx = np.empty((21,x.shape[0]))
    ref_dy = np.empty((21,x.shape[0]))
    _ap.ap_ref_gradients(x, y, x.shape[0], ref_dx, ref_dy)
    return (ref_dx, ref_dy)

def ref_hessians(x, y):
    """
    Calculate the second derivatives of the Argyris basis functions at given
    reference points.

    Arguments:
    - `x` : 1-dimensional matrix of x-coordinates.
    - `y` : 1-dimensional matrix of y-coordinates.
    """
    check_evaluation_points(x, y)
    ref_dxx = np.empty((21,x.shape[0]))
    ref_dxy = np.empty((21,x.shape[0]))
    ref_dyy = np.empty((21,x.shape[0]))
    _ap.ap_ref_hessians(x, y, x.shape[0], ref_dxx, ref_dxy, ref_dyy)
    return (ref_dxx, ref_dxy, ref_dyy)

def physical_maps(x, y):
    """
    Calculate the Argyris change of basis matrices C, B, and b.

    Arguments:
    - `x` : x-coordinates of the triangle vertices.
    - `y` : y-coordinates of the triangle vertices.
    """
    assert x.shape == (3,) and y.shape == (3,)
    assert x.dtype == np.float64 and y.dtype == np.float64

    C = np.empty((21,21), dtype=np.float64)
    B = np.empty((2,2), dtype=np.float64)
    b = np.empty((2,), dtype=np.float64)
    _ap.ap_physical_maps(x, y, C, B, b)
    return (C, B, b)

def physical_values(C, ref_values):
    """
    Calculate the values of the Argyris basis functions on a physical element.

    Arguments:
    - `C`          : (21, 21) Argyris transformation matrix.
    - `ref_values` : (21, N) matrix of reference function values.
    """
    check_transformations(C)
    check_ref_values(ref_values)
    values = np.empty(ref_values.shape, dtype=np.float64)
    _ap.ap_physical_values(C, ref_values, ref_values.shape[1], values)
    return values

def physical_gradients(C, B, ref_dx, ref_dy):
    """
    Calculate the first derivatives of the Argyris basis functions on a physical
    element.

    Arguments:
    - `C`      : (21, 21) Argyris transformation matrix.
    - `B`      : (2, 2) Affine multiplier matrix.
    - `ref_dx` : (21, N) matrix of reference function x-derivative values.
    - `ref_dy` : (21, N) matrix of reference function y-derivative values.
    """
    check_transformations(C, B)
    check_ref_values(ref_dx, ref_dy)
    dx = np.empty(ref_dx.shape, dtype=np.float64)
    dy = np.empty(ref_dy.shape, dtype=np.float64)
    _ap.ap_physical_gradients(C, B, ref_dx, ref_dy, ref_dx.shape[1], dx, dy)
    return (dx, dy)

def physical_hessians(C, B, ref_dxx, ref_dxy, ref_dyy):
    """
    Calculate the second derivatives of the Argyris basis functions on a physical
    element.

    Arguments:
    - `C`       : (21, 21) Argyris transformation matrix.
    - `B`       : (2, 2) Affine multiplier matrix.
    - `ref_dxx` : (21, N) matrix of reference function xx-derivative values.
    - `ref_dxy` : (21, N) matrix of reference function xy-derivative values.
    - `ref_dyy` : (21, N) matrix of reference function yy-derivative values.
    """
    check_transformations(C, B)
    check_ref_values(ref_dxx, ref_dxy, ref_dyy)
    dxx = np.empty(ref_dxx.shape, dtype=np.float64)
    dxy = np.empty(ref_dxx.shape, dtype=np.float64)
    dyy = np.empty(ref_dxx.shape, dtype=np.float64)
    _ap.ap_physical_hessians(C, B, ref_dxx, ref_dxy, ref_dyy, ref_dxx.shape[1],
                             dxx, dxy, dyy)
    return (dxx, dxy, dyy)

def matrix_mass(C, B, ref_values, weights):
    """
    Calculate the local mass matrix on a physical triangle.

    Arguments:
    - `C`          : (21, 21) Argyris transformation matrix.
    - `B`          : (2, 2) Affine multiplier matrix.
    - `ref_values` : (21, N) matrix of values of the reference functions at
                     quadrature points.
    - `weights`    : (21,) matrix of the weights corresponding to the
                     quadrature points.
    """
    check_transformations(C, B)
    check_ref_values(ref_values, weights=weights)
    mass = np.empty((21,21), dtype=np.float64)
    _ap.ap_matrix_mass(C, B, ref_values, weights, ref_values.shape[1], mass)
    return mass

def matrix_betaplane(C, B, ref_values, ref_dx, ref_dy, weights):
    """
    Calculate the local betaplane matrix on a physical triangle.

    Arguments:
    - `C`          : (21, 21) Argyris transformation matrix.
    - `B`          : (2, 2) Affine multiplier matrix.
    - `ref_values` : (21, N) matrix of reference function values at quadrature
                     points.
    - `ref_dx`     : (21, N) matrix of reference function x-derivative values at
                     quadrature points.
    - `ref_dy`     : (21, N) matrix of reference function y-derivative values at
                     quadrature points.
    - `weights`    : (21,) matrix of the weights corresponding to the quadrature
                     points.
    """
    check_transformations(C, B)
    check_ref_values(ref_values, ref_dx, ref_dy, weights=weights)
    betaplane = np.empty((21,21), dtype=np.float64)
    _ap.ap_matrix_betaplane(C, B, ref_values, ref_dx, ref_dy, weights,
                            ref_values.shape[1], betaplane)
    return betaplane

def matrix_stiffness(C, B, ref_dx, ref_dy, weights):
    """
    Calculate the local stiffness matrix on a physical triangle.

    Arguments:
    - `C`       : (21, 21) Argyris transformation matrix.
    - `B`       : (2, 2) Affine multiplier matrix.
    - `ref_dx`  : (21, N) matrix of reference function x-derivative values at
                  quadrature points.
    - `ref_dy`  : (21, N) matrix of reference function y-derivative values at
                  quadrature points.
    - `weights` : (21,) matrix of the weights corresponding to the quadrature
                  points.
    """
    check_transformations(C, B)
    check_ref_values(ref_dx, ref_dy, weights=weights)
    stiffness = np.empty((21,21), dtype=np.float64)
    _ap.ap_matrix_stiffness(C, B, ref_dx, ref_dy, weights, ref_dx.shape[1],
                            stiffness)
    return stiffness

def matrix_biharmonic(C, B, ref_dxx, ref_dxy, ref_dyy, weights):
    """
    Calculate the local biharmonic matrix on a physical triangle.

    Arguments:
    - `C`       : (21, 21) Argyris transformation matrix.
    - `B`       : (2, 2) Affine multiplier matrix.
    - `ref_dxx` : (21, N) matrix of reference function xx-derivative values at
                  quadrature points.
    - `ref_dxy` : (21, N) matrix of reference function xy-derivative values at
                  quadrature points.
    - `ref_dyy` : (21, N) matrix of reference function yy-derivative values at
                  quadrature points.
    - `weights` : (21,) matrix of the weights corresponding to the quadrature
                  points.
    """
    check_transformations(C, B)
    check_ref_values(ref_dxx, ref_dxy, ref_dyy, weights=weights)
    biharmonic = np.empty((21,21), dtype=np.float64)
    _ap.ap_matrix_biharmonic(C, B, ref_dxx, ref_dxy, ref_dyy, weights,
                             ref_dxx.shape[1], biharmonic)
    return biharmonic

def check_evaluation_points(x, y):
    """
    Assure that the provided points have the correct shape and type.
    """
    assert x.ndim  == y.ndim == 1
    assert x.shape == y.shape
    assert x.dtype == y.dtype == np.float64

def check_transformations(*args):
    """
    Assure that the C and B transformations have the correct shape and type.
    """
    assert args[0].shape == (21,21)
    assert args[0].dtype == np.float64
    if len(args) == 2:
        assert args[1].shape == (2,2)
        assert args[1].dtype == np.float64

def check_ref_values(*args, **kwargs):
    """
    Assure that the reference values (be them derivatives or function values)
    have the correct shape and type.
    """
    for ref_values in args:
        assert ref_values.shape == (21, args[0].shape[1])
        assert ref_values.dtype == np.float64
    if 'weights' in kwargs.keys():
        assert kwargs['weights'].shape == (args[0].shape[1],)
        assert kwargs['weights'].dtype == np.float64

def get_quad_points():
    """
    Return the tuple (x, y, w) of quadrature data. Taken from Dr. Burkhardt's
    Dunvant program: order 37, degree 13.
    """
    points = np.array(
     [[0.333333333333333333333333333333, 0.333333333333333333333333333333],
      [0.950275662924105565450352089520, 0.024862168537947217274823955239],
      [0.024862168537947217274823955239, 0.950275662924105565450352089520],
      [0.024862168537947217274823955239, 0.024862168537947217274823955239],
      [0.171614914923835347556304795551, 0.414192542538082326221847602214],
      [0.414192542538082326221847602214, 0.171614914923835347556304795551],
      [0.414192542538082326221847602214, 0.414192542538082326221847602214],
      [0.539412243677190440263092985511, 0.230293878161404779868453507244],
      [0.230293878161404779868453507244, 0.539412243677190440263092985511],
      [0.230293878161404779868453507244, 0.230293878161404779868453507244],
      [0.772160036676532561750285570113, 0.113919981661733719124857214943],
      [0.113919981661733719124857214943, 0.772160036676532561750285570113],
      [0.113919981661733719124857214943, 0.113919981661733719124857214943],
      [0.009085399949835353883572964740, 0.495457300025082323058213517632],
      [0.495457300025082323058213517632, 0.009085399949835353883572964740],
      [0.495457300025082323058213517632, 0.495457300025082323058213517632],
      [0.062277290305886993497083640527, 0.468861354847056503251458179727],
      [0.468861354847056503251458179727, 0.062277290305886993497083640527],
      [0.468861354847056503251458179727, 0.468861354847056503251458179727],
      [0.022076289653624405142446876931, 0.851306504174348550389457672223],
      [0.022076289653624405142446876931, 0.126617206172027096933163647918],
      [0.851306504174348550389457672223, 0.022076289653624405142446876931],
      [0.851306504174348550389457672223, 0.126617206172027096933163647918],
      [0.126617206172027096933163647918, 0.022076289653624405142446876931],
      [0.126617206172027096933163647918, 0.851306504174348550389457672223],
      [0.018620522802520968955913511549, 0.689441970728591295496647976487],
      [0.018620522802520968955913511549, 0.291937506468887771754472382212],
      [0.689441970728591295496647976487, 0.018620522802520968955913511549],
      [0.689441970728591295496647976487, 0.291937506468887771754472382212],
      [0.291937506468887771754472382212, 0.018620522802520968955913511549],
      [0.291937506468887771754472382212, 0.689441970728591295496647976487],
      [0.096506481292159228736516560903, 0.635867859433872768286976979827],
      [0.096506481292159228736516560903, 0.267625659273967961282458816185],
      [0.635867859433872768286976979827, 0.096506481292159228736516560903],
      [0.635867859433872768286976979827, 0.267625659273967961282458816185],
      [0.267625659273967961282458816185, 0.096506481292159228736516560903],
      [0.267625659273967961282458816185, 0.635867859433872768286976979827]]);

    w = np.array(
      [0.051739766065744133555179145422,
       0.008007799555564801597804123460,
       0.008007799555564801597804123460,
       0.008007799555564801597804123460,
       0.046868898981821644823226732071,
       0.046868898981821644823226732071,
       0.046868898981821644823226732071,
       0.046590940183976487960361770070,
       0.046590940183976487960361770070,
       0.046590940183976487960361770070,
       0.031016943313796381407646220131,
       0.031016943313796381407646220131,
       0.031016943313796381407646220131,
       0.010791612736631273623178240136,
       0.010791612736631273623178240136,
       0.010791612736631273623178240136,
       0.032195534242431618819414482205,
       0.032195534242431618819414482205,
       0.032195534242431618819414482205,
       0.015445834210701583817692900053,
       0.015445834210701583817692900053,
       0.015445834210701583817692900053,
       0.015445834210701583817692900053,
       0.015445834210701583817692900053,
       0.015445834210701583817692900053,
       0.017822989923178661888748319485,
       0.017822989923178661888748319485,
       0.017822989923178661888748319485,
       0.017822989923178661888748319485,
       0.017822989923178661888748319485,
       0.017822989923178661888748319485,
       0.037038683681384627918546472190,
       0.037038683681384627918546472190,
       0.037038683681384627918546472190,
       0.037038683681384627918546472190,
       0.037038683681384627918546472190,
       0.037038683681384627918546472190])*0.5;
    quad_x = np.copy(points[:,0])
    quad_y = np.copy(points[:,1])
    return (quad_x, quad_y, w)
