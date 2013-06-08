#! /usr/bin/env python
import itertools as it
import numpy as np
from sage.all import *
import ap.symbolic.symbolic as symb
import ap.symbolic.inner_products as ip
import ap.numeric as nm

def test_triangle(xs, ys):
    """
    Test the numeric implementation against the symbolic one on a given
    triangle.

    Arguments:
    - `xs`: x-coordinates of vertices.
    - `ys`: y-coordinates of vertices.
    """
    if len(xs) != 3 or len(ys) != 3:
        raise ValueError("List of corners must have three entries.")

    symbolic_interpolation = symb.InterpolatedBasisFunctions(xs, ys)

    # calculate all the numeric quantities.
    (quad_x, quad_y, quad_weights) = nm.get_quad_points()
    ref_values = nm.ref_values(quad_x, quad_y)
    (ref_dx, ref_dy) = nm.ref_gradients(quad_x, quad_y)
    (ref_dxx, ref_dxy, ref_dyy) = nm.ref_hessians(quad_x, quad_y)
    (C, B, b) = nm.physical_maps(xs, ys)
    physical_values = nm.physical_values(C, ref_values)
    (dx, dy) = nm.physical_gradients(C, B, ref_dx, ref_dy)
    (dxx, dxy, dyy) = nm.physical_hessians(C, B, ref_dxx, ref_dxy, ref_dyy)

    mass_symbolic = ip.mass(symbolic_interpolation)
    mass_numeric  = nm.matrix_mass(C, B, ref_values, quad_weights)
    print "max relative error for mass matrix:", \
        relative_error(mass_symbolic, mass_numeric)

    stiffness_symbolic = ip.stiffness(symbolic_interpolation)
    stiffness_numeric  = nm.matrix_stiffness(C, B, ref_dx, ref_dy, quad_weights)
    print "max relative error for stiffness matrix:", \
        relative_error(stiffness_symbolic, stiffness_numeric)

    betaplane_symbolic = ip.betaplane(symbolic_interpolation)
    betaplane_numeric  = nm.matrix_betaplane(C, B, ref_values, ref_dx, ref_dy, quad_weights)
    print "max relative error for betaplane matrix:", \
        relative_error(betaplane_symbolic, betaplane_numeric)

    biharmonic_symbolic = ip.biharmonic(symbolic_interpolation)
    biharmonic_numeric  = nm.matrix_biharmonic(C, B, ref_dxx, ref_dxy, ref_dyy,
                                               quad_weights)
    print "max relative error for biharmonic matrix:", \
        relative_error(biharmonic_symbolic, biharmonic_numeric)

    symbolic_values = np.zeros(ref_values.shape, dtype=np.float64)*np.nan
    symbolic_dx     = np.zeros(ref_values.shape, dtype=np.float64)*np.nan
    symbolic_dy     = np.zeros(ref_values.shape, dtype=np.float64)*np.nan
    symbolic_dxx    = np.zeros(ref_values.shape, dtype=np.float64)*np.nan
    symbolic_dxy    = np.zeros(ref_values.shape, dtype=np.float64)*np.nan
    symbolic_dyy    = np.zeros(ref_values.shape, dtype=np.float64)*np.nan

    var('x, y')
    ref_to_physical = np.array([[xs[1] - xs[0], xs[2] - xs[0]],
                                [ys[1] - ys[0], ys[2] - ys[0]]])
    affine_shift = np.array([xs[0], ys[0]])
    map_quad    = np.dot(ref_to_physical, np.vstack((quad_x, quad_y)))
    map_quad_x = map_quad[0,:] + affine_shift[0]
    map_quad_y = map_quad[1,:] + affine_shift[1]

    for i in range(ref_values.shape[0]):
        for j in range(ref_values.shape[1]):
            symbolic_values[i,j] = \
            symbolic_interpolation.physical_polynomials[i].subs({x : map_quad_x[j],
                                                                 y : map_quad_y[j]})
            symbolic_dx[i,j] = \
            symbolic_interpolation.physical_polynomials_x[i].subs({x : map_quad_x[j],
                                                                   y : map_quad_y[j]})
            symbolic_dy[i,j] = \
            symbolic_interpolation.physical_polynomials_y[i].subs({x : map_quad_x[j],
                                                                   y : map_quad_y[j]})
            symbolic_dxx[i,j] = \
            symbolic_interpolation.physical_polynomials_xx[i].subs({x : map_quad_x[j],
                                                                    y : map_quad_y[j]})
            symbolic_dxx[i,j] = \
            symbolic_interpolation.physical_polynomials_xx[i].subs({x : map_quad_x[j],
                                                                    y : map_quad_y[j]})
            symbolic_dxy[i,j] = \
            symbolic_interpolation.physical_polynomials_xy[i].subs({x : map_quad_x[j],
                                                                    y : map_quad_y[j]})
            symbolic_dyy[i,j] = \
            symbolic_interpolation.physical_polynomials_yy[i].subs({x : map_quad_x[j],
                                                                    y : map_quad_y[j]})

    print "max relative error for values:", relative_error(symbolic_values,
                                                               physical_values)
    print "max relative error for dx:",  relative_error(symbolic_dx, dx)
    print "max relative error for dy:",  relative_error(symbolic_dy, dy)
    print "max relative error for dxx:", relative_error(symbolic_dxx, dxx)
    print "max relative error for dxy:", relative_error(symbolic_dxy, dxy)
    print "max relative error for dyy:", relative_error(symbolic_dyy, dyy)

def relative_error(true, approx):
    """
    Determine the error between two matricies by computing the 2-norm of the
    difference and dividing by the 2-norm of the 'true' value.
    """
    return np.linalg.norm(true - approx, ord=2)/np.linalg.norm(true, ord=2)

def run():
    """
    Run enough tests to verify the numeric part of the code.
    """
    test_triangle(np.array([1.0, 2.0, 1.5]), np.array([1.0, 1.0, 2.0]))
    test_triangle(np.array([4.0, 6.0, 5.1]), np.array([2.0, 2.0, 1.0]))
