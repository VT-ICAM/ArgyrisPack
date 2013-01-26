#! /usr/bin/env python

# Compute the matricies of inner products. Each function takes a
# InterpolatedBasisFunctions as an argument and returns the appropriate matrix
# of inner products.
import operator as op
import numpy as np
from sage.all import *
import ap.symbolic.symbolic as symb

def integrate_simplex(f):
    var('x, y')
    return float(integrate(integrate(f, y, 0, 1 - x), x, 0, 1))

def mass(basis_functions):
    return np.array([[basis_functions.jacobian*
                      integrate_simplex(basis_functions.ref_polynomials[i]*
                                        basis_functions.ref_polynomials[j])
                                        for j in range(21)]
                     for i in range(21)])

def beta_plane(basis_functions):
    return np.array([[basis_functions.jacobian*
                      integrate_simplex(basis_functions.ref_polynomials[i]*
                                        basis_functions.ref_polynomials_x[j])
                                        for j in range(21)]
                     for i in range(21)])

def stiffness(basis_functions):
    return np.array([[basis_functions.jacobian*
                      integrate_simplex(basis_functions.ref_polynomials_x[i]*
                                        basis_functions.ref_polynomials_x[j]+
                                        basis_functions.ref_polynomials_y[i]*
                                        basis_functions.ref_polynomials_y[j])
                                        for j in range(21)]
                     for i in range(21)])

def biharmonic(basis_functions):
    return np.array([[basis_functions.jacobian*
                      integrate_simplex((basis_functions.ref_polynomials_xx[i]+
                                         basis_functions.ref_polynomials_yy[i])*
                                        (basis_functions.ref_polynomials_xx[j]+
                                         basis_functions.ref_polynomials_yy[j]))
                                         for j in range(21)]
                     for i in range(21)])

def mass_stabilized(basis_functions):
    return np.array([[basis_functions.jacobian*
                      integrate_simplex(basis_functions.ref_polynomials_x[i]*
                                        basis_functions.ref_polynomials[j])
                                        for j in range(21)]
                     for i in range(21)])


def beta_plane_stabilized(basis_functions):
    return np.array([[basis_functions.jacobian*
                      integrate_simplex(basis_functions.ref_polynomials_x[i]*
                                        basis_functions.ref_polynomials_x[j])
                                        for j in range(21)]
                     for i in range(21)])

def stiffness_stabilized(basis_functions):
    return np.array([[basis_functions.jacobian*
                      integrate_simplex(basis_functions.ref_polynomials_xx[i]*
                                        basis_functions.ref_polynomials_x[j]+
                                        basis_functions.ref_polynomials_xy[i]*
                                        basis_functions.ref_polynomials_y[j])
                                        for j in range(21)]
                     for i in range(21)])

def biharmonic_stabilized(basis_functions):
    return np.array([[basis_functions.jacobian*
                      integrate_simplex((basis_functions.ref_polynomials_xxx[i]+
                                         basis_functions.ref_polynomials_xyy[i])*
                                        (basis_functions.ref_polynomials_xx[j]+
                                         basis_functions.ref_polynomials_yy[j]))
                                         for j in range(21)]
                     for i in range(21)])
