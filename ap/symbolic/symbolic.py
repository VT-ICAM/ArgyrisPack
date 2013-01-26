#! /usr/bin/env python

import operator as op
import numpy as np
from sage.all import *

class InterpolatedBasisFunctions(SageObject):
    """
    Class for the
    """
    def __init__(self, xs, ys):
        var('x, y')
        xs = map(Rational, xs)
        ys = map(Rational, ys)

        self.corners = zip(xs, ys)
        self.polynomials = argyris_polynomials(xs, ys)
        self.ref_polynomials = []
        self.ref_polynomials_x = []
        self.ref_polynomials_y = []
        self.ref_polynomials_xx = []
        self.ref_polynomials_xy = []
        self.ref_polynomials_yy = []
        self.ref_polynomials_xxx = []
        self.ref_polynomials_xxy = []
        self.ref_polynomials_xyy = []

        # Compute the change of coordinates from the specified triangle to the
        # reference simplex.
        local_to_global = matrix([[xs[1] - xs[0], xs[2] - xs[0]],
                                  [ys[1] - ys[0], ys[2] - ys[0]]])
        affine_shift = matrix([[xs[0]], [ys[0]]])
        ref_values = local_to_global*matrix([[x],[y]]) + affine_shift
        xref = ref_values[0,0]
        yref = ref_values[1,0]
        self.jacobian = float(local_to_global.det())

        for polynomial in self.polynomials:
            self.ref_polynomials.append(
                polynomial.subs({x : xref, y : yref}))
            self.ref_polynomials_x.append(
                diff(polynomial,x).subs({x : xref, y : yref}))
            self.ref_polynomials_y.append(
                diff(polynomial,y).subs({x : xref, y : yref}))
            self.ref_polynomials_xx.append(
                diff(polynomial,x,2).subs({x : xref, y : yref}))
            self.ref_polynomials_xy.append(
                diff(diff(polynomial,x),y).subs({x : xref, y : yref}))
            self.ref_polynomials_yy.append(
                diff(polynomial,y,2).subs({x : xref, y : yref}))
            self.ref_polynomials_xxx.append(
                diff(polynomial,x,3).subs({x : xref, y : yref}))
            self.ref_polynomials_xxy.append(
                diff(diff(polynomial,x,2),y).subs({x : xref, y : yref}))
            self.ref_polynomials_xyy.append(
                diff(diff(polynomial,y,2),x).subs({x : xref, y : yref}))

def argyris_polynomials(xs, ys):
    """
    Calculate (symbolically) the Argyris basis functions for some given
    triangle.
    """
    var('x,y')

    constants = [var('c' + str(num)) for num in range(1,22)]
    monic_basis = [x**m * y**(n-m) for n in range(6)
                   for m in [n-k for k in range(n+1)]]
    test_polynomial = sum(map(op.mul, constants, monic_basis))

    constraints = constraint_system(test_polynomial, xs, ys)
    constraint_matrix = matrix([get_coefficients(constraint, constants)
                                 for constraint in constraints])

    constraint_inverse = constraint_matrix.inverse()
    coefficients = []

    for i in range(21):
        rhs = matrix(SR, 21, 1)
        rhs[i, 0] = 1
        coefficients.append(constraint_inverse*rhs)

    # form the finite element polynomials by recombining the coefficients with
    # the monic basis.
    return [sum(map(op.mul, c, monic_basis))[0] for c in coefficients]


def get_coefficients(linear_equation, variables):
    """
    Extract the constants associated with variables in some linear equation.
    """
    coefficients = []
    for variable in variables:
        coefficients.append(linear_equation.coefficient(variable))

    return coefficients

def constraint_system(f, xs, ys):
    var('x, y')

    edge_vector = {
        1 : matrix([[xs[1] - xs[0]], [ys[1] - ys[0]]]),
        2 : matrix([[xs[2] - xs[0]], [ys[2] - ys[0]]]),
        3 : matrix([[xs[2] - xs[1]], [ys[2] - ys[1]]])
    }

    edge_midpoint = {
        1 : matrix([[(xs[0] + xs[1])/2], [(ys[0] + ys[1])/2]]),
        2 : matrix([[(xs[0] + xs[2])/2], [(ys[0] + ys[2])/2]]),
        3 : matrix([[(xs[1] + xs[2])/2], [(ys[1] + ys[2])/2]])
    }

    edge_normal = {i : matrix([[0,-1],[1,0]])*
                   edge_vector[i]/sqrt(edge_vector[i][0,0]**2 + edge_vector[i][1,0]**2)
                   for i in edge_vector.keys()}

    return [
    f.subs({x : xs[0], y : ys[0]}),
    f.subs({x : xs[1], y : ys[1]}),
    f.subs({x : xs[2], y : ys[2]}),
    diff(f,x).subs({x : xs[0], y : ys[0]}),
    diff(f,y).subs({x : xs[0], y : ys[0]}),
    diff(f,x).subs({x : xs[1], y : ys[1]}),
    diff(f,y).subs({x : xs[1], y : ys[1]}),
    diff(f,x).subs({x : xs[2], y : ys[2]}),
    diff(f,y).subs({x : xs[2], y : ys[2]}),
    diff(f,x,2).subs({x : xs[0], y : ys[0]}),
    diff(diff(f,x),y).subs({x : xs[0], y : ys[0]}),
    diff(f,y,2).subs({x : xs[0], y : ys[0]}),
    diff(f,x,2).subs({x : xs[1], y : ys[1]}),
    diff(diff(f,x),y).subs({x : xs[1], y : ys[1]}),
    diff(f,y,2).subs({x : xs[1], y : ys[1]}),
    diff(f,x,2).subs({x : xs[2], y : ys[2]}),
    diff(diff(f,x),y).subs({x : xs[2], y : ys[2]}),
    diff(f,y,2).subs({x : xs[2], y : ys[2]}),
    ((diff(f,x).subs({x : edge_midpoint[1][0,0], y : edge_midpoint[1][1,0]})
      *edge_normal[1][0,0]) +
    (diff(f,y).subs({x : edge_midpoint[1][0,0], y : edge_midpoint[1][1,0]})
      *edge_normal[1][1,0])),
    ((diff(f,x).subs({x : edge_midpoint[2][0,0], y : edge_midpoint[2][1,0]})
      *edge_normal[2][0,0]) +
     (diff(f,y).subs({x : edge_midpoint[2][0,0], y : edge_midpoint[2][1,0]})
      *edge_normal[2][1,0])),
    ((diff(f,x).subs({x : edge_midpoint[3][0,0], y : edge_midpoint[3][1,0]})
        *edge_normal[3][0,0]) +
    (diff(f,y).subs({x : edge_midpoint[3][0,0], y : edge_midpoint[3][1,0]})
        *edge_normal[3][1,0]))]
