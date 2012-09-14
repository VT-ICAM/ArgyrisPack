import numpy as np
import operator as op

def argyris_polynomials(xs, ys):
    """
    Calculate (symbolically) the Argyris basis functions for some given
    triangle.
    """
    var('x,y')

    edge_vector = {
        1 : matrix([[xs[1] - xs[0]], [ys[1] - ys[0]]]),
        2 : matrix([[xs[2] - xs[0]], [ys[2] - ys[0]]]),
        3 : matrix([[xs[2] - xs[1]], [ys[2] - ys[1]]])
    }

    edge_normal = {i : matrix([[0,-1],[1,0]])*
                   edge_vector[i]/sqrt(edge_vector[i][0,0]^2 + edge_vector[i][1,0]^2)
                   for i in edge_vector.keys()}

    edge_midpoint = {
        1 : matrix([[(xs[0] + xs[1])/2], [(ys[0] + ys[1])/2]]),
        2 : matrix([[(xs[0] + xs[2])/2], [(ys[0] + ys[2])/2]]),
        3 : matrix([[(xs[1] + xs[2])/2], [(ys[1] + ys[2])/2]])
    }

    # I kept this as a function for debugging purposes.
    constraint_system = lambda f : [
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

    # set up the system of linear equations resulting from applying the
    # constraints to some arbitrary order 5 polynomial.
    constants = [var('c' + str(num)) for num in [1..21]]
    monic_basis = [x^m*y^(n-m) for n in [0..5] for m in [n-k for k in [0..n]]]
    test_polynomial = sum(map(op.mul, constants, monic_basis))

    constraints = constraint_system(test_polynomial)
    constraint_matrix = matrix([get_coefficients(constraint, constants)
                                 for constraint in constraints])

    constraint_inverse = constraint_matrix.inverse()
    coefficients = []

    for i in range(1,22):
        rhs = matrix(SR, 21, 1)
        rhs[i-1,0] = 1
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
