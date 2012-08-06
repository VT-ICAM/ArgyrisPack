void ap_local_functions(double* restrict x, double* restrict y, int quad_points,
                        double* restrict ref_functions)
{
        double monomials[21*quad_points];
        int i;

        /* stuff for dgemm */
        int i_twentyone = 21;

#include "function_coefficients.h"

        /*
         * Rows in the monomial matrix correspond to monomials (x, y,
         * x^2, etc) while columns correspond to quadrature points. The
         * monomial basis spans quintic polynomials of dimension 2
         * (hence 21 rows).
         */
        for (i = 0; i < quad_points; i++) {
                monomials[ORDER(0,i,21,quad_points)]  = 1.0;
                monomials[ORDER(1,i,21,quad_points)]  = x[i];
                monomials[ORDER(2,i,21,quad_points)]  = y[i];
                monomials[ORDER(3,i,21,quad_points)]  = x[i]*x[i];
                monomials[ORDER(4,i,21,quad_points)]  = x[i]*y[i];
                monomials[ORDER(5,i,21,quad_points)]  = y[i]*y[i];
                monomials[ORDER(6,i,21,quad_points)]  = x[i]*x[i]*x[i];
                monomials[ORDER(7,i,21,quad_points)]  = x[i]*x[i]*y[i];
                monomials[ORDER(8,i,21,quad_points)]  = x[i]*y[i]*y[i];
                monomials[ORDER(9,i,21,quad_points)]  = y[i]*y[i]*y[i];
                monomials[ORDER(10,i,21,quad_points)] = x[i]*x[i]*x[i]*x[i];
                monomials[ORDER(11,i,21,quad_points)] = x[i]*x[i]*x[i]*y[i];
                monomials[ORDER(12,i,21,quad_points)] = x[i]*x[i]*y[i]*y[i];
                monomials[ORDER(13,i,21,quad_points)] = x[i]*y[i]*y[i]*y[i];
                monomials[ORDER(14,i,21,quad_points)] = y[i]*y[i]*y[i]*y[i];
                monomials[ORDER(15,i,21,quad_points)] = x[i]*x[i]*x[i]*x[i]*x[i];
                monomials[ORDER(16,i,21,quad_points)] = x[i]*x[i]*x[i]*x[i]*y[i];
                monomials[ORDER(17,i,21,quad_points)] = x[i]*x[i]*x[i]*y[i]*y[i];
                monomials[ORDER(18,i,21,quad_points)] = x[i]*x[i]*y[i]*y[i]*y[i];
                monomials[ORDER(19,i,21,quad_points)] = x[i]*y[i]*y[i]*y[i]*y[i];
                monomials[ORDER(20,i,21,quad_points)] = y[i]*y[i]*y[i]*y[i]*y[i];
        }

        DGEMM_WRAPPER(i_twentyone, quad_points, i_twentyone, coefficients,
                      monomials, ref_functions);
}
