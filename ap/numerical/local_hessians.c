void ap_local_hessians(double* restrict x, double* restrict y, LAPACKINDEX num_points,
                       double* restrict ref_dxx, double* restrict ref_dxy,
                       double* restrict ref_dyy)
{
        double monomials[10*num_points];
        int i;

        /* stuff for dgemm */
        LAPACKINDEX i_twentyone = 21;
        LAPACKINDEX i_ten = 10;

#include "coefficients_hessians.h"

        /*
         * Rows in the monomial matrix correspond to monomials (x, y, x^2, etc)
         * while columns correspond to quadrature points. The monomial basis
         * spans cubic polynomials of dimension 2 (hence 10 rows).
         */
        for (i = 0; i < num_points; i++) {
            monomials[ORDER(0, i, 10, num_points)] = 1.0;
            monomials[ORDER(1, i, 10, num_points)] = x[i];
            monomials[ORDER(2, i, 10, num_points)] = y[i];
            monomials[ORDER(3, i, 10, num_points)] = x[i]*x[i];
            monomials[ORDER(4, i, 10, num_points)] = x[i]*y[i];
            monomials[ORDER(5, i, 10, num_points)] = y[i]*y[i];
            monomials[ORDER(6, i, 10, num_points)] = x[i]*x[i]*x[i];
            monomials[ORDER(7, i, 10, num_points)] = x[i]*x[i]*y[i];
            monomials[ORDER(8, i, 10, num_points)] = x[i]*y[i]*y[i];
            monomials[ORDER(9, i, 10, num_points)] = y[i]*y[i]*y[i];
        }

        DGEMM_WRAPPER(i_twentyone, num_points, i_ten, coefficients_dxx,
                      monomials, ref_dxx);
        DGEMM_WRAPPER(i_twentyone, num_points, i_ten, coefficients_dxy,
                      monomials, ref_dxy);
        DGEMM_WRAPPER(i_twentyone, num_points, i_ten, coefficients_dyy,
                      monomials, ref_dyy);
}
