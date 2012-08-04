#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "global_hessians.c"

void multiply_by_diagonal(const int rows, const int cols,
                          double* restrict diagonal, double* restrict matrix);

void ap_matrix_biharmonic(double* restrict C, double* restrict B,
                          double* restrict Th,
                          double* restrict ref_dxx, double* restrict ref_dxy,
                          double* restrict ref_dyy, double* restrict weights,
                          LAPACKINDEX quad_points, double* restrict biharmonic)
{

        int i;
        double dxx[21*quad_points];
        double dxy[21*quad_points];
        double dyy[21*quad_points];
        double weights_scaled[quad_points];

        /* stuff for LAPACK */
        LAPACKINDEX i_twentyone = 21;

        const double jacobian = fabs(B[ORDER(0, 0, 2, 2)]*B[ORDER(1, 1, 2, 2)] -
                                     B[ORDER(0, 1, 2, 2)]*B[ORDER(1, 0, 2, 2)]);

        ap_global_hessians(C, Th, ref_dxx, ref_dxy, ref_dyy, quad_points, dxx,
                           dxy, dyy);

        /* Reassign dxx and dyy to be values of the laplacian. */
        for (i = 0; i < 21*quad_points; i++) {
                dxx[i] += dyy[i];
        }
        memcpy(dyy, dxx, sizeof(double)*(21*quad_points));

        for (i = 0; i < quad_points; i++) {
                weights_scaled[i] = weights[i]*jacobian;
        }

        /*
         * Scale one set of laplacian values by the weights (themselves scaled
         * by the Jacobian) and then calculate the matrix of inner products.
         */
        multiply_by_diagonal(21, quad_points, weights_scaled, dxx);
        DGEMM_WRAPPER_NT(i_twentyone, i_twentyone, quad_points, dxx, dyy,
                         biharmonic);
}

void multiply_by_diagonal(const int rows, const int cols,
                          double* restrict diagonal, double* restrict matrix) {
        /*
         * Multiply a matrix by a diagonal matrix (represented as a flat array
         * of n values). The diagonal structure is exploited (so that we use
         * n^2 multiplications instead of n^3); the product looks something
         * like
         *
         *     sage: A = matrix([[a,b,c], [d,e,f],[g,h,i]]);
         *     sage: w = matrix([[z1,0,0],[0,z2,0],[0,0,z3]]);
         *     sage: A * w
         *      => [a*z1 b*z2 c*z3]
         *         [d*z1 e*z2 f*z3]
         *         [g*z1 h*z2 i*z3]
         */
        int i, j;
        /* traverse the matrix in the correct order. */
#ifdef USE_COL_MAJOR
        for (j = 0; j < cols; j++)
                for (i = 0; i < rows; i++)
#else
        for (i = 0; i < rows; i++)
                for (j = 0; j < cols; j++)
#endif
                        matrix[ORDER(i, j, rows, cols)] *= diagonal[j];
}
