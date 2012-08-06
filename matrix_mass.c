#include <string.h>
#include <stdlib.h>

#include "global_functions.c"

void multiply_by_diagonal(const int rows, const int cols,
                          double* restrict diagonal, double* restrict matrix);

void ap_matrix_mass(double* restrict C, double* restrict B,
                    double* restrict ref_functions, double* restrict weights,
                    LAPACKINDEX quad_points, double* restrict mass)
{

        int i;
        double functions[21*quad_points];
        double functions_scaled[21*quad_points];
        double weights_scaled[quad_points];

        /* stuff for DGEMM. */
        LAPACKINDEX i_twentyone = 21;

        const double jacobian = fabs(B[ORDER(0, 0, 2, 2)]*B[ORDER(1, 1, 2, 2)] -
                                     B[ORDER(0, 1, 2, 2)]*B[ORDER(1, 0, 2, 2)]);

        ap_global_functions(C, ref_functions, quad_points, functions);

        /* scale the weights by the jacobian. */
        for (i = 0; i < quad_points; i++) {
                weights_scaled[i] = weights[i]*jacobian;
        }

        memcpy(functions_scaled, functions, sizeof(double)*(21*quad_points));

        /*
         * scale the first set of function values by the weights and
         * determinant. Then perform matrix multiplication.
         */
        multiply_by_diagonal(21, quad_points, weights_scaled, functions_scaled);
        DGEMM_WRAPPER_NT(i_twentyone, i_twentyone, quad_points,
                         functions_scaled, functions, mass);
}

void multiply_by_diagonal(const int rows, const int cols,
                          double* restrict diagonal, double* restrict matrix) {
        /*
         * Multiply a matrix by a diagonal matrix (represented as a flat array
         * of n values). The diagonal structure is exploited; the product looks
         * something like
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
