void ap_diagonal_multiply(const int rows, const int cols,
                          double* restrict matrix, double* restrict diagonal) {
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
