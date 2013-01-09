void ap_matrix_betaplane(double* restrict C, double* restrict B,
                         double* restrict ref_values,
                         double* restrict ref_dx, double* restrict ref_dy,
                         double* restrict weights,
                         LAPACKINDEX quad_points, double* restrict betaplane)
{

        int i;
        double values[21*quad_points];
        double dx[21*quad_points];
        double dy[21*quad_points];
        double weights_scaled[quad_points];

        /* stuff for DGEMM. */
        LAPACKINDEX i_twentyone = 21;

        const double jacobian = fabs(B[ORDER(0, 0, 2, 2)]*B[ORDER(1, 1, 2, 2)] -
                                     B[ORDER(0, 1, 2, 2)]*B[ORDER(1, 0, 2, 2)]);

        ap_global_gradients(C, B, ref_dx, ref_dy, quad_points, dx, dy);
        ap_global_functions(C, ref_values, quad_points, values);

        /* scale the weights by the jacobian. */
        for (i = 0; i < quad_points; i++) {
                weights_scaled[i] = weights[i]*jacobian;
        }

        /*
         * scale the function values by the weights and determinant. Then
         * perform matrix multiplication.
         */
        multiply_by_diagonal(21, quad_points, weights_scaled, values);

        DGEMM_WRAPPER_NT(i_twentyone, i_twentyone, quad_points, values, dx,
                         betaplane);
}
