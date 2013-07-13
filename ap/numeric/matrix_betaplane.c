void ap_matrix_betaplane(double* restrict C, double* restrict B,
                         double* restrict ref_values,
                         double* restrict ref_dx, double* restrict ref_dy,
                         double* restrict weights,
                         LAPACKINDEX num_points, double* restrict betaplane)
{
        int i;
        double values[21*num_points];
        double dx[21*num_points];
        double dy[21*num_points];
        double weights_scaled[num_points];

        /* stuff for DGEMM. */
        LAPACKINDEX i_twentyone = 21;

        const double jacobian = fabs(B[ORDER(0, 0, 2, 2)]*B[ORDER(1, 1, 2, 2)] -
                                     B[ORDER(0, 1, 2, 2)]*B[ORDER(1, 0, 2, 2)]);

        ap_physical_gradients(C, B, ref_dx, ref_dy, num_points, dx, dy);
        ap_physical_values(C, ref_values, num_points, values);

        /* scale the weights by the jacobian. */
        for (i = 0; i < num_points; i++) {
                weights_scaled[i] = weights[i]*jacobian;
        }

        /*
         * scale the function values by the weights and determinant. Then
         * perform matrix multiplication.
         */
        ap_diagonal_multiply(21, num_points, values, weights_scaled);

        DGEMM_WRAPPER_NT(i_twentyone, i_twentyone, num_points, values, dx,
                         betaplane);
}
