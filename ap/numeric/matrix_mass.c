void ap_matrix_mass(double* restrict C, double* restrict B,
                    double* restrict ref_values, double* restrict weights,
                    LAPACKINDEX num_points, double* restrict mass)
{
        int i;
        double function_values[21*num_points];
        double function_values_scaled[21*num_points];
        double weights_scaled[num_points];

        /* stuff for DGEMM. */
        LAPACKINDEX i_twentyone = 21;

        const double jacobian = fabs(B[ORDER(0, 0, 2, 2)]*B[ORDER(1, 1, 2, 2)] -
                                     B[ORDER(0, 1, 2, 2)]*B[ORDER(1, 0, 2, 2)]);

        ap_physical_values(C, ref_values, num_points, function_values);

        /* scale the weights by the jacobian. */
        for (i = 0; i < num_points; i++) {
                weights_scaled[i] = weights[i]*jacobian;
        }

        memcpy(function_values_scaled, function_values,
               sizeof(double)*(21*num_points));

        /*
         * scale the first set of function values by the weights and
         * determinant. Then perform matrix multiplication.
         */
        ap_diagonal_multiply(21, num_points, function_values_scaled, weights_scaled);
        DGEMM_WRAPPER_NT(i_twentyone, i_twentyone, num_points,
                         function_values_scaled, function_values, mass);
}
