void ap_matrix_stiffness(double* restrict C, double* restrict B,
                         double* restrict ref_dx, double* restrict ref_dy,
                         double* restrict weights,
                         LAPACKINDEX num_points, double* restrict stiffness)
{

        int i;
        double dx[21*num_points];
        double dx_scaled[21*num_points];
        double dy[21*num_points];
        double dy_scaled[21*num_points];
        double weights_scaled[num_points];

        /* stuff for DGEMM. */
        LAPACKINDEX i_twentyone = 21;

        const double jacobian = fabs(B[ORDER(0, 0, 2, 2)]*B[ORDER(1, 1, 2, 2)] -
                                     B[ORDER(0, 1, 2, 2)]*B[ORDER(1, 0, 2, 2)]);

        ap_global_gradients(C, B, ref_dx, ref_dy, num_points, dx, dy);

        /* scale the weights by the jacobian. */
        for (i = 0; i < num_points; i++) {
                weights_scaled[i] = weights[i]*jacobian;
        }

        memcpy(dx_scaled, dx, sizeof(double)*(21*num_points));
        memcpy(dy_scaled, dy, sizeof(double)*(21*num_points));

        /*
         * scale the first set of gradient values by the weights and
         * determinant. Then perform matrix multiplication.
         */
        multiply_by_diagonal(21, num_points, weights_scaled, dx_scaled);
        multiply_by_diagonal(21, num_points, weights_scaled, dy_scaled);

        DGEMM_WRAPPER_NT(i_twentyone, i_twentyone, num_points, dx_scaled, dx,
                         stiffness);
        DGEMM_WRAPPER_NT_ADD_C(i_twentyone, i_twentyone, num_points, dy_scaled,
                               dy, stiffness);
}
