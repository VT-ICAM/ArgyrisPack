void ap_global_gradients(double* restrict C, double* restrict B,
                         double* restrict ref_dx, double* restrict ref_dy,
                         LAPACKINDEX num_points,
                         double* restrict dx, double* restrict dy)
{
        double dx_unmapped[21*num_points];
        double dy_unmapped[21*num_points];
        int i;

        /* stuff for DGEMM */
        LAPACKINDEX i_twentyone = 21;

        /* Calculate the global-to-local mapping. */
        const double B_det_inv = 1/(B[ORDER(0, 0, 2, 2)]*B[ORDER(1, 1, 2, 2)] -
                                    B[ORDER(0, 1, 2, 2)]*B[ORDER(1, 0, 2, 2)]);

        const double B_inv00 = B_det_inv*B[ORDER(1, 1, 2, 2)];
        const double B_inv01 = -B_det_inv*B[ORDER(0, 1, 2, 2)];
        const double B_inv10 = -B_det_inv*B[ORDER(1, 0, 2, 2)];
        const double B_inv11 = B_det_inv*B[ORDER(0, 0, 2, 2)];

        /*
         * Perform the transformation using B inverse. This is equivalent to
         * putting the reference values in long columns side-by-side and
         * multiplying by B_inv.
         */
        for (i = 0; i < 21*num_points; i++) {
                dx_unmapped[i] = B_inv00*ref_dx[i] + B_inv10*ref_dy[i];
                dy_unmapped[i] = B_inv01*ref_dx[i] + B_inv11*ref_dy[i];
        }

        /* perform the transformation using the C matrix. */
        DGEMM_WRAPPER(i_twentyone, num_points, i_twentyone, C, dx_unmapped, dx);
        DGEMM_WRAPPER(i_twentyone, num_points, i_twentyone, C, dy_unmapped, dy);
}
