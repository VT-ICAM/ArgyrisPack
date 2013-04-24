void ap_matrix_biharmonic(double* restrict C, double* restrict B,
                          double* restrict Th,
                          double* restrict ref_dxx, double* restrict ref_dxy,
                          double* restrict ref_dyy, double* restrict weights,
                          LAPACKINDEX num_points, double* restrict biharmonic)
{

        int i;
        double dxx[21*num_points];
        double dxy[21*num_points];
        double dyy[21*num_points];
        double weights_scaled[num_points];

        /* stuff for LAPACK */
        LAPACKINDEX i_twentyone = 21;

        const double jacobian = fabs(B[ORDER(0, 0, 2, 2)]*B[ORDER(1, 1, 2, 2)] -
                                     B[ORDER(0, 1, 2, 2)]*B[ORDER(1, 0, 2, 2)]);

        ap_physical_hessians(C, Th, ref_dxx, ref_dxy, ref_dyy, num_points, dxx,
                           dxy, dyy);

        /* Reassign dxx and dyy to be values of the laplacian. */
        for (i = 0; i < 21*num_points; i++) {
                dxx[i] += dyy[i];
        }
        memcpy(dyy, dxx, sizeof(double)*(21*num_points));

        for (i = 0; i < num_points; i++) {
                weights_scaled[i] = weights[i]*jacobian;
        }

        /*
         * Scale one set of laplacian values by the weights (themselves scaled
         * by the Jacobian) and then calculate the matrix of inner products.
         */
        multiply_by_diagonal(21, num_points, weights_scaled, dxx);
        DGEMM_WRAPPER_NT(i_twentyone, i_twentyone, num_points, dxx, dyy,
                         biharmonic);
}
