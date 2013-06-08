void ap_physical_hessians(double* restrict C, double* restrict B,
                          double* restrict ref_dxx, double* restrict ref_dxy,
                          double* restrict ref_dyy, LAPACKINDEX num_points,
                          double* restrict dxx, double* restrict dxy,
                          double* restrict dyy)
{
        double dxx_unmapped[21*num_points];
        double dxy_unmapped[21*num_points];
        double dyy_unmapped[21*num_points];
        int i;

        /* stuff for DGEMM */
        LAPACKINDEX i_twentyone = 21;

        /*
         * There is an extra 3x3 matrix (corresponding, in Dominguez's notation,
         * to Theta transpose) due to application of the chain rule. It's
         * entries depend on the original affine transformation from reference
         * to physical coordinates (hence the dependence on B).
         */
        const double t = (B[ORDER(0, 0, 2, 2)]*B[ORDER(1, 1, 2, 2)]
                        - B[ORDER(0, 1, 2, 2)]*B[ORDER(1, 0, 2, 2)])*
                         (B[ORDER(0, 0, 2, 2)]*B[ORDER(1, 1, 2, 2)]
                        - B[ORDER(0, 1, 2, 2)]*B[ORDER(1, 0, 2, 2)]);

        const double map00 = B[ORDER(1, 1, 2, 2)]*B[ORDER(1, 1, 2, 2)]/t;
        const double map01 = -2.0*B[ORDER(1, 0, 2, 2)]*B[ORDER(1, 1, 2, 2)]/t;
        const double map02 = B[ORDER(1, 0, 2, 2)]*B[ORDER(1, 0, 2, 2)]/t;

        const double map10 = -1.0*B[ORDER(0, 1, 2, 2)]*B[ORDER(1, 1, 2, 2)]/t;
        const double map11 = (B[ORDER(0, 0, 2, 2)]*B[ORDER(1, 1, 2, 2)]
                            + B[ORDER(0, 1, 2, 2)]*B[ORDER(1, 0, 2, 2)])/t;
        const double map12 = -1.0*B[ORDER(0, 0, 2, 2)]*B[ORDER(1, 0, 2, 2)]/t;

        const double map20 = B[ORDER(0, 1, 2, 2)]*B[ORDER(0, 1, 2, 2)]/t;
        const double map21 = -2.0*B[ORDER(0, 0, 2, 2)]*B[ORDER(0, 1, 2, 2)]/t;
        const double map22 = B[ORDER(0, 0, 2, 2)]*B[ORDER(0, 0, 2, 2)]/t;

        /*
         * Perform the transformation. This is equivalent to putting the
         * reference values in long columns side-by-side and multiplying by
         * (Theta inverse) transpose.
         */
        for (i = 0; i < i_twentyone*num_points; i++) {
                dxx_unmapped[i] = ref_dxx[i]*map00
                                + ref_dxy[i]*map01
                                + ref_dyy[i]*map02;
                dxy_unmapped[i] = ref_dxx[i]*map10
                                + ref_dxy[i]*map11
                                + ref_dyy[i]*map12;
                dyy_unmapped[i] = ref_dxx[i]*map20
                                + ref_dxy[i]*map21
                                + ref_dyy[i]*map22;
        }

        /* perform the transformation using the C matrix. */
        DGEMM_WRAPPER(i_twentyone, num_points, i_twentyone, C, dxx_unmapped,
                      dxx);
        DGEMM_WRAPPER(i_twentyone, num_points, i_twentyone, C, dxy_unmapped,
                      dxy);
        DGEMM_WRAPPER(i_twentyone, num_points, i_twentyone, C, dyy_unmapped,
                      dyy);
}
