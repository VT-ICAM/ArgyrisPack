void ap_global_hessians(double* restrict C, double* restrict Th,
                        double* restrict ref_dxx, double* restrict ref_dxy,
                        double* restrict ref_dyy, LAPACKINDEX quad_points,
                        double* restrict dxx, double* restrict dxy,
                        double* restrict dyy)
{
        double dxx_unmapped[21*quad_points];
        double dxy_unmapped[21*quad_points];
        double dyy_unmapped[21*quad_points];
        int i;

        /* stuff for DGEMM */
        LAPACKINDEX i_twentyone = 21;

        /* Calculate Th_inv entries. */
        const double Th_det =
                Th[ORDER(0, 0, 3, 3)]*(Th[ORDER(1, 1, 3, 3)]*Th[ORDER(2, 2, 3, 3)] -
                                    Th[ORDER(1, 2, 3, 3)]*Th[ORDER(2, 1, 3, 3)]) +
                Th[ORDER(0, 1, 3, 3)]*(Th[ORDER(1, 2, 3, 3)]*Th[ORDER(2, 0, 3, 3)] -
                                    Th[ORDER(2, 2, 3, 3)]*Th[ORDER(1, 0, 3, 3)]) +
                Th[ORDER(0, 2, 3, 3)]*(Th[ORDER(1, 0, 3, 3)]*Th[ORDER(2, 1, 3, 3)] -
                                    Th[ORDER(1, 1, 3, 3)]*Th[ORDER(2, 0, 3, 3)]);

        const double Th_inv00 = (Th[ORDER(1, 1, 3, 3)]*Th[ORDER(2, 2, 3, 3)] -
                                 Th[ORDER(1, 2, 3, 3)]*Th[ORDER(2, 1, 3, 3)])/Th_det;
        const double Th_inv01 = (Th[ORDER(0, 2, 3, 3)]*Th[ORDER(2, 1, 3, 3)] -
                                 Th[ORDER(0, 1, 3, 3)]*Th[ORDER(2, 2, 3, 3)])/Th_det;
        const double Th_inv02 = (Th[ORDER(0, 1, 3, 3)]*Th[ORDER(1, 2, 3, 3)] -
                                 Th[ORDER(0, 2, 3, 3)]*Th[ORDER(1, 1, 3, 3)])/Th_det;

        const double Th_inv10 = (Th[ORDER(1, 2, 3, 3)]*Th[ORDER(2, 0, 3, 3)] -
                                 Th[ORDER(1, 0, 3, 3)]*Th[ORDER(2, 2, 3, 3)])/Th_det;
        const double Th_inv11 = (Th[ORDER(0, 0, 3, 3)]*Th[ORDER(2, 2, 3, 3)] -
                                 Th[ORDER(0, 2, 3, 3)]*Th[ORDER(2, 0, 3, 3)])/Th_det;
        const double Th_inv12 = (Th[ORDER(0, 2, 3, 3)]*Th[ORDER(1, 0, 3, 3)] -
                                 Th[ORDER(0, 0, 3, 3)]*Th[ORDER(1, 2, 3, 3)])/Th_det;

        const double Th_inv20 = (Th[ORDER(1, 0, 3, 3)]*Th[ORDER(2, 1, 3, 3)] -
                                 Th[ORDER(1, 1, 3, 3)]*Th[ORDER(2, 0, 3, 3)])/Th_det;
        const double Th_inv21 = (Th[ORDER(2, 0, 3, 3)]*Th[ORDER(0, 1, 3, 3)] -
                                 Th[ORDER(0, 0, 3, 3)]*Th[ORDER(2, 1, 3, 3)])/Th_det;
        const double Th_inv22 = (Th[ORDER(0, 0, 3, 3)]*Th[ORDER(1, 1, 3, 3)] -
                                 Th[ORDER(0, 1, 3, 3)]*Th[ORDER(1, 0, 3, 3)])/Th_det;

        /*
         * Perform the transformation using Th inverse. This is equivalent to
         * putting the reference values in long columns side-by-side and
         * multiplying by Th_inv.
         */
        for (i = 0; i < i_twentyone*quad_points; i++) {
                dxx_unmapped[i] = ref_dxx[i]*Th_inv00 + ref_dxy[i]*Th_inv10 +
                        ref_dyy[i]*Th_inv20;
                dxy_unmapped[i] = ref_dxx[i]*Th_inv01 + ref_dxy[i]*Th_inv11 +
                        ref_dyy[i]*Th_inv21;
                dyy_unmapped[i] = ref_dxx[i]*Th_inv02 + ref_dxy[i]*Th_inv12 +
                        ref_dyy[i]*Th_inv22;
        }

        /* perform the transformation using the C matrix. */
        DGEMM_WRAPPER(i_twentyone, quad_points, i_twentyone, C, dxx_unmapped,
                      dxx);
        DGEMM_WRAPPER(i_twentyone, quad_points, i_twentyone, C, dxy_unmapped,
                      dxy);
        DGEMM_WRAPPER(i_twentyone, quad_points, i_twentyone, C, dyy_unmapped,
                      dyy);
}
