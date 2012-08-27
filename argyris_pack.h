/* Note that LAPACKINDEX is not 'int' only for building MEX files. */
#ifndef LAPACKINDEX
#define LAPACKINDEX int
#endif

void ap_local_functions(double* restrict x, double* restrict y, LAPACKINDEX quad_points,
                        double* restrict ref_functions);

void ap_local_gradients(double* restrict x, double* restrict y, LAPACKINDEX quad_points,
                        double* restrict ref_dx, double* restrict ref_dy);

void ap_local_hessians(double* restrict x, double* restrict y, LAPACKINDEX quad_points,
                       double* restrict ref_dxx, double* restrict ref_dxy,
                       double* restrict ref_dyy);

void ap_global_maps(double* restrict x, double* restrict y, double* restrict C,
                    double* restrict B, double* restrict b, double* restrict Th);

void ap_global_functions(double* restrict C, double* restrict ref_values,
                         LAPACKINDEX quad_points, double* restrict values);

void ap_global_gradients(double* restrict C, double* restrict B,
                         double* restrict ref_dx, double* restrict ref_dy,
                         LAPACKINDEX quad_points,
                         double* restrict dx, double* restrict dy);

void ap_global_hessians(double* restrict C, double* restrict Th,
                        double* restrict ref_dxx, double* restrict ref_dxy,
                        double* restrict ref_dyy, LAPACKINDEX quad_points,
                        double* restrict dxx, double* restrict dxy,
                        double* restrict dyy);

void ap_matrix_mass(double* restrict C, double* restrict B,
                    double* restrict ref_functions, double* restrict weights,
                    LAPACKINDEX quad_points, double* restrict mass);

void ap_matrix_stiffness(double* restrict C, double* restrict B,
                         double* restrict ref_dx, double* restrict ref_dy,
                         double* restrict weights,
                         LAPACKINDEX quad_points, double* restrict stiffness);

void ap_matrix_biharmonic(double* restrict C, double* restrict B,
                          double* restrict Th,
                          double* restrict ref_dxx, double* restrict ref_dxy,
                          double* restrict ref_dyy, double* restrict weights,
                          LAPACKINDEX quad_points, double* restrict biharmonic);

#undef LAPACKINDEX
