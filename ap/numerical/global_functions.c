void ap_global_functions(double* restrict C, double* restrict ref_functions,
                         LAPACKINDEX num_points, double* restrict functions)
{
        LAPACKINDEX i_twentyone = 21;
        DGEMM_WRAPPER(i_twentyone, num_points, i_twentyone, C, ref_functions,
                      functions);
}
