void ap_functions(double* restrict C, double* restrict ref_values,
                  LAPACKINDEX quad_points, double* restrict values)
{
        int i_twentyone = 21;
        DGEMM_WRAPPER(i_twentyone, quad_points, i_twentyone, C, ref_values,
                      values);
}
