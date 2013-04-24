void ap_physical_values(double* restrict C, double* restrict ref_values,
                        LAPACKINDEX num_points, double* restrict values)
{
        LAPACKINDEX i_twentyone = 21;
        DGEMM_WRAPPER(i_twentyone, num_points, i_twentyone, C, ref_values,
                      values);
}
