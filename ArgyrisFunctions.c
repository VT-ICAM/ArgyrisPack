void evaluateArgyrisFunctions(double *C, double *rArgyrisFunctions,
                              double *argyrisFunctions,
                              LAPACKINDEX quadPoints, LAPACKINDEX rows)
{
    // stuff for DGEMM
    char *chn = "N";
    double one = 1.0, zero = 0.0;

    // perform the transformation using the C matrix.
    dgemm(chn, chn, &rows, &quadPoints, &rows, &one, C, &rows,
          rArgyrisFunctions, &rows, &zero, argyrisFunctions, &rows);
}
