#define columnOrder(row,col,entriesPerColumn) (col) * (entriesPerColumn) + (row)

void evaluateArgyrisGradients(double * restrict C, double * restrict B,
     double * restrict rArgyrisGradientX, double * restrict rArgyrisGradientY,
     double * restrict argyrisGradientX, double * restrict argyrisGradientY,
     LAPACKINDEX quadPoints, LAPACKINDEX rows)
{
    double BdeterminantInverse;
    double Binverse00, Binverse01, Binverse10, Binverse11;
    double argyrisGradientXUnMapped[rows * quadPoints];
    double argyrisGradientYUnMapped[rows * quadPoints];
    // stuff for DGEMM
    char *chn = "N";
    double one = 1.0, zero = 0.0;
    int index;

    BdeterminantInverse = 1/(B[columnOrder(0,0,2)] * B[columnOrder(1,1,2)] -
        B[columnOrder(0,1,2)] * B[columnOrder(1,0,2)]);

    Binverse00 = BdeterminantInverse * B[columnOrder(1,1,2)];
    Binverse01 = -BdeterminantInverse * B[columnOrder(0,1,2)];
    Binverse10 = -BdeterminantInverse * B[columnOrder(1,0,2)];
    Binverse11 = BdeterminantInverse * B[columnOrder(0,0,2)];

    // perform the transformation using B inverse. This is equivalent to
    // putting the reference values in long columns side-by-side and
    // multiplying by inverseB.
    for (index = 0; index < rows * quadPoints; index++) {
        argyrisGradientXUnMapped[index] = Binverse00 * rArgyrisGradientX[index]
        + Binverse10 * rArgyrisGradientY[index];

        argyrisGradientYUnMapped[index] = Binverse01 * rArgyrisGradientX[index]
        + Binverse11 * rArgyrisGradientY[index];
    }

    // perform the transformation using the C matrix.
    dgemm(chn, chn, &rows, &quadPoints, &rows, &one, C, &rows,
          argyrisGradientXUnMapped, &rows, &zero, argyrisGradientX, &rows);
    dgemm(chn, chn, &rows, &quadPoints, &rows, &one, C, &rows,
          argyrisGradientYUnMapped, &rows, &zero, argyrisGradientY, &rows);
}
