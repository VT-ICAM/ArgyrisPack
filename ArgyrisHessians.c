#define columnOrder(row,col,entriesPerColumn) (col) * (entriesPerColumn) + (row) 
void evaluateArgyrisHessians(double *C, double *Th,
     double *rArgyrisXX, double *rArgyrisXY,
     double *rArgyrisYY, double *argyrisXX,
     double *argyrisXY,  double *argyrisYY,
     LAPACKINDEX quadPoints, LAPACKINDEX rows)
{
    double ThDeterminant;
    double ThInverse00, ThInverse01, ThInverse02;
    double ThInverse10, ThInverse11, ThInverse12;
    double ThInverse20, ThInverse21, ThInverse22;
    double argyrisXXUnMapped[rows * quadPoints];
    double argyrisXYUnMapped[rows * quadPoints];
    double argyrisYYUnMapped[rows * quadPoints];
    int index;

    // stuff for DGEMM
    char *chn = "N";
    double one = 1.0, zero = 0.0;

    // Calculate ThInverse entries.
    ThDeterminant =
    Th[columnOrder(0, 0, 3)] * (Th[columnOrder(1, 1, 3)]*Th[columnOrder(2, 2, 3)] -
                                Th[columnOrder(1, 2, 3)]*Th[columnOrder(2, 1, 3)]) +
    Th[columnOrder(0, 1, 3)] * (Th[columnOrder(1, 2, 3)]*Th[columnOrder(2, 0, 3)] -
                                Th[columnOrder(2, 2, 3)]*Th[columnOrder(1, 0, 3)]) +
    Th[columnOrder(0, 2, 3)] * (Th[columnOrder(1, 0, 3)]*Th[columnOrder(2, 1, 3)] -
                                Th[columnOrder(1, 1, 3)]*Th[columnOrder(2, 0, 3)]);

    ThInverse00 = (Th[columnOrder(1, 1, 3)]*Th[columnOrder(2, 2, 3)] -
        Th[columnOrder(1, 2, 3)]*Th[columnOrder(2, 1, 3)])/ThDeterminant;
    ThInverse01 = (Th[columnOrder(0, 2, 3)]*Th[columnOrder(2, 1, 3)] -
        Th[columnOrder(0, 1, 3)]*Th[columnOrder(2, 2, 3)])/ThDeterminant;
    ThInverse02 = (Th[columnOrder(0, 1, 3)]*Th[columnOrder(1, 2, 3)] -
        Th[columnOrder(0, 2, 3)]*Th[columnOrder(1, 1, 3)])/ThDeterminant;

    ThInverse10 = (Th[columnOrder(1, 2, 3)]*Th[columnOrder(2, 0, 3)] -
        Th[columnOrder(1, 0, 3)]*Th[columnOrder(2, 2, 3)])/ThDeterminant;
    ThInverse11 = (Th[columnOrder(0, 0, 3)]*Th[columnOrder(2, 2, 3)] -
        Th[columnOrder(0, 2, 3)]*Th[columnOrder(2, 0, 3)])/ThDeterminant;
    ThInverse12 = (Th[columnOrder(0, 2, 3)]*Th[columnOrder(1, 0, 3)] -
        Th[columnOrder(0, 0, 3)]*Th[columnOrder(1, 2, 3)])/ThDeterminant;

    ThInverse20 = (Th[columnOrder(1, 0, 3)]*Th[columnOrder(2, 1, 3)] -
        Th[columnOrder(1, 1, 3)]*Th[columnOrder(2, 0, 3)])/ThDeterminant;
    ThInverse21 = (Th[columnOrder(2, 0, 3)]*Th[columnOrder(0, 1, 3)] -
        Th[columnOrder(0, 0, 3)]*Th[columnOrder(2, 1, 3)])/ThDeterminant;
    ThInverse22 = (Th[columnOrder(0, 0, 3)]*Th[columnOrder(1, 1, 3)] -
        Th[columnOrder(0, 1, 3)]*Th[columnOrder(1, 0, 3)])/ThDeterminant;

    // It does not matter what order we examine elements, or even if the
    // arrays are in column or row order. Therefore access everything as if it
    // was flat.
    for (index = 0; index < rows * quadPoints; index++) {
        argyrisXXUnMapped[index] = rArgyrisXX[index] * ThInverse00 +
            rArgyrisXY[index] * ThInverse10 + rArgyrisYY[index] * ThInverse20;

        argyrisXYUnMapped[index] = rArgyrisXX[index] * ThInverse01 +
            rArgyrisXY[index] * ThInverse11 + rArgyrisYY[index] * ThInverse21;

        argyrisYYUnMapped[index] = rArgyrisXX[index] * ThInverse02 +
            rArgyrisXY[index] * ThInverse12 + rArgyrisYY[index] * ThInverse22;
    }

    // perform the transformation using the C matrix.
    dgemm(chn, chn, &rows, &quadPoints, &rows, &one, C, &rows,
          argyrisXXUnMapped, &rows, &zero, argyrisXX, &rows);
    dgemm(chn, chn, &rows, &quadPoints, &rows, &one, C, &rows,
          argyrisXYUnMapped, &rows, &zero, argyrisXY, &rows);
    dgemm(chn, chn, &rows, &quadPoints, &rows, &one, C, &rows,
          argyrisYYUnMapped, &rows, &zero, argyrisYY, &rows);
}
