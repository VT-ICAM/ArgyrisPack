#include <string.h>
#include <stdlib.h>
#include "ArgyrisHessians.c"

void multiplyByDiagonal(const LAPACKINDEX rows, const LAPACKINDEX cols,
                        double* restrict diagonal, double* restrict matrix);

void localLaplacian(double* restrict C, double* restrict B, double* restrict Th,
                    double* restrict rDxx, double* restrict rDxy,
                    double* restrict rDyy, double* restrict weights,
                    double* restrict laplacian,
                    LAPACKINDEX quadPoints, LAPACKINDEX rows) {

    // global values.
    LAPACKINDEX i;
    double dxx[rows*quadPoints];
    double dxy[rows*quadPoints];
    double dyy[rows*quadPoints];
    // stuff for DGEMM
    char cN = 'N';
    char cT = 'T';
    double dZero = 0.0;
    double jacobian = fabs(B[columnOrder(0,0,2)] * B[columnOrder(1,1,2)] -
                           B[columnOrder(0,1,2)] * B[columnOrder(1,0,2)]);

    evaluateArgyrisHessians(C, Th, rDxx, rDxy, rDyy, dxx, dxy, dyy,
                            quadPoints, rows);
    // dxx := dxx + dyy and dyy := dxx + dyy
    for (i = 0; i < rows*quadPoints; i++) {
        dxx[i] += dyy[i];
    }
    memcpy(dyy, dxx, sizeof(double)*(rows*quadPoints));

    // scale the first set of laplacian values by the weights and determinant.
    multiplyByDiagonal(rows, quadPoints, weights, dxx);
    // perform matrix multiplication.
    dgemm(&cN, &cT, &rows, &rows, &quadPoints, &jacobian, dxx, &rows,
          dyy, &rows, &dZero, laplacian, &rows);
}

void multiplyByDiagonal(const LAPACKINDEX rows, const LAPACKINDEX cols,
                        double* restrict diagonal, double* restrict matrix) {
    int i, j;
    for (i = 0; i < cols; i++) {
        for (j = 0; j < rows; j++) {
            matrix[columnOrder(j,i,rows)] *= diagonal[i];
        }
    }
}
