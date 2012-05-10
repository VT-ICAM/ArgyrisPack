#if !defined(_WIN32)
#define dgemm dgemm_
#endif

#include "mex.h"
#include "blas.h"

#define columnOrder(row,col,entriesPerColumn) col * entriesPerColumn + row

void evaluateArgyrisHessians(double *C, double *Th,
     double *rArgyrisXX, double *rArgyrisXY,
     double *rArgyrisYY, double *argyrisXX,
     double *argyrisXY,          double *argyrisYY,
     mwSignedIndex quadPoints, mwSignedIndex rows)
{
    double ThDeterminant;
    double ThInverse00, ThInverse01, ThInverse02;
    double ThInverse10, ThInverse11, ThInverse12;
    double ThInverse20, ThInverse21, ThInverse22;
    double argyrisXXUnMapped[rows * quadPoints];
    double argyrisXYUnMapped[rows * quadPoints];
    double argyrisYYUnMapped[rows * quadPoints];

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
    for (int index = 0; index < rows * quadPoints; index++)
    {
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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *Th, *C;
    double *rArgyrisXX, *rArgyrisXY, *rArgyrisYY;
    double *argyrisXX, *argyrisXY, *argyrisYY;
    mwSignedIndex quadPoints, rows;

    // check input.
    if (nrhs != 5)
    {
        mexErrMsgIdAndTxt("ARGYRISKERNEL:ArgyrisHessians",
                          "Requires five arguements.");
    }
    // check output.
    if (nlhs != 3)
    {
        mexErrMsgIdAndTxt("ARGYRISKERNEL:ArgyrisHessians",
                          "Requires three outputs.");
    }

    C = mxGetPr(prhs[0]);
    Th = mxGetPr(prhs[1]);
    rArgyrisXX = mxGetPr(prhs[2]);
    rArgyrisXY = mxGetPr(prhs[3]);
    rArgyrisYY = mxGetPr(prhs[4]);

    quadPoints = mxGetN(prhs[2]);
    rows = mxGetM(prhs[0]);

    if (quadPoints != mxGetN(prhs[3]))
    {
        mexErrMsgIdAndTxt("ARGYRISKERNEL:ArgyrisGradients",
                          "Mismatch in number of quadrature points.");
    }
    if ((21 != mxGetN(prhs[0])) || (21 != mxGetM(prhs[0])))
    {
        mexErrMsgIdAndTxt("ARGYRISKERNEL:ArgyrisGradients",
                          "The C matrix must be the first arguement.");
    }
    if ((3 != mxGetN(prhs[1])) || (3 != mxGetM(prhs[1])))
    {
        mexErrMsgIdAndTxt("ARGYRISKERNEL:ArgyrisGradients",
                          "The Th matrix must be the second arguement.");
    }
    for (int i = 2; i < 5; ++i)
    {
        if (21 != mxGetM(prhs[2]))
        {
            mexErrMsgIdAndTxt("ARGYRISKERNEL:ArgyrisGradients",
              "There should be 21 basis functions corresponding to 21 rows.");
        }
    }

    // check output.
    plhs[0] = mxCreateDoubleMatrix(rows, quadPoints, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(rows, quadPoints, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(rows, quadPoints, mxREAL);
    argyrisXX = mxGetPr(plhs[0]);
    argyrisXY = mxGetPr(plhs[1]);
    argyrisYY = mxGetPr(plhs[2]);

    // generate output.
    evaluateArgyrisHessians(C, Th, rArgyrisXX, rArgyrisXY, rArgyrisYY,
                            argyrisXX, argyrisXY, argyrisYY, quadPoints, rows);
}
