#if !defined(_WIN32)
#define dgemm dgemm_
#endif

#include "mex.h"
#include "blas.h"

#define columnOrder(row,col,entriesPerColumn) col * entriesPerColumn + row

void evaluateArgyrisGradients(double *C, double *B,
     double *referenceArgyrisGradientX, double *referenceArgyrisGradientY,
     double *argyrisGradientX,          double *argyrisGradientY,
     mwSignedIndex quadPoints, mwSignedIndex rows)
{
    double BdeterminantInverse;
    double Binverse00, Binverse01, Binverse10, Binverse11;
    double argyrisGradientXUnMapped[rows * quadPoints];
    double argyrisGradientYUnMapped[rows * quadPoints];
    // stuff for DGEMM
    char *chn = "N";
    double one = 1.0, zero = 0.0;

    BdeterminantInverse = 1/(B[columnOrder(0,0,2)] * B[columnOrder(1,1,2)] -
        B[columnOrder(0,1,2)] * B[columnOrder(1,0,2)]);

    Binverse00 = BdeterminantInverse * B[columnOrder(1,1,2)];
    Binverse01 = -BdeterminantInverse * B[columnOrder(0,1,2)];
    Binverse10 = -BdeterminantInverse * B[columnOrder(1,0,2)];
    Binverse11 = BdeterminantInverse * B[columnOrder(0,0,2)];

    // perform the transformation using B inverse. This is equivalent to
    // putting the reference values in long columns side-by-side and
    // multiplying by inverseB.
    for (int index = 0; index < rows * quadPoints; index++)
    {
        argyrisGradientXUnMapped[index] = Binverse00 * referenceArgyrisGradientX[index]
        + Binverse10 * referenceArgyrisGradientY[index];

        argyrisGradientYUnMapped[index] = Binverse01 * referenceArgyrisGradientX[index]
        + Binverse11 * referenceArgyrisGradientY[index];
    }

    // perform the transformation using the C matrix.
    dgemm(chn, chn, &rows, &quadPoints, &rows, &one, C, &rows,
          argyrisGradientXUnMapped, &rows, &zero, argyrisGradientX, &rows);
    dgemm(chn, chn, &rows, &quadPoints, &rows, &one, C, &rows,
          argyrisGradientYUnMapped, &rows, &zero, argyrisGradientY, &rows);
 }

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *B, *C;
    double *referenceArgyrisGradientX, *referenceArgyrisGradientY;
    double *argyrisGradientX, *argyrisGradientY;
    mwSignedIndex quadPoints, rows;

    // check input.
    if (nrhs != 4)
    {
        mexErrMsgIdAndTxt("ARGYRISKERNEL:ArgyrisGradients",
                          "Requires four arguements.");
    }
    C = mxGetPr(prhs[0]);
    B = mxGetPr(prhs[1]);
    referenceArgyrisGradientX = mxGetPr(prhs[2]);
    referenceArgyrisGradientY = mxGetPr(prhs[3]);
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
    if ((2 != mxGetN(prhs[1])) || (2 != mxGetM(prhs[1])))
    {
        mexErrMsgIdAndTxt("ARGYRISKERNEL:ArgyrisGradients",
                          "The B matrix must be the second arguement.");
    }
    if (21 != mxGetM(prhs[2]))
    {
        mexErrMsgIdAndTxt("ARGYRISKERNEL:ArgyrisGradients",
          "There should be 21 basis functions corresponding to 21 rows.");
    }
    if (21 != mxGetM(prhs[3]))
    {
        mexErrMsgIdAndTxt("ARGYRISKERNEL:ArgyrisGradients",
          "There should be 21 basis functions corresponding to 21 rows.");
    }

    // check output.
    plhs[0] = mxCreateDoubleMatrix(rows, quadPoints, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(rows, quadPoints, mxREAL);
    argyrisGradientX = mxGetPr(plhs[0]);
    argyrisGradientY = mxGetPr(plhs[1]);

    // generate output.
    evaluateArgyrisGradients(C, B, referenceArgyrisGradientX,
         referenceArgyrisGradientY, argyrisGradientX,
                             argyrisGradientY, quadPoints, rows);
}
