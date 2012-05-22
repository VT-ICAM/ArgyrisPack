#include "mex.h"
#include "blas.h"

#include "ArgyrisGradients.c"

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
