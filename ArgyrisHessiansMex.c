#include "mex.h"
#include "blas.h"

#include "ArgyrisHessians.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *Th, *C;
    double *rArgyrisXX, *rArgyrisXY, *rArgyrisYY;
    double *argyrisXX, *argyrisXY, *argyrisYY;
    mwSignedIndex quadPoints, rows;
    int i;

    // check input.
    if (nrhs != 5) {
        mexErrMsgIdAndTxt("ARGYRISKERNEL:ArgyrisHessians",
                          "Requires five arguements.");
    }
    // check output.
    if (nlhs != 3) {
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

    if (quadPoints != mxGetN(prhs[3])) {
        mexErrMsgIdAndTxt("ARGYRISKERNEL:ArgyrisGradients",
                          "Mismatch in number of quadrature points.");
    }
    if ((21 != mxGetN(prhs[0])) || (21 != mxGetM(prhs[0]))) {
        mexErrMsgIdAndTxt("ARGYRISKERNEL:ArgyrisGradients",
                          "The C matrix must be the first arguement.");
    }
    if ((3 != mxGetN(prhs[1])) || (3 != mxGetM(prhs[1]))) {
        mexErrMsgIdAndTxt("ARGYRISKERNEL:ArgyrisGradients",
                          "The Th matrix must be the second arguement.");
    }
    for (i = 2; i < 5; ++i) {
        if (21 != mxGetM(prhs[2])) {
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
