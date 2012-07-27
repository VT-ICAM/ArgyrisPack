#include "mex.h"
#include "blas.h"

#include "gradients.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *B, *C;
    double *refdx, *refdy;
    double *dx, *dy;
    mwSignedIndex quadPoints, rows;

    // check input.
    if (nrhs != 4) {
        mexErrMsgIdAndTxt("ARGYRISPACK:apGradientsMex",
                          "Requires four arguements.");
    }
    C = mxGetPr(prhs[0]);
    B = mxGetPr(prhs[1]);
    refdx = mxGetPr(prhs[2]);
    refdy = mxGetPr(prhs[3]);
    quadPoints = mxGetN(prhs[2]);
    rows = mxGetM(prhs[0]);

    if (21 != mxGetN(prhs[0]) || 21 != mxGetM(prhs[0])) {
        mexErrMsgIdAndTxt("ARGYRISPACK:apGradientsMex",
                          "The C matrix must be the first arguement.");
    }
    if ((2 != mxGetN(prhs[1])) || (2 != mxGetM(prhs[1]))) {
        mexErrMsgIdAndTxt("ARGYRISPACK:apGradientsMex",
                          "The B matrix must be the second arguement.");
    }
    if (quadPoints != mxGetN(prhs[2]) || quadPoints != mxGetN(prhs[3])) {
        mexErrMsgIdAndTxt("ARGYRISPACK:apGradientsMex",
                          "Mismatch in number of quadrature points.");
    }
    if (21 != mxGetM(prhs[2]) || 21 != mxGetM(prhs[3])) {
        mexErrMsgIdAndTxt("ARGYRISPACK:apGradientsMex",
          "There should be 21 basis functions corresponding to 21 rows.");
    }

    // check output.
    if (nlhs != 2) {
        mexErrMsgIdAndTxt("ARGYRISPACK:apGradientsMex",
                          "Requires two outputs.");
    }
    plhs[0] = mxCreateDoubleMatrix(rows, quadPoints, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(rows, quadPoints, mxREAL);
    dx = mxGetPr(plhs[0]);
    dy = mxGetPr(plhs[1]);

    // generate output.
    ap_gradients(C, B, refdx, refdy, quadPoints, rows, dx, dy);
}
