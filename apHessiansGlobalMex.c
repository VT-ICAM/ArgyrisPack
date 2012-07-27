#include "mex.h"
#include "blas.h"

#include "ArgyrisHessians.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *Th, *C;
    double *refdxx, *refdxy, *refdyy;
    double *dxx, *dxy, *dyy;
    mwSignedIndex quadPoints, rows;
    int i;

    // check input.
    if (nrhs != 5) {
        mexErrMsgIdAndTxt("ARGYRISPACK:ArgyrisHessiansMex",
                          "Requires five arguements.");
    }

    C = mxGetPr(prhs[0]);
    Th = mxGetPr(prhs[1]);
    refdxx = mxGetPr(prhs[2]);
    refdxy = mxGetPr(prhs[3]);
    refdyy = mxGetPr(prhs[4]);

    quadPoints = mxGetN(prhs[2]);
    rows = mxGetM(prhs[0]);

    if ((21 != mxGetN(prhs[0])) || (21 != mxGetM(prhs[0]))) {
        mexErrMsgIdAndTxt("ARGYRISPACK:ArgyrisHessiansMex",
                          "The C matrix must be the first arguement.");
    }
    if ((3 != mxGetN(prhs[1])) || (3 != mxGetM(prhs[1]))) {
        mexErrMsgIdAndTxt("ARGYRISPACK:ArgyrisHessiansMex",
                          "The Th matrix must be the second arguement.");
    }
    for (i = 2; i < 5; ++i) {
        if (21 != mxGetM(prhs[i])) {
            mexErrMsgIdAndTxt("ARGYRISPACK:ArgyrisHessiansMex",
              "There should be 21 basis functions corresponding to 21 rows.");
        }
        if (quadPoints != mxGetN(prhs[i])) {
            mexErrMsgIdAndTxt("ARGYRISPACK:ArgyrisHessiansMex",
                              "Mismatch in number of quadrature points.");
        }
    }

    // check output.
    if (nlhs != 3) {
        mexErrMsgIdAndTxt("ARGYRISPACK:ArgyrisHessiansMex",
                          "Requires three outputs.");
    }
    plhs[0] = mxCreateDoubleMatrix(rows, quadPoints, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(rows, quadPoints, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(rows, quadPoints, mxREAL);
    dxx = mxGetPr(plhs[0]);
    dxy = mxGetPr(plhs[1]);
    dyy = mxGetPr(plhs[2]);

    // generate output.
    evaluateArgyrisHessians(C, Th, refdxx, refdxy, refdyy, dxx, dxy, dyy,
                            quadPoints, rows);
}
