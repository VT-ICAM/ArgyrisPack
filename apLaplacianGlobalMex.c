#include "mex.h"
#include "blas.h"

#include "laplacian.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *C, *B, *Th;
    double *rArgyrisXX, *rArgyrisXY, *rArgyrisYY;
    double* weights;
    double* laplacian;
    mwSignedIndex quadPoints, rows;

    // check input.
    if (nrhs != 7) {
        mexErrMsgIdAndTxt("ARGYRISPACK:ArgyrisLaplacian",
                          "Requires seven arguements.");
    }
    // check output.
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("ARGYRISPACK:ArgyrisLaplacian",
                          "Requires one output.");
    }

    C = mxGetPr(prhs[0]);
    B = mxGetPr(prhs[1]);
    Th = mxGetPr(prhs[2]);
    rArgyrisXX = mxGetPr(prhs[3]);
    rArgyrisXY = mxGetPr(prhs[4]);
    rArgyrisYY = mxGetPr(prhs[5]);
    weights = mxGetPr(prhs[6]);

    quadPoints = mxGetN(prhs[3]);
    rows = mxGetM(prhs[3]);

    if (quadPoints != mxGetN(prhs[4])) {
        mexErrMsgIdAndTxt("ARGYRISPACK:Laplacian",
                          "Mismatch in number of quadrature points.");
    }
    if (quadPoints != mxGetN(prhs[5])) {
        mexErrMsgIdAndTxt("ARGYRISPACK:Laplacian",
                          "Mismatch in number of quadrature points.");
    }
    if ((21 != mxGetN(prhs[0])) || (21 != mxGetM(prhs[0]))) {
        mexErrMsgIdAndTxt("ARGYRISPACK:Laplacian",
                          "The C matrix must be the first arguement.");
    }
    if ((2 != mxGetN(prhs[1])) || (2 != mxGetM(prhs[1]))) {
        mexErrMsgIdAndTxt("ARGYRISPACK:Laplacian",
                          "The B matrix must be the second arguement.");
    }
    if ((3 != mxGetN(prhs[2])) || (3 != mxGetM(prhs[2]))) {
        mexErrMsgIdAndTxt("ARGYRISPACK:Laplacian",
                          "The Th matrix must be the third arguement.");
    }
    if (quadPoints != mxGetM(prhs[6]) * mxGetN(prhs[6])) {
        mexErrMsgIdAndTxt("ARGYRISPACK:Laplacian",
                          "Mismatch in number of quadrature weights.");
    }

    // output.
    plhs[0]   = mxCreateDoubleMatrix(rows, rows, mxREAL);
    laplacian = mxGetPr(plhs[0]);

    // generate output.
    localLaplacian(C, B, Th, rArgyrisXX, rArgyrisXY, rArgyrisYY,
                   weights, laplacian, quadPoints, rows);
}
