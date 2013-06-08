#include "mex.h"
#include "blas.h"
#include "order_logic.h"
#include "physical_hessians.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
        double *B, *C;
        double *refdxx, *refdxy, *refdyy;
        double *dxx, *dxy, *dyy;
        mwSignedIndex quadPoints;
        int i;

        /* check input. */
        if (nrhs != 5) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apPhysicalHessiansMex",
                                  "Requires five arguments.");
        }
        if ((21 != mxGetN(prhs[0])) || (21 != mxGetM(prhs[0]))) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apPhysicalHessiansMex",
                                  "The C matrix must be the first argument.");
        }
        if ((2 != mxGetN(prhs[1])) || (2 != mxGetM(prhs[1]))) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apPhysicalHessiansMex",
                                  "The B matrix must be the second argument.");
        }
        for (i = 2; i < 5; i++) {
                if (21 != mxGetM(prhs[i])) {
                        mexErrMsgIdAndTxt("ARGYRISPACK:apPhysicalHessiansMex",
                                          "There should be 21 basis functions "
                                          "corresponding to 21 rows.");
                }
                if (mxGetN(prhs[2]) != mxGetN(prhs[i])) {
                        mexErrMsgIdAndTxt("ARGYRISPACK:apPhysicalHessiansMex",
                                          "Mismatch in number of quadrature "
                                          "points.");
                }
        }

        C = mxGetPr(prhs[0]);
        B = mxGetPr(prhs[1]);
        refdxx = mxGetPr(prhs[2]);
        refdxy = mxGetPr(prhs[3]);
        refdyy = mxGetPr(prhs[4]);
        quadPoints = mxGetN(prhs[2]);

        /* check output. */
        if (nlhs != 3) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apPhysicalHessiansMex",
                                  "Requires three outputs.");
        }

        plhs[0] = mxCreateDoubleMatrix(21, quadPoints, mxREAL);
        plhs[1] = mxCreateDoubleMatrix(21, quadPoints, mxREAL);
        plhs[2] = mxCreateDoubleMatrix(21, quadPoints, mxREAL);
        dxx = mxGetPr(plhs[0]);
        dxy = mxGetPr(plhs[1]);
        dyy = mxGetPr(plhs[2]);

        ap_physical_hessians(C, B, refdxx, refdxy, refdyy, quadPoints, dxx, dxy,
                             dyy);
}
