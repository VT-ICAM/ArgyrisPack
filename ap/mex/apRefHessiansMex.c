#include "mex.h"
#include "blas.h"

#include "order_logic.h"
#include "ref_hessians.c"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
        double *x, *y;
        double *refdxx, *refdxy, *refdyy;
        mwSignedIndex numPoints;
        int i;

        /* check input. */
        if (nrhs != 2) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apRefHessiansMex",
                                  "Requires two arguments.");
        }
        for (i = 0; i < 2; i++) {
                if (mxGetN(prhs[i]) != 1 && mxGetM(prhs[i]) != 1) {
                        mexErrMsgIdAndTxt("ARGYRISPACK:apRefHessiansMex",
                                          "Inputs must be vectors.");
                }
        }
        if (mxGetM(prhs[0])*mxGetN(prhs[0]) != (mxGetM(prhs[1])*mxGetN(prhs[1]))) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apRefHessiansMex",
                                  "Inputs must have same size.");
        }

        x = mxGetPr(prhs[0]);
        y = mxGetPr(prhs[1]);
        numPoints = mxGetN(prhs[1])*mxGetM(prhs[1]);

        /* check output. */
        if (3 != nlhs) {
                mexErrMsgIdAndTxt("ARGYRISPACK:apRefHessiansMex",
                                  "Requires exactly three outputs.");
        }

        plhs[0] = mxCreateDoubleMatrix(21, numPoints, mxREAL);
        plhs[1] = mxCreateDoubleMatrix(21, numPoints, mxREAL);
        plhs[2] = mxCreateDoubleMatrix(21, numPoints, mxREAL);
        refdxx = mxGetPr(plhs[0]);
        refdxy = mxGetPr(plhs[1]);
        refdyy = mxGetPr(plhs[2]);

        ap_ref_hessians(x, y, numPoints, refdxx, refdxy, refdyy);
}
